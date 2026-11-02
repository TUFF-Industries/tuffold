//! Proteome Discoverer 1.x `.msf` (checked on 1.4) → one batch per search
//! node. The same rules as 2.x/3.x, on 1.x's own tables
//! (`docs/pd-format.md`):
//! - Each search node with peptides (`Peptides`, `Peptides_decoy`) is one
//!   engine, named as PD names the node (`Mascot`, `Sequest HT`; with its
//!   node number when a workflow ran two of one kind). Ranks are
//!   `SearchEngineRank`.
//! - Scores are the engine's own (`PeptideScores` named by
//!   `ProcessingNodeScores`). Percolator's q-values and PEPs
//!   (`CustomDataPeptides`) and PD's confidence levels are not imported.
//! - Residue modifications sit at **0-based** positions
//!   (`PeptidesAminoAcidModifications`); terminal ones are listed apart
//!   (`PeptidesTerminalModifications`), position types 1 and 3 N-terminal,
//!   2 and 4 C-terminal. Fixed modifications are listed per PSM too.
//! - 1.x stores no theoretical mass, so every PSM's mass is checked against
//!   its spectrum's precursor MH+ (`SpectrumHeaders.Mass`), allowing
//!   isotope peaks −1..+3: one off by more than [`MASS_TOLERANCE`] is left
//!   out and counted, and the file is refused when over
//!   [`MAX_MASS_MISSES`] are.
//! - Decoy peptides map to their target's protein id; they import as
//!   `DECOY_<accession>`, flagged explicitly. 1.x stores only proteins some
//!   target matched and leaves other PSMs unlinked (a third of the decoys
//!   in the sample): those are linked by sequence to the stored proteins
//!   containing their peptide (decoys: reversed), and a decoy still without
//!   one is filed under a placeholder `DECOY_unstored:<peptide>`, because a
//!   PSM without proteins would count as a target. A protein's accession is the
//!   first word of its FASTA title line (`ProteinAnnotations`); proteins
//!   sharing one are named by the whole line.
//! - Residues PD writes as `X` (with a "MappingL"-style substitution) or
//!   other non-standard letters can't be mass-checked; those PSMs come in
//!   unchecked and are counted, as in 2.x.

use std::collections::{HashMap, HashSet};

use rusqlite::Connection;
use tuffold_core::chemistry::{PROTON, peptide_mass};
use tuffold_core::{ModSite, Modification, PeptideForm, Protein, Psm, Score, Spectrum};

use super::{DECOY_PREFIX, Draft, Notes, PdError, PdImport, batches, spectrum_key};
use crate::xml::run_name;

/// Daltons from the nearest isotope peak of the precursor. Wider than any
/// Orbitrap search tolerance, far below a misread modification.
const MASS_TOLERANCE: f64 = 0.05;
/// Share of PSMs allowed to miss their precursor before the file is refused.
const MAX_MASS_MISSES: f64 = 0.05;
const ISOTOPE_SPACING: f64 = 1.003_354_835;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Side {
    Target,
    Decoy,
}

impl Side {
    fn suffix(self) -> &'static str {
        match self {
            Self::Target => "",
            Self::Decoy => "_decoy",
        }
    }
}

/// (node number, peptide id).
type PeptideKey = (i64, i64);

struct SpectrumRow {
    file: i64,
    first: Option<i64>,
    last: Option<i64>,
    rt_minutes: Option<f64>,
    charge: Option<i64>,
    /// Precursor MH+.
    mh: Option<f64>,
    mz: Option<f64>,
}

struct ModificationType {
    name: String,
    delta: f64,
    position_type: i64,
    unimod: Option<i64>,
}

pub(super) fn import(conn: &Connection, version: Option<&str>) -> Result<PdImport, PdError> {
    let runs = runs(conn)?;
    let spectra = spectra(conn)?;
    let engines = engines(conn)?;
    let types = modification_types(conn)?;
    let score_names = score_names(conn)?;
    let mut notes = Notes::default();
    let mut drafts = Vec::new();
    let mut checked = 0usize;
    let mut used_proteins: HashSet<i64> = HashSet::new();
    let mut stored = StoredSequences::read(conn)?;
    let mut links = Links::default();
    for side in [Side::Target, Side::Decoy] {
        let scores = scores(conn, side, &score_names)?;
        let modifications = modifications(conn, side, &types)?;
        let proteins = protein_links(conn, side)?;
        let suffix = side.suffix();
        let mut statement = conn.prepare(&format!(
            "SELECT ProcessingNodeNumber, PeptideID, SpectrumID, SearchEngineRank, Sequence \
             FROM Peptides{suffix}"
        ))?;
        let mut rows = statement.query([])?;
        while let Some(row) = rows.next()? {
            let key: PeptideKey = (row.get(0)?, row.get(1)?);
            let spectrum_id: i64 = row.get(2)?;
            let rank: Option<i64> = row.get(3)?;
            let sequence = row
                .get::<_, Option<String>>(4)?
                .unwrap_or_default()
                .to_ascii_uppercase();
            let Some(engine) = engines.get(&key.0) else {
                return Err(PdError::Layout(format!(
                    "peptides from processing node {}, which ProcessingNodes doesn't list",
                    key.0
                )));
            };
            let Some(spectrum) = spectra.get(&spectrum_id) else {
                notes.without_scan += 1;
                continue;
            };
            let Some(run) = runs.get(&spectrum.file) else {
                notes.without_file += 1;
                continue;
            };
            let Some(charge) = spectrum
                .charge
                .and_then(|charge| i16::try_from(charge).ok())
                .filter(|&charge| charge > 0)
            else {
                notes.without_charge += 1;
                continue;
            };
            if sequence.is_empty() || !sequence.bytes().all(|byte| byte.is_ascii_uppercase()) {
                return Err(PdError::Layout(format!(
                    "peptide {key:?} has sequence {sequence:?}"
                )));
            }
            let peptide = peptide_form(&sequence, modifications.get(&key), &types)?;
            let deltas: f64 = peptide.modifications.iter().map(|m| m.mass_delta).sum();
            let calculated_mh = peptide_mass(&sequence).map(|mass| mass + deltas + PROTON);
            match calculated_mh {
                Some(mh) => {
                    checked += 1;
                    if !fits_precursor(mh, spectrum.mh) {
                        notes.mass_misses += 1;
                        continue;
                    }
                }
                None => notes.unchecked_mass += 1,
            }
            let Some(spectrum_key) = spectrum_key(spectrum.first, spectrum.last) else {
                notes.without_scan += 1;
                continue;
            };
            let protein_ids = match proteins.get(&key) {
                Some(ids) => ids.iter().map(|&id| (side == Side::Decoy, id)).collect(),
                None => links.find(side, &peptide.sequence, &mut stored),
            };
            used_proteins.extend(protein_ids.iter().map(|&(_, id)| id).filter(|&id| id > 0));
            let z = f64::from(charge);
            let experimental_mz = spectrum.mz.filter(|mz| mz.is_finite() && *mz > 0.0);
            drafts.push(Draft {
                engine: engine.clone(),
                run: run.clone(),
                key: spectrum_key.clone(),
                spectrum: Spectrum {
                    ms_run: run.clone(),
                    key: spectrum_key.clone(),
                    source_id: None,
                    precursor_mz: experimental_mz,
                    precursor_charge: Some(charge),
                    rt_seconds: spectrum
                        .rt_minutes
                        .filter(|rt| rt.is_finite() && *rt >= 0.0)
                        .map(|rt| rt * 60.0),
                },
                proteins: protein_ids,
                psm: Psm {
                    ms_run: run.clone(),
                    spectrum: spectrum_key,
                    rank: rank
                        .and_then(|rank| u32::try_from(rank).ok())
                        .filter(|&rank| rank > 0)
                        .unwrap_or(1),
                    charge,
                    calculated_mz: calculated_mh.map(|mh| (mh - PROTON) / z + PROTON),
                    experimental_mz,
                    pass_threshold: None,
                    scores: scores.get(&key).cloned().unwrap_or_default(),
                    evidence: Vec::new(),
                    peptide,
                },
            });
        }
    }
    if notes.mass_misses as f64 > MAX_MASS_MISSES * checked.max(1) as f64 {
        return Err(PdError::Layout(format!(
            "{} of {checked} PSMs miss their precursor mass by more than {MASS_TOLERANCE} Da, \
             so the modifications were probably misread",
            notes.mass_misses
        )));
    }
    if drafts.is_empty() {
        return Err(PdError::Layout("no PSMs".into()));
    }
    let mut proteins = proteins(conn, &used_proteins)?;
    for (sequence, id) in &links.placeholders {
        proteins.insert(
            (true, *id),
            Protein {
                accession: format!("{DECOY_PREFIX}unstored:{sequence}"),
                description: Some(
                    "decoy match whose protein Proteome Discoverer did not store".into(),
                ),
                sequence: None,
                is_decoy: true,
                decoy_explicit: true,
            },
        );
    }
    let mut lines = notes.lines();
    lines.extend(links.notes());
    Ok(PdImport {
        batches: batches(drafts, &proteins, version),
        notes: lines,
    })
}

/// Stored protein sequences (I as L), for PSMs PD linked to no protein.
/// 1.x stores only proteins that some target matched, and links a decoy
/// PSM only when its target protein is stored.
struct StoredSequences {
    proteins: Vec<(i64, Vec<u8>)>,
    /// First [`SEED`] residues → (protein index, position).
    seeds: HashMap<[u8; SEED], Vec<(usize, usize)>>,
}

const SEED: usize = 5;

impl StoredSequences {
    fn read(conn: &Connection) -> Result<Self, PdError> {
        let mut statement = conn.prepare("SELECT ProteinID, Sequence FROM Proteins")?;
        let proteins: Vec<(i64, Vec<u8>)> = statement
            .query_map([], |row| {
                let sequence: Option<String> = row.get(1)?;
                Ok((
                    row.get::<_, i64>(0)?,
                    leucine(sequence.unwrap_or_default().as_bytes()),
                ))
            })?
            .collect::<rusqlite::Result<_>>()?;
        Ok(Self {
            proteins,
            seeds: HashMap::new(),
        })
    }

    /// Ids of the stored proteins containing `peptide` (I and L alike).
    fn containing(&mut self, peptide: &[u8]) -> Vec<i64> {
        if peptide.is_empty() {
            return Vec::new(); // a damaged row; `windows(0)` would panic
        }
        if peptide.len() < SEED {
            return self
                .proteins
                .iter()
                .filter(|(_, sequence)| sequence.windows(peptide.len()).any(|w| w == peptide))
                .map(|(id, _)| *id)
                .collect();
        }
        if self.seeds.is_empty() {
            for (index, (_, sequence)) in self.proteins.iter().enumerate() {
                for (position, window) in sequence.windows(SEED).enumerate() {
                    let seed: [u8; SEED] = window.try_into().expect("seed length");
                    self.seeds.entry(seed).or_default().push((index, position));
                }
            }
        }
        let seed: [u8; SEED] = peptide[..SEED].try_into().expect("seed length");
        let mut found: Vec<i64> = self
            .seeds
            .get(&seed)
            .into_iter()
            .flatten()
            .filter(|&&(index, position)| {
                self.proteins[index]
                    .1
                    .get(position..position + peptide.len())
                    == Some(peptide)
            })
            .map(|&(index, _)| self.proteins[index].0)
            .collect();
        found.sort_unstable();
        found.dedup();
        found
    }
}

fn leucine(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|&residue| match residue.to_ascii_uppercase() {
            b'I' => b'L',
            other => other,
        })
        .collect()
}

/// Proteins found for PSMs PD left unlinked, and what was done.
#[derive(Default)]
struct Links {
    by_sequence: [usize; 2],
    /// Decoy peptide → its placeholder protein's (negative) id.
    placeholders: HashMap<String, i64>,
    placeholder_psms: usize,
    unlinked_targets: usize,
}

impl Links {
    /// A target PSM's proteins are the stored ones containing its peptide;
    /// a decoy's, the decoys of those containing it reversed (whole, or
    /// but its last residue). A decoy with none gets a placeholder decoy
    /// protein, so it still counts as one; a target with none stays
    /// without protein (it counts as a target either way).
    fn find(
        &mut self,
        side: Side,
        sequence: &str,
        stored: &mut StoredSequences,
    ) -> Vec<(bool, i64)> {
        let peptide = leucine(sequence.as_bytes());
        let ids = match side {
            Side::Target => stored.containing(&peptide),
            Side::Decoy => {
                let reversed: Vec<u8> = peptide.iter().rev().copied().collect();
                let mut ids = stored.containing(&reversed);
                if peptide.len() > SEED {
                    let but_last: Vec<u8> =
                        peptide[..peptide.len() - 1].iter().rev().copied().collect();
                    ids.extend(stored.containing(&but_last));
                }
                ids.sort_unstable();
                ids.dedup();
                ids
            }
        };
        let decoy = side == Side::Decoy;
        if !ids.is_empty() {
            self.by_sequence[usize::from(decoy)] += 1;
            return ids.into_iter().map(|id| (decoy, id)).collect();
        }
        if !decoy {
            self.unlinked_targets += 1;
            return Vec::new();
        }
        self.placeholder_psms += 1;
        let next = -1 - i64::try_from(self.placeholders.len()).unwrap_or(i64::MAX - 1);
        let id = *self
            .placeholders
            .entry(sequence.to_string())
            .or_insert(next);
        vec![(true, id)]
    }

    fn notes(&self) -> Vec<String> {
        let mut lines = Vec::new();
        let [targets, decoys] = self.by_sequence;
        if targets + decoys > 0 {
            lines.push(format!(
                "{} PSMs (targets {targets}, decoys {decoys}) that Proteome Discoverer linked \
                 to no protein were linked by sequence to the stored proteins containing them \
                 (decoys: reversed).",
                targets + decoys
            ));
        }
        if self.placeholder_psms > 0 {
            lines.push(format!(
                "{} decoy PSMs come from proteins the file doesn't store; each decoy peptide \
                 is filed under a placeholder decoy protein (DECOY_unstored:<peptide>) so they \
                 still count as decoys.",
                self.placeholder_psms
            ));
        }
        if self.unlinked_targets > 0 {
            lines.push(format!(
                "{} target PSMs map to no stored protein and have none.",
                self.unlinked_targets
            ));
        }
        lines
    }
}

/// Whether a peptide's MH+ lands within [`MASS_TOLERANCE`] of an isotope
/// peak (−1..+3) of the precursor. An unknown precursor can't be checked
/// and passes.
fn fits_precursor(calculated_mh: f64, precursor_mh: Option<f64>) -> bool {
    let Some(precursor) = precursor_mh.filter(|mh| mh.is_finite()) else {
        return true;
    };
    (-1..=3).any(|isotope| {
        (precursor - calculated_mh - f64::from(isotope) * ISOTOPE_SPACING).abs() <= MASS_TOLERANCE
    })
}

/// File id → run name.
fn runs(conn: &Connection) -> Result<HashMap<i64, String>, PdError> {
    let mut statement = conn.prepare("SELECT FileID, FileName FROM FileInfos")?;
    let runs = statement
        .query_map([], |row| {
            Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?))
        })?
        .map(|row| row.map(|(id, file)| (id, run_name(&file))))
        .collect::<rusqlite::Result<HashMap<_, _>>>()?;
    if runs.is_empty() {
        return Err(PdError::Layout(
            "no spectra files listed (FileInfos)".into(),
        ));
    }
    Ok(runs)
}

fn spectra(conn: &Connection) -> Result<HashMap<i64, SpectrumRow>, PdError> {
    let mut statement = conn.prepare(
        "SELECT s.SpectrumID, m.FileID, s.FirstScan, s.LastScan, s.RetentionTime, s.Charge, \
         s.Mass, m.Mass FROM SpectrumHeaders s JOIN MassPeaks m ON m.MassPeakID = s.MassPeakID",
    )?;
    let spectra = statement
        .query_map([], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                SpectrumRow {
                    file: row.get(1)?,
                    first: row.get(2)?,
                    last: row.get(3)?,
                    rt_minutes: row.get(4)?,
                    charge: row.get(5)?,
                    mh: row.get(6)?,
                    mz: row.get(7)?,
                },
            ))
        })?
        .collect::<rusqlite::Result<_>>()?;
    Ok(spectra)
}

/// Node number → engine name, for nodes that identified peptides.
fn engines(conn: &Connection) -> Result<HashMap<i64, String>, PdError> {
    let mut statement = conn.prepare(
        "SELECT ProcessingNodeNumber, FriendlyName FROM ProcessingNodes \
         WHERE ProcessingNodeNumber IN (SELECT ProcessingNodeNumber FROM Peptides \
         UNION SELECT ProcessingNodeNumber FROM Peptides_decoy)",
    )?;
    let nodes: Vec<(i64, String)> = statement
        .query_map([], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                row.get::<_, Option<String>>(1)?.unwrap_or_default(),
            ))
        })?
        .collect::<rusqlite::Result<_>>()?;
    let mut per_name: HashMap<&str, usize> = HashMap::new();
    for (_, name) in &nodes {
        *per_name.entry(name.as_str()).or_default() += 1;
    }
    Ok(nodes
        .iter()
        .map(|(number, name)| {
            let name = match name.as_str() {
                "" => format!("Proteome Discoverer node {number}"),
                name if per_name[name] > 1 => format!("{name} ({number})"),
                name => name.to_string(),
            };
            (*number, name)
        })
        .collect())
}

fn modification_types(conn: &Connection) -> Result<HashMap<i64, ModificationType>, PdError> {
    let mut statement = conn.prepare(
        "SELECT AminoAcidModificationID, ModificationName, DeltaMass, PositionType, \
         UnimodAccession FROM AminoAcidModifications",
    )?;
    let types = statement
        .query_map([], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                ModificationType {
                    name: row.get::<_, Option<String>>(1)?.unwrap_or_default(),
                    delta: row.get(2)?,
                    position_type: row.get::<_, Option<i64>>(3)?.unwrap_or(0),
                    unimod: row.get(4)?,
                },
            ))
        })?
        .collect::<rusqlite::Result<_>>()?;
    Ok(types)
}

/// Score id → name.
fn score_names(conn: &Connection) -> Result<HashMap<i64, String>, PdError> {
    let mut statement =
        conn.prepare("SELECT ScoreID, FriendlyName, ScoreName FROM ProcessingNodeScores")?;
    let names = statement
        .query_map([], |row| {
            let friendly: Option<String> = row.get(1)?;
            let name: Option<String> = row.get(2)?;
            Ok((row.get::<_, i64>(0)?, friendly.or(name).unwrap_or_default()))
        })?
        .collect::<rusqlite::Result<_>>()?;
    Ok(names)
}

fn scores(
    conn: &Connection,
    side: Side,
    names: &HashMap<i64, String>,
) -> Result<HashMap<PeptideKey, Vec<Score>>, PdError> {
    let mut statement = conn.prepare(&format!(
        "SELECT ProcessingNodeNumber, PeptideID, ScoreID, ScoreValue FROM PeptideScores{}",
        side.suffix()
    ))?;
    let mut scores: HashMap<PeptideKey, Vec<Score>> = HashMap::new();
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let key: PeptideKey = (row.get(0)?, row.get(1)?);
        let id: i64 = row.get(2)?;
        let Some(value) = row
            .get::<_, Option<f64>>(3)?
            .filter(|value| value.is_finite())
        else {
            continue;
        };
        let name = names.get(&id).ok_or_else(|| {
            PdError::Layout(format!("score id {id} isn't in ProcessingNodeScores"))
        })?;
        scores.entry(key).or_default().push(Score {
            name: name.clone(),
            accession: None,
            value,
        });
    }
    Ok(scores)
}

/// Peptide → (modification type, site).
fn modifications(
    conn: &Connection,
    side: Side,
    types: &HashMap<i64, ModificationType>,
) -> Result<HashMap<PeptideKey, Vec<(i64, ModSite)>>, PdError> {
    let suffix = side.suffix();
    let mut found: HashMap<PeptideKey, Vec<(i64, ModSite)>> = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT ProcessingNodeNumber, PeptideID, AminoAcidModificationID, Position \
         FROM PeptidesAminoAcidModifications{suffix}"
    ))?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let key: PeptideKey = (row.get(0)?, row.get(1)?);
        let id: i64 = row.get(2)?;
        let position: i64 = row.get(3)?;
        let residue = u32::try_from(position + 1)
            .map_err(|_| PdError::Layout(format!("modification position {position}")))?;
        found
            .entry(key)
            .or_default()
            .push((id, ModSite::Residue(residue)));
    }
    let mut statement = conn.prepare(&format!(
        "SELECT ProcessingNodeNumber, PeptideID, TerminalModificationID \
         FROM PeptidesTerminalModifications{suffix}"
    ))?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let key: PeptideKey = (row.get(0)?, row.get(1)?);
        let id: i64 = row.get(2)?;
        let site = match types.get(&id).map(|kind| kind.position_type) {
            Some(1 | 3) => ModSite::NTerm,
            Some(2 | 4) => ModSite::CTerm,
            other => {
                return Err(PdError::Layout(format!(
                    "terminal modification {id} has position type {other:?}"
                )));
            }
        };
        found.entry(key).or_default().push((id, site));
    }
    Ok(found)
}

fn peptide_form(
    sequence: &str,
    found: Option<&Vec<(i64, ModSite)>>,
    types: &HashMap<i64, ModificationType>,
) -> Result<PeptideForm, PdError> {
    let length = u32::try_from(sequence.len()).unwrap_or(u32::MAX);
    let mut modifications = Vec::new();
    for &(id, site) in found.into_iter().flatten() {
        let kind = types
            .get(&id)
            .ok_or_else(|| PdError::Layout(format!("modification {id} isn't listed")))?;
        if let ModSite::Residue(residue) = site
            && !(1..=length).contains(&residue)
        {
            return Err(PdError::Layout(format!(
                "{} at residue {residue} doesn't fit {sequence}",
                kind.name
            )));
        }
        modifications.push(Modification {
            site,
            mass_delta: kind.delta,
            name: Some(kind.name.clone()).filter(|name| !name.is_empty()),
            accession: kind
                .unimod
                .filter(|&id| id > 0)
                .map(|id| format!("UNIMOD:{id}")),
        });
    }
    PeptideForm::new(sequence, modifications).map_err(|error| PdError::Layout(error.to_string()))
}

fn protein_links(conn: &Connection, side: Side) -> Result<HashMap<PeptideKey, Vec<i64>>, PdError> {
    let mut statement = conn.prepare(&format!(
        "SELECT ProcessingNodeNumber, PeptideID, ProteinID FROM PeptidesProteins{}",
        side.suffix()
    ))?;
    let mut links: HashMap<PeptideKey, Vec<i64>> = HashMap::new();
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        links
            .entry((row.get(0)?, row.get(1)?))
            .or_default()
            .push(row.get(2)?);
    }
    Ok(links)
}

/// (id, accession, description, FASTA title line, sequence).
type StoredProtein = (i64, String, Option<String>, String, Option<String>);

/// The proteins PSMs name, target and decoy, under unique accessions.
fn proteins(
    conn: &Connection,
    used: &HashSet<i64>,
) -> Result<HashMap<(bool, i64), Protein>, PdError> {
    let mut statement = conn.prepare(
        "SELECT p.ProteinID, p.Sequence, a.Description FROM Proteins p \
         LEFT JOIN ProteinAnnotations a ON a.ProteinID = p.ProteinID",
    )?;
    let mut raw: Vec<StoredProtein> = Vec::new();
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        if !used.contains(&id) {
            continue;
        }
        let sequence: Option<String> = row.get(1)?;
        let title = row
            .get::<_, Option<String>>(2)?
            .map(|title| {
                title
                    .lines()
                    .next()
                    .unwrap_or("")
                    .trim_start_matches('>')
                    .trim()
                    .to_string()
            })
            .filter(|title| !title.is_empty())
            .ok_or_else(|| PdError::Layout(format!("protein {id} has no FASTA title line")))?;
        let (accession, description) = match title.split_once(char::is_whitespace) {
            Some((accession, rest)) => (accession.to_string(), Some(rest.trim().to_string())),
            None => (title.clone(), None),
        };
        raw.push((
            id,
            accession,
            description.filter(|d| !d.is_empty()),
            title,
            sequence,
        ));
    }
    let mut uses: HashMap<&str, usize> = HashMap::new();
    for (_, accession, ..) in &raw {
        *uses.entry(accession.as_str()).or_default() += 1;
    }
    let shared: HashSet<String> = uses
        .into_iter()
        .filter(|&(_, count)| count > 1)
        .map(|(accession, _)| accession.to_string())
        .collect();
    let mut proteins = HashMap::new();
    for (id, accession, description, title, sequence) in raw {
        let name = if shared.contains(&accession) {
            title
        } else {
            accession
        };
        proteins.insert(
            (true, id),
            Protein {
                accession: format!("{DECOY_PREFIX}{name}"),
                description: None,
                sequence: None,
                is_decoy: true,
                decoy_explicit: true,
            },
        );
        proteins.insert(
            (false, id),
            Protein {
                accession: name,
                description,
                sequence: sequence
                    .filter(|sequence| !sequence.is_empty())
                    .map(|sequence| sequence.to_ascii_uppercase()),
                is_decoy: false,
                decoy_explicit: true,
            },
        );
    }
    Ok(proteins)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn precursors_fit_within_an_isotope_peak() {
        // Monoisotopic, +1 and +2 isotope picks, a −1 pick.
        assert!(fits_precursor(1000.0, Some(1000.02)));
        assert!(fits_precursor(1000.0, Some(1001.0034)));
        assert!(fits_precursor(1000.0, Some(1002.0067)));
        assert!(fits_precursor(1000.0, Some(998.9966)));
        // Half a dalton off, or a missing modification (+16): no.
        assert!(!fits_precursor(1000.0, Some(1000.5)));
        assert!(!fits_precursor(1000.0, Some(1015.9949)));
        assert!(fits_precursor(1000.0, None));
    }
}
