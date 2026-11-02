//! Thermo Proteome Discoverer 2.x/3.x results (`.pdResult`, and the
//! per-search `.msf` files of those versions) → one [`ImportBatch`] per
//! search node.
//!
//! Both files are SQLite databases with `TargetPsms` and `DecoyPsms`
//! tables of the same layout; see `docs/pd-format.md` for what was checked
//! and how. Reading only, and nothing is guessed:
//! - Each identifying search node is one engine (`Sequest HT`, or
//!   `Sequest HT (6)` when a workflow ran several of the same type). Ranks
//!   are the engine's own (`SearchEngineRank`), not PD's.
//! - Score columns are those `ProcessingNodeScores` names, plus `DeltaCn`
//!   and `DeltaScore`, when both PSM tables have them; PD's `Confidence`
//!   and any validator output are not imported (Tuffold validates).
//! - Modification positions are 1-based; position types 1 and 3 are
//!   N-terminal, 2 and 4 C-terminal. Every PSM's mass (sequence plus
//!   modifications) must match PD's own theoretical mass, or the file is
//!   refused: that check catches any misread modification.
//! - PD files list decoys under their target's accession, in their own
//!   tables; they are imported as `DECOY_<accession>`, flagged explicitly.
//!   Accessions shared by different proteins (FASTA headers with no unique
//!   first word) use the full FASTA title line instead.
//! - Each spectra file is one MS run (and biosample until grouped); a
//!   spectrum is keyed by its scan number, or by its scan range when it
//!   summed several.
//! - PD 1.x `.msf` (a `Peptides` table) has its own layout, read by [`v1`]
//!   under the same rules.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use rusqlite::{Connection, OpenFlags, Row};
use tuffold_core::chemistry::{PROTON, peptide_mass};
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, Spectrum, SpectrumKey,
};

use crate::xml::run_name;

mod v1;

/// PD reports theoretical MH+ to this many decimals; a misread
/// modification is off by at least a dalton.
const MASS_TOLERANCE: f64 = 0.01;
/// Score columns that aren't in `ProcessingNodeScores` but are engine
/// output all the same.
const GENERIC_SCORES: [(&str, &str); 2] = [("DeltaCn", "ΔCn"), ("DeltaScore", "ΔScore")];
const DECOY_PREFIX: &str = "DECOY_";

#[derive(Debug, thiserror::Error)]
pub enum PdError {
    #[error("not a Proteome Discoverer result: {0}")]
    NotPd(String),
    #[error("cannot read the Proteome Discoverer file: {0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error(
        "{0}. The file's layout differs from the Proteome Discoverer versions Tuffold knows \
         (1.4, 2.x–3.x), so it was not imported rather than guessed at"
    )]
    Layout(String),
}

/// One batch per search node, and notes on what was left out.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PdImport {
    pub batches: Vec<ImportBatch>,
    pub notes: Vec<String>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Side {
    Target,
    Decoy,
}

impl Side {
    fn psms(self) -> &'static str {
        match self {
            Self::Target => "TargetPsms",
            Self::Decoy => "DecoyPsms",
        }
    }

    fn proteins(self) -> &'static str {
        match self {
            Self::Target => "TargetProteins",
            Self::Decoy => "DecoyProteins",
        }
    }
}

struct Draft {
    engine: String,
    run: String,
    key: SpectrumKey,
    spectrum: Spectrum,
    /// Evidence is filled in from these (decoy?, protein id) when batched.
    proteins: Vec<(bool, i64)>,
    psm: Psm,
}

/// What every PSM row is read against.
struct Context<'a> {
    runs: &'a HashMap<i64, String>,
    scores: &'a [(String, String)],
    engine_names: &'a HashMap<String, String>,
    modifications: &'a HashMap<i64, FoundModification>,
}

pub fn import_pd(path: &Path) -> Result<PdImport, PdError> {
    let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)
        .map_err(|error| PdError::NotPd(error.to_string()))?;
    let tables = table_names(&conn).map_err(|error| PdError::NotPd(error.to_string()))?;
    let pd1 = tables.contains("Peptides") && tables.contains("SpectrumHeaders");
    if !tables.contains("TargetPsms") && !pd1 {
        return Err(PdError::NotPd("no TargetPsms or Peptides table".into()));
    }
    let version: Option<String> = conn
        .query_row(
            "SELECT SoftwareVersion FROM SchemaInfo WHERE Kind = 'Result' LIMIT 1",
            [],
            |row| row.get(0),
        )
        .ok();
    if !tables.contains("TargetPsms") {
        return v1::import(&conn, version.as_deref());
    }
    let runs = runs(&conn)?;
    let scores = score_columns(&conn, &tables)?;
    let engine_names = engine_names(&conn, &tables)?;
    let modifications = modification_table(&conn)?;
    let mut notes = Notes::default();

    let context = Context {
        runs: &runs,
        scores: &scores,
        engine_names: &engine_names,
        modifications: &modifications,
    };
    let mut drafts = Vec::new();
    for side in [Side::Target, Side::Decoy] {
        if !tables.contains(side.psms()) {
            continue;
        }
        let proteins = proteins_of(&conn, side, &tables)?;
        let psm_mods = psm_modifications(&conn, side, &modifications)?;
        read_psms(
            &conn,
            side,
            &context,
            &proteins,
            &psm_mods,
            &mut drafts,
            &mut notes,
        )?;
    }
    if drafts.is_empty() {
        return Err(PdError::Layout("no PSMs".into()));
    }
    let proteins = all_proteins(&conn, &tables)?;
    Ok(PdImport {
        batches: batches(drafts, &proteins, version.as_deref()),
        notes: notes.lines(),
    })
}

fn table_names(conn: &Connection) -> rusqlite::Result<HashSet<String>> {
    let mut statement = conn.prepare("SELECT name FROM sqlite_master WHERE type = 'table'")?;
    let names = statement
        .query_map([], |row| row.get(0))?
        .collect::<rusqlite::Result<_>>()?;
    Ok(names)
}

fn columns(conn: &Connection, table: &str) -> Result<HashSet<String>, PdError> {
    let mut statement = conn.prepare(&format!("PRAGMA table_info(\"{table}\")"))?;
    let names = statement
        .query_map([], |row| row.get::<_, String>(1))?
        .collect::<rusqlite::Result<_>>()?;
    Ok(names)
}

/// Spectra file id → run name.
fn runs(conn: &Connection) -> Result<HashMap<i64, String>, PdError> {
    let mut statement = conn.prepare("SELECT FileID, FileName FROM WorkflowInputFiles")?;
    let rows = statement.query_map([], |row| {
        Ok((row.get::<_, i64>(0)?, row.get::<_, String>(1)?))
    })?;
    let mut runs = HashMap::new();
    for row in rows {
        let (id, file) = row?;
        runs.insert(id, run_name(&file));
    }
    if runs.is_empty() {
        return Err(PdError::Layout(
            "no spectra files listed (WorkflowInputFiles)".into(),
        ));
    }
    Ok(runs)
}

/// (column, display name) of every engine score both PSM tables carry.
fn score_columns(
    conn: &Connection,
    tables: &HashSet<String>,
) -> Result<Vec<(String, String)>, PdError> {
    let target = columns(conn, "TargetPsms")?;
    let decoy = if tables.contains("DecoyPsms") {
        columns(conn, "DecoyPsms")?
    } else {
        target.clone()
    };
    let mut scores: Vec<(String, String)> = Vec::new();
    if tables.contains("ProcessingNodeScores") {
        let mut statement = conn.prepare(
            "SELECT ScoreName, COALESCE(DisplayName, ScoreName) FROM ProcessingNodeScores \
             ORDER BY IsMainScore DESC, ScoreID",
        )?;
        for row in statement.query_map([], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
        })? {
            let (column, name) = row?;
            if !scores.iter().any(|(known, _)| *known == column) {
                scores.push((column, name));
            }
        }
    }
    scores.extend(
        GENERIC_SCORES
            .iter()
            .map(|(column, name)| (column.to_string(), name.to_string())),
    );
    scores.retain(|(column, _)| {
        target.contains(column) && decoy.contains(column) && is_identifier(column)
    });
    Ok(scores)
}

/// Column names are spliced into SQL, so only plain identifiers are used.
fn is_identifier(name: &str) -> bool {
    !name.is_empty()
        && name
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || ch == '_')
}

/// Engine name per identifying node: its type, or its node name when one
/// workflow ran several nodes of that type.
fn engine_names(
    conn: &Connection,
    tables: &HashSet<String>,
) -> Result<HashMap<String, String>, PdError> {
    let mut nodes: Vec<(String, String)> = Vec::new();
    for side in [Side::Target, Side::Decoy] {
        if !tables.contains(side.psms()) {
            continue;
        }
        let mut statement = conn.prepare(&format!(
            "SELECT DISTINCT IdentifyingNodeName, IdentifyingNodeTypeName FROM {}",
            side.psms()
        ))?;
        for row in statement.query_map([], |row| {
            Ok((
                row.get::<_, Option<String>>(0)?.unwrap_or_default(),
                row.get::<_, Option<String>>(1)?.unwrap_or_default(),
            ))
        })? {
            let node = row?;
            if !nodes.contains(&node) {
                nodes.push(node);
            }
        }
    }
    let mut per_type: HashMap<&str, usize> = HashMap::new();
    for (_, kind) in &nodes {
        *per_type.entry(kind.as_str()).or_default() += 1;
    }
    Ok(nodes
        .iter()
        .map(|(node, kind)| {
            let name = if kind.is_empty() || per_type[kind.as_str()] > 1 {
                node
            } else {
                kind
            };
            (
                node.clone(),
                if name.is_empty() {
                    "Proteome Discoverer".into()
                } else {
                    name.clone()
                },
            )
        })
        .collect())
}

struct FoundModification {
    name: String,
    delta: f64,
    position_type: i64,
    unimod: Option<i64>,
}

fn modification_table(conn: &Connection) -> Result<HashMap<i64, FoundModification>, PdError> {
    let mut statement = conn.prepare(
        "SELECT ModificationID, Name, DeltaMonoisotopicMass, PositionType, UnimodAccession \
         FROM FoundModifications",
    )?;
    let rows = statement.query_map([], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            FoundModification {
                name: row.get::<_, Option<String>>(1)?.unwrap_or_default(),
                delta: row.get(2)?,
                position_type: row.get(3)?,
                unimod: row.get(4)?,
            },
        ))
    })?;
    rows.collect::<rusqlite::Result<_>>().map_err(PdError::from)
}

type PsmKey = (i64, i64);

/// (workflow, peptide) → modifications at their sites.
fn psm_modifications(
    conn: &Connection,
    side: Side,
    table: &HashMap<i64, FoundModification>,
) -> Result<HashMap<PsmKey, Vec<(i64, i64)>>, PdError> {
    let link = format!("{}FoundModifications", side.psms());
    let mut statement = conn.prepare(&format!(
        "SELECT {0}WorkflowID, {0}PeptideID, FoundModificationsModificationID, Position FROM {link}",
        side.psms()
    ))?;
    let mut out: HashMap<PsmKey, Vec<(i64, i64)>> = HashMap::new();
    for row in statement.query_map([], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            row.get::<_, i64>(1)?,
            row.get::<_, i64>(2)?,
            row.get::<_, i64>(3)?,
        ))
    })? {
        let (workflow, peptide, modification, position) = row?;
        if !table.contains_key(&modification) {
            return Err(PdError::Layout(format!(
                "modification {modification} is not listed"
            )));
        }
        out.entry((workflow, peptide))
            .or_default()
            .push((modification, position));
    }
    Ok(out)
}

/// (workflow, peptide) → the protein ids (UniqueSequenceID) it maps to.
fn proteins_of(
    conn: &Connection,
    side: Side,
    tables: &HashSet<String>,
) -> Result<HashMap<PsmKey, Vec<i64>>, PdError> {
    let link = format!("{}{}", side.proteins(), side.psms());
    if !tables.contains(&link) {
        return Ok(HashMap::new());
    }
    let mut statement = conn.prepare(&format!(
        "SELECT {psms}WorkflowID, {psms}PeptideID, {proteins}UniqueSequenceID FROM {link}",
        psms = side.psms(),
        proteins = side.proteins()
    ))?;
    let mut out: HashMap<PsmKey, Vec<i64>> = HashMap::new();
    for row in statement.query_map([], |row| {
        Ok((
            row.get::<_, i64>(0)?,
            row.get::<_, i64>(1)?,
            row.get::<_, i64>(2)?,
        ))
    })? {
        let (workflow, peptide, protein) = row?;
        out.entry((workflow, peptide)).or_default().push(protein);
    }
    Ok(out)
}

/// Every protein, target and decoy, under a unique accession.
fn all_proteins(
    conn: &Connection,
    tables: &HashSet<String>,
) -> Result<HashMap<(bool, i64), Protein>, PdError> {
    // (decoy?, id, accession, description, sequence, FASTA title lines)
    type Raw = (
        bool,
        i64,
        String,
        Option<String>,
        Option<String>,
        Option<String>,
    );
    let mut raw: Vec<Raw> = Vec::new();
    for side in [Side::Target, Side::Decoy] {
        if !tables.contains(side.proteins()) {
            continue;
        }
        let mut statement = conn.prepare(&format!(
            "SELECT UniqueSequenceID, Accession, Description, Sequence, FastaTitleLines FROM {}",
            side.proteins()
        ))?;
        for row in statement.query_map([], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                row.get::<_, Option<String>>(1)?.unwrap_or_default(),
                row.get::<_, Option<String>>(2)?,
                row.get::<_, Option<String>>(3)?,
                row.get::<_, Option<String>>(4)?,
            ))
        })? {
            let (id, accession, description, sequence, title) = row?;
            raw.push((
                side == Side::Decoy,
                id,
                accession,
                description,
                sequence,
                title,
            ));
        }
    }
    // Accessions shared by different proteins of one side get their title line.
    let mut uses: HashMap<(bool, &str), HashSet<i64>> = HashMap::new();
    for (decoy, id, accession, ..) in &raw {
        uses.entry((*decoy, accession.as_str()))
            .or_default()
            .insert(*id);
    }
    let shared: HashSet<(bool, String)> = uses
        .into_iter()
        .filter(|(_, ids)| ids.len() > 1)
        .map(|((decoy, accession), _)| (decoy, accession.to_string()))
        .collect();
    let mut proteins = HashMap::new();
    for (decoy, id, accession, description, sequence, title) in raw {
        let mut name = accession.clone();
        if shared.contains(&(decoy, accession)) {
            let first_title = title
                .as_deref()
                .and_then(|title| title.lines().next())
                .map(|line| line.trim_start_matches('>').trim().to_string())
                .filter(|line| !line.is_empty());
            name = first_title.ok_or_else(|| {
                PdError::Layout(format!(
                    "several proteins share accession {name} and have no title"
                ))
            })?;
        }
        let accession = if decoy {
            format!("{DECOY_PREFIX}{name}")
        } else {
            name
        };
        proteins.insert(
            (decoy, id),
            Protein {
                accession,
                description,
                sequence: sequence
                    .filter(|sequence| !sequence.is_empty())
                    .map(|s| s.to_ascii_uppercase()),
                is_decoy: decoy,
                decoy_explicit: true,
            },
        );
    }
    Ok(proteins)
}

fn read_psms(
    conn: &Connection,
    side: Side,
    context: &Context<'_>,
    proteins: &HashMap<PsmKey, Vec<i64>>,
    psm_mods: &HashMap<PsmKey, Vec<(i64, i64)>>,
    drafts: &mut Vec<Draft>,
    notes: &mut Notes,
) -> Result<(), PdError> {
    let score_sql: String = context
        .scores
        .iter()
        .map(|(column, _)| format!(", \"{column}\""))
        .collect();
    let mut statement = conn.prepare(&format!(
        "SELECT WorkflowID, PeptideID, IdentifyingNodeName, Sequence, Charge, MassOverCharge, \
         TheoreticalMass, RetentionTime, FirstScan, LastScan, SpectrumFileId, SearchEngineRank\
         {score_sql} FROM {}",
        side.psms()
    ))?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let key: PsmKey = (row.get(0)?, row.get(1)?);
        let node: String = row.get::<_, Option<String>>(2)?.unwrap_or_default();
        let sequence: String = row
            .get::<_, Option<String>>(3)?
            .unwrap_or_default()
            .to_ascii_uppercase();
        let charge: Option<i64> = row.get(4)?;
        let mz: Option<f64> = row.get(5)?;
        let theoretical: Option<f64> = row.get(6)?;
        let rt_minutes: Option<f64> = row.get(7)?;
        let (first, last): (Option<i64>, Option<i64>) = (row.get(8)?, row.get(9)?);
        let file: Option<i64> = row.get(10)?;
        let rank: Option<i64> = row.get(11)?;
        let Some(run) = file.and_then(|file| context.runs.get(&file)) else {
            notes.without_file += 1;
            continue;
        };
        let Some(charge) = charge
            .and_then(|charge| i16::try_from(charge).ok())
            .filter(|&c| c > 0)
        else {
            notes.without_charge += 1;
            continue;
        };
        if sequence.is_empty() || !sequence.bytes().all(|byte| byte.is_ascii_uppercase()) {
            return Err(PdError::Layout(format!(
                "PSM {key:?} has sequence {sequence:?}"
            )));
        }
        let peptide = peptide_form(&sequence, psm_mods.get(&key), context.modifications)?;
        check_mass(&peptide, theoretical, notes)?;
        let scores = score_values(row, context.scores, 12)?;
        let protein_ids = proteins
            .get(&key)
            .into_iter()
            .flatten()
            .map(|&protein| (side == Side::Decoy, protein))
            .collect();
        let Some(key) = spectrum_key(first, last) else {
            notes.without_scan += 1;
            continue;
        };
        let z = f64::from(charge);
        let experimental_mz = mz.filter(|mz| mz.is_finite() && *mz > 0.0);
        drafts.push(Draft {
            engine: context
                .engine_names
                .get(&node)
                .cloned()
                .unwrap_or_else(|| node.clone()),
            run: run.clone(),
            key: key.clone(),
            spectrum: Spectrum {
                ms_run: run.clone(),
                key: key.clone(),
                source_id: None,
                precursor_mz: experimental_mz,
                precursor_charge: Some(charge),
                rt_seconds: rt_minutes
                    .filter(|rt| rt.is_finite() && *rt >= 0.0)
                    .map(|rt| rt * 60.0),
            },
            proteins: protein_ids,
            psm: Psm {
                ms_run: run.clone(),
                spectrum: key,
                rank: rank
                    .and_then(|rank| u32::try_from(rank).ok())
                    .filter(|&rank| rank > 0)
                    .unwrap_or(1),
                charge,
                calculated_mz: theoretical.map(|mh| (mh - PROTON) / z + PROTON),
                experimental_mz,
                pass_threshold: None,
                scores,
                evidence: Vec::new(),
                peptide,
            },
        });
    }
    Ok(())
}

/// A spectrum's key: its scan number, or its scan range when it summed
/// several.
fn spectrum_key(first: Option<i64>, last: Option<i64>) -> Option<SpectrumKey> {
    match (first, last) {
        (Some(first), Some(last)) if first == last => Some(u32::try_from(first).map_or_else(
            |_| SpectrumKey::Native(format!("scan {first}")),
            SpectrumKey::Scan,
        )),
        (Some(first), Some(last)) => Some(SpectrumKey::Native(format!("scans {first}-{last}"))),
        _ => None,
    }
}

fn score_values(
    row: &Row<'_>,
    scores: &[(String, String)],
    first: usize,
) -> Result<Vec<Score>, PdError> {
    let mut values = Vec::new();
    for (index, (_, name)) in scores.iter().enumerate() {
        let value: Option<f64> = row.get(first + index)?;
        if let Some(value) = value.filter(|value| value.is_finite()) {
            values.push(Score {
                name: name.clone(),
                accession: None,
                value,
            });
        }
    }
    Ok(values)
}

fn peptide_form(
    sequence: &str,
    found: Option<&Vec<(i64, i64)>>,
    table: &HashMap<i64, FoundModification>,
) -> Result<PeptideForm, PdError> {
    let length = i64::try_from(sequence.len()).unwrap_or(i64::MAX);
    let mut modifications = Vec::new();
    for &(id, position) in found.into_iter().flatten() {
        let modification = &table[&id];
        let site = match modification.position_type {
            0 if (1..=length).contains(&position) => {
                ModSite::Residue(u32::try_from(position).expect("checked range"))
            }
            1 | 3 => ModSite::NTerm,
            2 | 4 => ModSite::CTerm,
            other => {
                return Err(PdError::Layout(format!(
                    "{} at {position} (position type {other}) doesn't fit {sequence}",
                    modification.name
                )));
            }
        };
        modifications.push(Modification {
            site,
            mass_delta: modification.delta,
            name: Some(modification.name.clone()).filter(|name| !name.is_empty()),
            accession: modification
                .unimod
                .filter(|&id| id > 0)
                .map(|id| format!("UNIMOD:{id}")),
        });
    }
    PeptideForm::new(sequence, modifications).map_err(|error| PdError::Layout(error.to_string()))
}

/// Sequence plus modifications must give PD's theoretical MH+.
fn check_mass(
    peptide: &PeptideForm,
    theoretical: Option<f64>,
    notes: &mut Notes,
) -> Result<(), PdError> {
    let (Some(theoretical), Some(mass)) = (theoretical, peptide_mass(&peptide.sequence)) else {
        notes.unchecked_mass += 1;
        return Ok(());
    };
    let deltas: f64 = peptide.modifications.iter().map(|m| m.mass_delta).sum();
    let error = mass + deltas + PROTON - theoretical;
    if error.is_nan() || error.abs() > MASS_TOLERANCE {
        return Err(PdError::Layout(format!(
            "{} weighs {:.4} Da less than PD's theoretical MH+ {theoretical:.4}",
            peptide.proforma(),
            -error
        )));
    }
    Ok(())
}

fn batches(
    drafts: Vec<Draft>,
    proteins: &HashMap<(bool, i64), Protein>,
    version: Option<&str>,
) -> Vec<ImportBatch> {
    // Only spectra files PSMs name: a .pdResult also lists its input .msf
    // files among WorkflowInputFiles.
    let mut run_names: Vec<String> = drafts.iter().map(|draft| draft.run.clone()).collect();
    run_names.sort();
    run_names.dedup();
    let mut by_engine: BTreeMap<String, Vec<Draft>> = BTreeMap::new();
    for draft in drafts {
        by_engine
            .entry(draft.engine.clone())
            .or_default()
            .push(draft);
    }

    by_engine
        .into_iter()
        .map(|(engine, drafts)| {
            let mut spectra: Vec<Spectrum> = Vec::new();
            let mut seen: HashMap<(String, SpectrumKey), usize> = HashMap::new();
            let mut used: BTreeMap<String, Protein> = BTreeMap::new();
            let mut psms = Vec::with_capacity(drafts.len());
            for mut draft in drafts {
                match seen.get(&(draft.run.clone(), draft.key.clone())) {
                    Some(&index) => spectra[index].merge(&draft.spectrum),
                    None => {
                        seen.insert((draft.run.clone(), draft.key.clone()), spectra.len());
                        spectra.push(draft.spectrum);
                    }
                }
                for id in &draft.proteins {
                    let Some(protein) = proteins.get(id) else {
                        continue;
                    };
                    if !draft
                        .psm
                        .evidence
                        .iter()
                        .any(|e| e.accession == protein.accession)
                    {
                        draft.psm.evidence.push(Evidence {
                            accession: protein.accession.clone(),
                            start: None,
                            end: None,
                            pre: None,
                            post: None,
                        });
                    }
                    used.entry(protein.accession.clone())
                        .or_insert_with(|| protein.clone());
                }
                psms.push(draft.psm);
            }
            ImportBatch {
                software: Software {
                    name: engine,
                    version: version.map(|version| format!("Proteome Discoverer {version}")),
                },
                biosamples: run_names
                    .iter()
                    .map(|run| BioSample {
                        name: (*run).clone(),
                        category: None,
                    })
                    .collect(),
                ms_runs: run_names
                    .iter()
                    .map(|run| MsRun {
                        name: (*run).clone(),
                        biosample: (*run).clone(),
                        spectra_file: None,
                    })
                    .collect(),
                spectra,
                proteins: used.into_values().collect(),
                psms,
            }
        })
        .collect()
}

#[derive(Default)]
struct Notes {
    without_file: usize,
    without_charge: usize,
    without_scan: usize,
    unchecked_mass: usize,
    /// PD 1.x: PSMs whose mass misses their precursor's.
    mass_misses: usize,
}

impl Notes {
    fn lines(&self) -> Vec<String> {
        let mut lines = vec![
            "Proteome Discoverer's own confidences and validator results were not imported: \
             Tuffold validates the search engines' scores itself."
                .to_string(),
        ];
        for (count, what) in [
            (self.without_file, "PSMs naming no listed spectra file"),
            (self.without_charge, "PSMs without a precursor charge"),
            (self.without_scan, "PSMs without a scan number"),
        ] {
            if count > 0 {
                lines.push(format!("{count} {what} were left out."));
            }
        }
        if self.unchecked_mass > 0 {
            lines.push(format!(
                "{} PSMs with non-standard residues could not be mass-checked.",
                self.unchecked_mass
            ));
        }
        if self.mass_misses > 0 {
            lines.push(format!(
                "{} PSMs whose mass misses their precursor's (allowing isotope peaks) were \
                 left out.",
                self.mass_misses
            ));
        }
        lines
    }
}
