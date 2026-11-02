//! Read-side queries behind the browsing views: protein groups at the saved
//! thresholds (with spectral counts per biosample and coverage), a group's
//! peptides, a peptide's spectra, and one spectrum in full.
//!
//! Counts follow the thresholds: a group's spectra are the target spectra
//! passing the PSM FDR whose peptide is unique or razor to the group;
//! coverage uses every peptide of the group (shared ones too) on the
//! leader's sequence, and is `None` when that sequence is unknown.

use std::collections::HashMap;

use rusqlite::{Connection, OptionalExtension, params};
use tuffold_core::coverage::{SequenceCoverage, covered_fraction, map_peptide_coverage};
use tuffold_core::digest::observable_peptides;
use tuffold_core::{ModSite, Modification, PeptideForm};
use tuffold_validate::{Level, SiteScore, Thresholds, format_sites};

use crate::peaks::{self, Peaks};
use crate::{StoreError, validation};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Biosample {
    pub id: i64,
    pub name: String,
    /// Its group for comparisons (e.g. "control"), when assigned.
    pub category: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ProteinListRow {
    pub group_id: i64,
    pub accession: String,
    pub description: Option<String>,
    /// Indistinguishable proteins besides the leader.
    pub equivalents: usize,
    /// Distinct unique and razor peptide sequences.
    pub peptides: usize,
    /// Passing target spectra of the group's unique and razor peptides.
    pub spectra: usize,
    /// Parallel to [`ProteinList::biosamples`].
    pub spectra_by_sample: Vec<usize>,
    /// Passing spectra of the group's unique peptides only (Scaffold's
    /// exclusive spectrum count), per biosample.
    pub exclusive_spectra_by_sample: Vec<usize>,
    /// Distinct unique peptide sequences with a passing spectrum, per
    /// biosample.
    pub exclusive_peptides_by_sample: Vec<usize>,
    /// Distinct unique and razor peptide sequences with a passing spectrum,
    /// per biosample (emPAI's observed peptides).
    pub peptides_by_sample: Vec<usize>,
    pub q_value: Option<f64>,
    /// 1 − PEP; `None` when the group lost its competition.
    pub probability: Option<f64>,
    pub score: f64,
    /// Fraction of the leader's sequence covered; `None` without a sequence.
    pub coverage: Option<f64>,
    /// The leader's sequence length, when known.
    pub length: Option<usize>,
    /// The leader's observable tryptic peptides (for emPAI), when its
    /// sequence is known.
    pub observable_peptides: Option<usize>,
    /// Modifications on the group's unique and razor peptides: Unimod names,
    /// or the mass delta (`+15.995`) where a source gave none; sorted.
    pub modifications: Vec<String>,
    /// The user starred the group's leader.
    pub starred: bool,
    /// GO slim terms the leader's annotations reach, as `name (aspect)`,
    /// sorted; empty without a GO import.
    pub go: Vec<String>,
}

/// Protein groups passing the saved thresholds, most spectra first.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ProteinList {
    pub biosamples: Vec<Biosample>,
    pub rows: Vec<ProteinListRow>,
    /// The thresholds the rows pass (their mode decides what views show).
    pub thresholds: Thresholds,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroupMember {
    pub protein_id: i64,
    pub accession: String,
    pub description: Option<String>,
    /// `leader`, `equivalent`, or `subset`.
    pub role: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GroupPeptide {
    pub peptide_id: i64,
    pub proforma: String,
    pub sequence: String,
    /// `unique`, `razor`, or `shared`.
    pub role: String,
    pub q_value: Option<f64>,
    /// 1 − PEP.
    pub probability: Option<f64>,
    pub score: Option<f64>,
    /// Passing target spectra assigned this peptide.
    pub spectra: usize,
    /// Highest-scoring localized isoform ProForma from the best-localized
    /// passing spectrum of this peptide, when ambiguous PTM sites were scored.
    pub localized_peptide: Option<String>,
    /// Compact target site probability summary (e.g. `"S2: 99.0%, T4: 1.0%"`).
    pub site_probabilities: Option<String>,
    /// Best Ascore on the best-localized passing spectrum of this peptide.
    pub best_ascore: Option<f64>,
    /// Minimum assigned site probability on the best-localized passing spectrum.
    pub min_site_probability: Option<f64>,
    /// Target candidate site scores on the best-localized passing spectrum.
    pub site_scores: Vec<SiteScore>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GroupDetail {
    pub group_id: i64,
    pub members: Vec<GroupMember>,
    /// The leader's sequence, when known.
    pub sequence: Option<String>,
    pub peptides: Vec<GroupPeptide>,
    /// Where the peptides fall on `sequence` (empty without one).
    pub coverage: Vec<SequenceCoverage>,
    /// Other passing groups holding any of this group's peptides, most
    /// shared first (Scaffold's similarity view).
    pub similar: Vec<SimilarGroup>,
    /// Passing groups linked to this one through shared peptides, directly
    /// or via others, itself included (its protein cluster).
    pub cluster_size: usize,
}

/// A passing group sharing peptides with the one shown.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimilarGroup {
    pub group_id: i64,
    pub accession: String,
    pub description: Option<String>,
    /// Which of the shown group's peptides it holds.
    pub shared: Vec<i64>,
    /// Its own distinct unique and razor sequences.
    pub peptides: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumListRow {
    pub spectrum_id: i64,
    pub run: String,
    pub biosample: String,
    /// `scan=N`, `index=N`, or a native id.
    pub key: String,
    pub charge: Option<i64>,
    pub precursor_mz: Option<f64>,
    pub rt_seconds: Option<f64>,
    pub q_value: f64,
    /// 1 − PEP.
    pub probability: f64,
    pub score: f64,
    pub engines: u32,
    pub has_peaks: bool,
    /// Whether the spectrum passes the saved PSM FDR.
    pub passes: bool,
    pub localized_peptide: Option<String>,
    pub site_probabilities: Option<String>,
    pub best_ascore: Option<f64>,
    pub min_site_probability: Option<f64>,
    pub site_scores: Vec<SiteScore>,
}

/// The best spectra of a peptide, up to a limit, and how many it has in
/// all — a capped list always says so (`rows.len() < total`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SpectrumPage {
    pub rows: Vec<SpectrumListRow>,
    pub total: usize,
}

impl SpectrumPage {
    pub fn is_capped(&self) -> bool {
        self.rows.len() < self.total
    }
}

/// What is known about a spectrum's measured peaks.
#[derive(Debug, Clone, PartialEq)]
pub enum StoredPeaks {
    /// No spectra file was attached for this spectrum.
    Missing,
    Stored(Peaks),
    /// Stored but unreadable (damaged file).
    Unreadable,
}

/// One engine's rank-1 PSM for a spectrum.
#[derive(Debug, Clone, PartialEq)]
pub struct EnginePsm {
    pub engine: String,
    pub proforma: String,
    pub charge: i64,
    pub scores: Vec<(String, f64)>,
    /// From validation; `None` when the engine was not ranked or lost.
    pub q_value: Option<f64>,
    pub pep: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AssignedPeptide {
    pub peptide_id: i64,
    pub form: PeptideForm,
    pub is_decoy: bool,
    pub q_value: f64,
    /// 1 − PEP.
    pub probability: f64,
    pub score: f64,
    /// The charge to annotate with: the assigned PSM's, else the precursor's.
    pub charge: Option<u8>,
    pub localized_peptide: Option<String>,
    pub site_probabilities: Option<String>,
    pub best_ascore: Option<f64>,
    pub min_site_probability: Option<f64>,
    pub site_scores: Vec<SiteScore>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumDetail {
    pub spectrum_id: i64,
    pub run: String,
    pub key: String,
    pub precursor_mz: Option<f64>,
    pub precursor_charge: Option<i64>,
    pub rt_seconds: Option<f64>,
    pub peaks: StoredPeaks,
    /// The peptide validation assigned, when the spectrum was validated.
    pub assigned: Option<AssignedPeptide>,
    pub psms: Vec<EnginePsm>,
}

fn thresholds(conn: &Connection) -> Result<Thresholds, StoreError> {
    Ok(validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default())
}

fn count(value: i64) -> usize {
    usize::try_from(value).unwrap_or_default()
}

/// Every biosample, by name.
pub(crate) fn biosamples(conn: &Connection) -> Result<Vec<Biosample>, StoreError> {
    conn.prepare("SELECT id, name, category FROM biosamples ORDER BY name")?
        .query_map([], |row| {
            Ok(Biosample {
                id: row.get(0)?,
                name: row.get(1)?,
                category: row.get(2)?,
            })
        })?
        .collect::<Result<_, _>>()
        .map_err(Into::into)
}

/// One group's counts per biosample (see [`ProteinListRow`]).
struct SampleCounts {
    spectra: Vec<usize>,
    exclusive_spectra: Vec<usize>,
    exclusive_peptides: Vec<usize>,
    peptides: Vec<usize>,
}

impl SampleCounts {
    fn new(samples: usize) -> Self {
        Self {
            spectra: vec![0; samples],
            exclusive_spectra: vec![0; samples],
            exclusive_peptides: vec![0; samples],
            peptides: vec![0; samples],
        }
    }
}

pub(crate) fn protein_list(conn: &Connection) -> Result<ProteinList, StoreError> {
    let thresholds = thresholds(conn)?;
    let measure = validation::measure_column(thresholds.mode);
    let biosamples = biosamples(conn)?;
    let column: HashMap<i64, usize> = biosamples
        .iter()
        .enumerate()
        .map(|(i, b)| (b.id, i))
        .collect();
    let mut by_sample: HashMap<i64, SampleCounts> = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT gp.group_id, r.biosample_id, COUNT(*), SUM(gp.role = 'unique'),
                COUNT(DISTINCT CASE WHEN gp.role = 'unique' THEN pe.sequence END),
                COUNT(DISTINCT pe.sequence)
         FROM protein_groups g
         JOIN protein_group_peptides gp ON gp.group_id = g.id
         JOIN peptides pe ON pe.id = gp.peptide_id
         JOIN spectrum_results sr ON sr.peptide_id = gp.peptide_id
         JOIN spectra s ON s.id = sr.spectrum_id
         JOIN ms_runs r ON r.id = s.ms_run_id
         WHERE g.is_decoy = 0 AND g.{measure} <= ?2 AND g.peptides >= ?3
           AND gp.role IN ('unique', 'razor') AND sr.is_decoy = 0 AND sr.{measure} <= ?1
         GROUP BY gp.group_id, r.biosample_id"
    ))?;
    let mut rows = statement.query(params![
        thresholds.bound(Level::Psm),
        thresholds.bound(Level::Protein),
        thresholds.min_peptides
    ])?;
    while let Some(row) = rows.next()? {
        let counts = by_sample
            .entry(row.get(0)?)
            .or_insert_with(|| SampleCounts::new(biosamples.len()));
        if let Some(&index) = column.get(&row.get::<_, i64>(1)?) {
            counts.spectra[index] = count(row.get(2)?);
            counts.exclusive_spectra[index] = count(row.get(3)?);
            counts.exclusive_peptides[index] = count(row.get(4)?);
            counts.peptides[index] = count(row.get(5)?);
        }
    }
    let peptides = group_sequences(conn, &thresholds)?;
    let mut modifications: HashMap<i64, Vec<String>> = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT DISTINCT gp.group_id,
                COALESCE(m.name, printf('%+.3f', m.mass_delta))
         FROM protein_groups g
         JOIN protein_group_peptides gp ON gp.group_id = g.id
         JOIN modifications m ON m.peptide_id = gp.peptide_id
         WHERE g.is_decoy = 0 AND g.{measure} <= ?1 AND g.peptides >= ?2
           AND gp.role IN ('unique', 'razor')
         ORDER BY 1, 2"
    ))?;
    let mut rows = statement.query(params![
        thresholds.bound(Level::Protein),
        thresholds.min_peptides
    ])?;
    while let Some(row) = rows.next()? {
        modifications
            .entry(row.get(0)?)
            .or_default()
            .push(row.get(1)?);
    }
    let mut go: HashMap<String, Vec<String>> = HashMap::new();
    let mut statement = conn.prepare(
        "SELECT a.uniprot, t.name || ' (' || t.aspect || ')'
         FROM go_annotations a JOIN go_terms t ON t.id = a.go_id
         ORDER BY 1, 2",
    )?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        go.entry(row.get(0)?).or_default().push(row.get(1)?);
    }
    let starred: std::collections::HashSet<String> = conn
        .prepare("SELECT accession FROM starred_proteins")?
        .query_map([], |row| row.get(0))?
        .collect::<Result<_, _>>()?;
    let mut list = conn
        .prepare(&format!(
            "SELECT g.id, p.accession, p.description, p.sequence, g.peptides, g.q_value, g.score,
                    (SELECT COUNT(*) FROM protein_group_members m
                      WHERE m.group_id = g.id AND m.role = 'equivalent'),
                    g.pep
             FROM protein_groups g JOIN proteins p ON p.id = g.leader_protein_id
             WHERE g.is_decoy = 0 AND g.{measure} <= ?1 AND g.peptides >= ?2"
        ))?
        .query_map(
            params![thresholds.bound(Level::Protein), thresholds.min_peptides],
            |row| {
                let group_id: i64 = row.get(0)?;
                let sequence: Option<String> = row.get(3)?;
                Ok((
                    group_id,
                    row.get::<_, String>(1)?,
                    row.get(2)?,
                    sequence,
                    row.get::<_, i64>(4)?,
                    row.get(5)?,
                    row.get(6)?,
                    row.get::<_, i64>(7)?,
                    row.get::<_, Option<f64>>(8)?,
                ))
            },
        )?
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .map(
            |(
                group_id,
                accession,
                description,
                sequence,
                distinct,
                q_value,
                score,
                equivalents,
                pep,
            )| {
                let counts = by_sample
                    .remove(&group_id)
                    .unwrap_or_else(|| SampleCounts::new(biosamples.len()));
                let length = sequence.as_deref().map(|s| s.chars().count());
                let observable_peptides = sequence.as_deref().map(observable_peptides);
                let coverage = sequence.map(|sequence| {
                    let spans = map_peptide_coverage(
                        &sequence,
                        peptides
                            .get(&group_id)
                            .into_iter()
                            .flatten()
                            .map(|p| ("", p.as_str())),
                    );
                    covered_fraction(sequence.chars().count(), &spans)
                });
                let is_starred = starred.contains(&accession);
                let go_terms = tuffold_core::uniprot::uniprot_accession(&accession)
                    .and_then(|uniprot| go.get(uniprot).cloned())
                    .unwrap_or_default();
                ProteinListRow {
                    group_id,
                    accession,
                    description,
                    equivalents: count(equivalents),
                    peptides: count(distinct),
                    spectra: counts.spectra.iter().sum(),
                    spectra_by_sample: counts.spectra,
                    exclusive_spectra_by_sample: counts.exclusive_spectra,
                    exclusive_peptides_by_sample: counts.exclusive_peptides,
                    peptides_by_sample: counts.peptides,
                    q_value,
                    probability: pep.map(|pep| 1.0 - pep),
                    score,
                    coverage,
                    starred: is_starred,
                    go: go_terms,
                    length,
                    observable_peptides,
                    modifications: modifications.remove(&group_id).unwrap_or_default(),
                }
            },
        )
        .collect::<Vec<_>>();
    list.sort_by(|a, b| {
        b.spectra
            .cmp(&a.spectra)
            .then_with(|| a.accession.cmp(&b.accession))
    });
    Ok(ProteinList {
        biosamples,
        rows: list,
        thresholds,
    })
}

/// Passing target group id -> the distinct sequences of all its peptides.
fn group_sequences(
    conn: &Connection,
    thresholds: &Thresholds,
) -> Result<HashMap<i64, Vec<String>>, StoreError> {
    let mut groups: HashMap<i64, Vec<String>> = HashMap::new();
    let measure = validation::measure_column(thresholds.mode);
    let mut statement = conn.prepare(&format!(
        "SELECT DISTINCT gp.group_id, pe.sequence
         FROM protein_groups g
         JOIN protein_group_peptides gp ON gp.group_id = g.id
         JOIN peptides pe ON pe.id = gp.peptide_id
         WHERE g.is_decoy = 0 AND g.{measure} <= ?1 AND g.peptides >= ?2"
    ))?;
    let mut rows = statement.query(params![
        thresholds.bound(Level::Protein),
        thresholds.min_peptides
    ])?;
    while let Some(row) = rows.next()? {
        groups.entry(row.get(0)?).or_default().push(row.get(1)?);
    }
    Ok(groups)
}

pub(crate) fn group_detail(
    conn: &Connection,
    group_id: i64,
) -> Result<Option<GroupDetail>, StoreError> {
    let thresholds = thresholds(conn)?;
    let sequence: Option<Option<String>> = conn
        .query_row(
            "SELECT p.sequence FROM protein_groups g JOIN proteins p ON p.id = g.leader_protein_id
             WHERE g.id = ?1",
            [group_id],
            |row| row.get(0),
        )
        .optional()?;
    let Some(sequence) = sequence else {
        return Ok(None);
    };
    let members = conn
        .prepare(
            "SELECT p.id, p.accession, p.description, m.role
             FROM protein_group_members m JOIN proteins p ON p.id = m.protein_id
             WHERE m.group_id = ?1
             ORDER BY CASE m.role WHEN 'leader' THEN 0 WHEN 'equivalent' THEN 1 ELSE 2 END,
                      p.accession",
        )?
        .query_map([group_id], |row| {
            Ok(GroupMember {
                protein_id: row.get(0)?,
                accession: row.get(1)?,
                description: row.get(2)?,
                role: row.get(3)?,
            })
        })?
        .collect::<Result<_, _>>()?;
    let mut peptides: Vec<GroupPeptide> = conn
        .prepare(&format!(
            "SELECT pe.id, pe.proforma, pe.sequence, gp.role, r.q_value, r.score,
                    (SELECT COUNT(*) FROM spectrum_results sr
                      WHERE sr.peptide_id = pe.id AND sr.is_decoy = 0 AND sr.{measure} <= ?2),
                    r.pep
             FROM protein_group_peptides gp
             JOIN peptides pe ON pe.id = gp.peptide_id
             LEFT JOIN peptide_results r ON r.peptide_id = pe.id
             WHERE gp.group_id = ?1
             ORDER BY r.score DESC, pe.proforma",
            measure = validation::measure_column(thresholds.mode)
        ))?
        .query_map(params![group_id, thresholds.bound(Level::Psm)], |row| {
            Ok(GroupPeptide {
                peptide_id: row.get(0)?,
                proforma: row.get(1)?,
                sequence: row.get(2)?,
                role: row.get(3)?,
                q_value: row.get(4)?,
                probability: row.get::<_, Option<f64>>(7)?.map(|pep| 1.0 - pep),
                score: row.get(5)?,
                spectra: count(row.get(6)?),
                localized_peptide: None,
                site_probabilities: None,
                best_ascore: None,
                min_site_probability: None,
                site_scores: Vec::new(),
            })
        })?
        .collect::<Result<_, _>>()?;
    for peptide in &mut peptides {
        if let Some(loc) = best_peptide_localization(conn, peptide.peptide_id, &thresholds)? {
            peptide.localized_peptide = Some(loc.localized_peptide);
            peptide.site_probabilities = Some(loc.site_probabilities);
            peptide.best_ascore = Some(loc.best_ascore);
            peptide.min_site_probability = Some(loc.min_probability);
            peptide.site_scores = loc.sites;
        }
    }
    let coverage = sequence.as_deref().map_or_else(Vec::new, |sequence| {
        let ids: Vec<String> = peptides.iter().map(|p| p.peptide_id.to_string()).collect();
        map_peptide_coverage(
            sequence,
            ids.iter()
                .zip(&peptides)
                .map(|(id, peptide)| (id.as_str(), peptide.sequence.as_str())),
        )
    });
    let similar = similar_groups(conn, &thresholds, group_id)?;
    let cluster_size = cluster_size(conn, &thresholds, group_id)?;
    Ok(Some(GroupDetail {
        group_id,
        members,
        sequence,
        peptides,
        coverage,
        similar,
        cluster_size,
    }))
}

/// The SQL condition for a passing target group (alias `g`), with the
/// protein bound as `?2` and the peptide minimum as `?3`.
fn passing_group(thresholds: &Thresholds) -> String {
    format!(
        "g.is_decoy = 0 AND g.{} <= ?2 AND g.peptides >= ?3",
        validation::measure_column(thresholds.mode)
    )
}

fn similar_groups(
    conn: &Connection,
    thresholds: &Thresholds,
    group_id: i64,
) -> Result<Vec<SimilarGroup>, StoreError> {
    let mut statement = conn.prepare(&format!(
        "SELECT other.group_id, p.accession, p.description, g.peptides, other.peptide_id
         FROM protein_group_peptides mine
         JOIN protein_group_peptides other
           ON other.peptide_id = mine.peptide_id AND other.group_id != mine.group_id
         JOIN protein_groups g ON g.id = other.group_id
         JOIN proteins p ON p.id = g.leader_protein_id
         WHERE mine.group_id = ?1 AND {}
         ORDER BY other.group_id, other.peptide_id",
        passing_group(thresholds)
    ))?;
    let mut rows = statement.query(params![
        group_id,
        thresholds.bound(Level::Protein),
        thresholds.min_peptides
    ])?;
    let mut similar: Vec<SimilarGroup> = Vec::new();
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        if similar.last().is_none_or(|last| last.group_id != id) {
            similar.push(SimilarGroup {
                group_id: id,
                accession: row.get(1)?,
                description: row.get(2)?,
                shared: Vec::new(),
                peptides: count(row.get(3)?),
            });
        }
        if let Some(last) = similar.last_mut() {
            last.shared.push(row.get(4)?);
        }
    }
    similar.sort_by(|a, b| {
        b.shared
            .len()
            .cmp(&a.shared.len())
            .then_with(|| a.accession.cmp(&b.accession))
    });
    Ok(similar)
}

/// How many passing groups the shared-peptide graph connects to `group_id`
/// (breadth first, from the group itself).
fn cluster_size(
    conn: &Connection,
    thresholds: &Thresholds,
    group_id: i64,
) -> Result<usize, StoreError> {
    let mut statement = conn.prepare(&format!(
        "SELECT DISTINCT other.group_id
         FROM protein_group_peptides mine
         JOIN protein_group_peptides other
           ON other.peptide_id = mine.peptide_id AND other.group_id != mine.group_id
         JOIN protein_groups g ON g.id = other.group_id
         WHERE mine.group_id = ?1 AND {}",
        passing_group(thresholds)
    ))?;
    let mut seen = std::collections::HashSet::from([group_id]);
    let mut queue = vec![group_id];
    while let Some(current) = queue.pop() {
        let neighbours = statement
            .query_map(
                params![
                    current,
                    thresholds.bound(Level::Protein),
                    thresholds.min_peptides
                ],
                |row| row.get::<_, i64>(0),
            )?
            .collect::<Result<Vec<_>, _>>()?;
        for neighbour in neighbours {
            if seen.insert(neighbour) {
                queue.push(neighbour);
            }
        }
    }
    Ok(seen.len())
}

pub(crate) fn peptide_spectra(
    conn: &Connection,
    peptide_id: i64,
    limit: usize,
) -> Result<SpectrumPage, StoreError> {
    let thresholds = thresholds(conn)?;
    let total: i64 = conn.query_row(
        "SELECT COUNT(*) FROM spectrum_results WHERE peptide_id = ?1",
        [peptide_id],
        |row| row.get(0),
    )?;
    let limit = i64::try_from(limit).unwrap_or(i64::MAX);
    let mut rows = conn
        .prepare(
            "SELECT s.id, r.name, b.name, s.spectrum_key, s.precursor_charge, s.precursor_mz,
                    s.rt_seconds, sr.q_value, sr.score, sr.engines, s.peaks IS NOT NULL, sr.pep
             FROM spectrum_results sr
             JOIN spectra s ON s.id = sr.spectrum_id
             JOIN ms_runs r ON r.id = s.ms_run_id
             JOIN biosamples b ON b.id = r.biosample_id
             WHERE sr.peptide_id = ?1
             ORDER BY sr.score DESC, s.id
             LIMIT ?2",
        )?
        .query_map(params![peptide_id, limit], |row| {
            let q_value: f64 = row.get(7)?;
            let pep: f64 = row.get(11)?;
            Ok(SpectrumListRow {
                spectrum_id: row.get(0)?,
                run: row.get(1)?,
                biosample: row.get(2)?,
                key: row.get(3)?,
                charge: row.get(4)?,
                precursor_mz: row.get(5)?,
                rt_seconds: row.get(6)?,
                q_value,
                probability: 1.0 - pep,
                score: row.get(8)?,
                engines: row.get(9)?,
                has_peaks: row.get(10)?,
                passes: thresholds.passes(Level::Psm, Some(q_value), Some(pep)),
                localized_peptide: None,
                site_probabilities: None,
                best_ascore: None,
                min_site_probability: None,
                site_scores: Vec::new(),
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    for row in &mut rows {
        if let Some(loc) = spectrum_localization(conn, row.spectrum_id)? {
            row.localized_peptide = Some(loc.localized_peptide);
            row.site_probabilities = Some(loc.site_probabilities);
            row.best_ascore = Some(loc.best_ascore);
            row.min_site_probability = Some(loc.min_probability);
            row.site_scores = loc.sites;
        }
    }
    Ok(SpectrumPage {
        rows,
        total: usize::try_from(total).unwrap_or(0),
    })
}

pub(crate) fn spectrum_detail(
    conn: &Connection,
    spectrum_id: i64,
) -> Result<Option<SpectrumDetail>, StoreError> {
    let found = conn
        .query_row(
            "SELECT r.name, s.spectrum_key, s.precursor_mz, s.precursor_charge, s.rt_seconds,
                    s.peaks
             FROM spectra s JOIN ms_runs r ON r.id = s.ms_run_id WHERE s.id = ?1",
            [spectrum_id],
            |row| {
                Ok((
                    row.get::<_, String>(0)?,
                    row.get::<_, String>(1)?,
                    row.get::<_, Option<f64>>(2)?,
                    row.get::<_, Option<i64>>(3)?,
                    row.get::<_, Option<f64>>(4)?,
                    row.get::<_, Option<Vec<u8>>>(5)?,
                ))
            },
        )
        .optional()?;
    let Some((run, key, precursor_mz, precursor_charge, rt_seconds, blob)) = found else {
        return Ok(None);
    };
    let peaks = match blob {
        None => StoredPeaks::Missing,
        Some(blob) => peaks::decode(&blob).map_or(StoredPeaks::Unreadable, StoredPeaks::Stored),
    };
    let psms = engine_psms(conn, spectrum_id)?;
    let assigned = conn
        .query_row(
            "SELECT peptide_id, is_decoy, q_value, score, pep FROM spectrum_results
             WHERE spectrum_id = ?1",
            [spectrum_id],
            |row| {
                Ok((
                    row.get::<_, i64>(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get::<_, f64>(4)?,
                ))
            },
        )
        .optional()?
        .map(
            |(peptide_id, is_decoy, q_value, score, pep)| -> Result<_, StoreError> {
                let form = peptide_form(conn, peptide_id)?;
                let psm_charge = psm_charge(conn, spectrum_id, peptide_id)?;
                let charge = psm_charge
                    .or(precursor_charge)
                    .and_then(|charge| u8::try_from(charge).ok());
                let loc = spectrum_localization(conn, spectrum_id)?;
                let (
                    localized_peptide,
                    site_probabilities,
                    best_ascore,
                    min_site_probability,
                    site_scores,
                ) = match loc {
                    Some(loc) => (
                        Some(loc.localized_peptide),
                        Some(loc.site_probabilities),
                        Some(loc.best_ascore),
                        Some(loc.min_probability),
                        loc.sites,
                    ),
                    None => (None, None, None, None, Vec::new()),
                };
                Ok(AssignedPeptide {
                    peptide_id,
                    form,
                    is_decoy,
                    q_value,
                    probability: 1.0 - pep,
                    score,
                    charge,
                    localized_peptide,
                    site_probabilities,
                    best_ascore,
                    min_site_probability,
                    site_scores,
                })
            },
        )
        .transpose()?;
    Ok(Some(SpectrumDetail {
        spectrum_id,
        run,
        key,
        precursor_mz,
        precursor_charge,
        rt_seconds,
        peaks,
        assigned,
        psms,
    }))
}

pub(crate) struct StoredLocalization {
    pub localized_peptide: String,
    pub best_ascore: f64,
    pub min_probability: f64,
    pub sites: Vec<SiteScore>,
    pub site_probabilities: String,
}

pub(crate) fn spectrum_localization(
    conn: &Connection,
    spectrum_id: i64,
) -> Result<Option<StoredLocalization>, StoreError> {
    let summary: Option<(Option<String>, Option<f64>, Option<f64>)> = conn
        .query_row(
            "SELECT localized_peptide, best_ascore, min_probability
             FROM ptm_spectrum_summary
             WHERE spectrum_id = ?1 AND status = 'scored'",
            [spectrum_id],
            |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
        )
        .optional()?;
    let Some((Some(localized_peptide), Some(best_ascore), Some(min_probability))) = summary else {
        return Ok(None);
    };
    let sites: Vec<SiteScore> = conn
        .prepare_cached(
            "SELECT position, residue, mass_delta, modification_name, probability, ascore,
                    determining_matched, determining_total, assigned
             FROM ptm_site_scores
             WHERE spectrum_id = ?1 AND is_decoy = 0
             ORDER BY position, mass_delta",
        )?
        .query_map([spectrum_id], |row| {
            let residue_str: String = row.get(1)?;
            Ok(SiteScore {
                position: row.get(0)?,
                residue: residue_str.chars().next().unwrap_or('?'),
                mass_delta: row.get(2)?,
                modification_name: row.get(3)?,
                probability: row.get(4)?,
                ascore: row.get(5)?,
                determining_ions_matched: row.get(6)?,
                determining_ions_total: row.get(7)?,
                assigned: row.get::<_, i64>(8)? != 0,
                is_decoy: false,
            })
        })?
        .collect::<Result<_, _>>()?;
    let site_probabilities = format_sites(&sites);
    Ok(Some(StoredLocalization {
        localized_peptide,
        best_ascore,
        min_probability,
        sites,
        site_probabilities,
    }))
}

/// A peptide's localization is that of its passing spectrum with the highest
/// minimum site probability (then Ascore, spectrum score, lower spectrum id):
/// the best spectrum's, not a combined estimate over its spectra. `report`
/// loads the same choice in bulk (`best_localized_spectra`); keep them equal.
pub(crate) fn best_peptide_localization(
    conn: &Connection,
    peptide_id: i64,
    thresholds: &Thresholds,
) -> Result<Option<StoredLocalization>, StoreError> {
    let measure = validation::measure_column(thresholds.mode);
    let best_spec: Option<i64> = conn
        .query_row(
            &format!(
                "SELECT ps.spectrum_id
                 FROM spectrum_results sr
                 JOIN ptm_spectrum_summary ps ON ps.spectrum_id = sr.spectrum_id
                 WHERE sr.peptide_id = ?1 AND sr.is_decoy = 0 AND sr.{measure} <= ?2
                   AND ps.status = 'scored'
                 ORDER BY ps.min_probability DESC, ps.best_ascore DESC, sr.score DESC, sr.spectrum_id
                 LIMIT 1"
            ),
            params![peptide_id, thresholds.bound(Level::Psm)],
            |row| row.get(0),
        )
        .optional()?;
    match best_spec {
        Some(spec_id) => spectrum_localization(conn, spec_id),
        None => Ok(None),
    }
}

fn engine_psms(conn: &Connection, spectrum_id: i64) -> Result<Vec<EnginePsm>, StoreError> {
    let mut psms: Vec<(i64, EnginePsm)> = conn
        .prepare(
            "SELECT p.id, i.software, pe.proforma, p.charge, r.q_value, r.pep
             FROM psms p
             JOIN imports i ON i.id = p.import_id
             JOIN peptides pe ON pe.id = p.peptide_id
             LEFT JOIN psm_results r ON r.psm_id = p.id
             WHERE p.spectrum_id = ?1 AND p.rank = 1
             ORDER BY i.software, p.id",
        )?
        .query_map([spectrum_id], |row| {
            Ok((
                row.get(0)?,
                EnginePsm {
                    engine: row.get(1)?,
                    proforma: row.get(2)?,
                    charge: row.get(3)?,
                    scores: Vec::new(),
                    q_value: row.get(4)?,
                    pep: row.get(5)?,
                },
            ))
        })?
        .collect::<Result<_, _>>()?;
    let mut statement =
        conn.prepare_cached("SELECT name, value FROM psm_scores WHERE psm_id = ?1 ORDER BY name")?;
    for (psm_id, psm) in &mut psms {
        psm.scores = statement
            .query_map([*psm_id], |row| Ok((row.get(0)?, row.get(1)?)))?
            .collect::<Result<_, _>>()?;
    }
    Ok(psms.into_iter().map(|(_, psm)| psm).collect())
}

/// The charge of the rank-1 PSM that reported this peptide for the spectrum.
fn psm_charge(
    conn: &Connection,
    spectrum_id: i64,
    peptide_id: i64,
) -> Result<Option<i64>, StoreError> {
    Ok(conn
        .query_row(
            "SELECT charge FROM psms WHERE spectrum_id = ?1 AND peptide_id = ?2 AND rank = 1
             ORDER BY id LIMIT 1",
            [spectrum_id, peptide_id],
            |row| row.get(0),
        )
        .optional()?)
}

/// Rebuild a peptide form from its stored sequence and modifications.
pub(crate) fn peptide_form(conn: &Connection, peptide_id: i64) -> Result<PeptideForm, StoreError> {
    let sequence: String = conn.query_row(
        "SELECT sequence FROM peptides WHERE id = ?1",
        [peptide_id],
        |row| row.get(0),
    )?;
    let modifications = conn
        .prepare(
            "SELECT site_kind, position, mass_delta, name, accession FROM modifications
             WHERE peptide_id = ?1 ORDER BY ordinal",
        )?
        .query_map([peptide_id], |row| {
            let kind: String = row.get(0)?;
            let position: Option<u32> = row.get(1)?;
            let site = match (kind.as_str(), position) {
                ("n_term", _) => ModSite::NTerm,
                ("c_term", _) => ModSite::CTerm,
                (_, Some(position)) => ModSite::Residue(position),
                (_, None) => ModSite::Residue(0),
            };
            Ok(Modification {
                site,
                mass_delta: row.get(2)?,
                name: row.get(3)?,
                accession: row.get(4)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    PeptideForm::new(&sequence, modifications)
        .map_err(|error| StoreError::Damaged(format!("peptide {peptide_id}: {error}")))
}
