//! Reading validation inputs and writing validation results.
//!
//! Results replace one another wholesale: a validation run rewrites every
//! result table, a threshold change rewrites the protein tables, and an
//! import clears them all (the settings row survives, with `validated_at`
//! set to NULL).

use std::collections::HashMap;

use rusqlite::{Connection, OptionalExtension, Transaction, params};
use tuffold_validate::{
    DEFAULT_SITE_PROBABILITY, Input, Level, LocalizationSummary, MemberRole, Method, Note,
    NoteCode, PeptideRole, PeptideRow, ProteinInference, ProteinRow, PsmRow, SIBLINGS,
    SIBLINGS_IMPUTED, Settings, SpectrumLocalization, ThresholdMode, Thresholds, Validation,
    estimate_flr,
};

use crate::StoreError;

/// Protein inference results, children first (foreign keys).
const PROTEIN_TABLES: &[&str] = &[
    "protein_group_peptides",
    "protein_group_members",
    "protein_groups",
];
/// Every other result table.
const PSM_TABLES: &[&str] = &[
    "ptm_site_scores",
    "ptm_spectrum_summary",
    "peptide_results",
    "spectrum_results",
    "psm_results",
    "validation_notes",
    "validation_engines",
];

/// How one engine was ranked, as stored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EngineRecord {
    pub name: String,
    pub psms: usize,
    pub decoys: usize,
    /// `rescored`, `best_score`, or `unscored`.
    pub method: String,
    /// A sentence for the user (what was used, or why not).
    pub detail: String,
    pub features: Vec<String>,
    /// Targets passing 1% FDR on this engine alone.
    pub accepted: usize,
}

/// How peptides were scored, as stored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PeptideScoring {
    /// `pooled` (over the peptide's spectra) or `best_spectrum`.
    pub method: String,
    /// A sentence for the user (what was used, or why not).
    pub detail: String,
    /// The features pooling combined (empty unless pooled).
    pub features: Vec<String>,
}

/// Targets and decoys passing a level's threshold.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LevelCounts {
    pub targets: usize,
    pub decoys: usize,
    /// The lowest q-value any target reaches. Target-decoy estimates cannot
    /// go below `1 / targets`, so a threshold under this passes nothing,
    /// however good the data; front ends say so.
    pub lowest_q: Option<f64>,
}

impl LevelCounts {
    /// The decoy estimate of the FDR of what passes, (decoys + 1) / targets
    /// as for q-values, at most 1. `None` when no target passes.
    pub fn estimated_fdr(&self) -> Option<f64> {
        (self.targets > 0).then(|| ((self.decoys + 1) as f64 / self.targets as f64).min(1.0))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidationSummary {
    pub settings: Settings,
    /// `None` while results are stale (data imported since) or never made.
    pub validated_at: Option<i64>,
    pub engines: Vec<EngineRecord>,
    /// How peptides were scored; `None` when no engine was ranked (or the
    /// file was validated before this was recorded).
    pub peptide_scoring: Option<PeptideScoring>,
    pub notes: Vec<Note>,
    pub spectra: LevelCounts,
    pub peptides: LevelCounts,
    pub protein_groups: LevelCounts,
    pub localization: LocalizationSummary,
}

/// A validated peptide, for protein inference at a new threshold.
#[derive(Debug, Clone, PartialEq)]
pub struct StoredPeptide {
    pub id: i64,
    pub sequence: String,
    pub is_decoy: bool,
    /// The score its protein group ranks by (`group_score`, else `score`).
    pub score: f64,
    pub q_value: f64,
    pub pep: f64,
}

/// Everything protein inference needs, read back from validation results.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct InferenceRows {
    pub peptides: Vec<StoredPeptide>,
    pub peptide_proteins: Vec<(i64, i64)>,
    pub proteins: Vec<ProteinRow>,
}

pub(crate) fn input(conn: &Connection) -> Result<Input, StoreError> {
    let mut input = Input::default();
    let mut engines: HashMap<String, usize> = HashMap::new();
    let mut statement = conn.prepare(
        "SELECT p.id, p.spectrum_id, p.peptide_id, i.software, p.charge, p.calculated_mz,
                p.experimental_mz, s.ms_run_id, s.rt_seconds
         FROM psms p JOIN imports i ON i.id = p.import_id JOIN spectra s ON s.id = p.spectrum_id
         WHERE p.rank = 1 ORDER BY p.id",
    )?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let software = text(row, 3)?;
        let engine = match engines.get(software) {
            Some(&engine) => engine,
            None => {
                engines.insert(software.to_string(), input.engines.len());
                input.engines.push(software.to_string());
                input.engines.len() - 1
            }
        };
        input.psms.push(PsmRow {
            id: row.get(0)?,
            spectrum_id: row.get(1)?,
            peptide_id: row.get(2)?,
            engine,
            charge: row.get(4)?,
            calculated_mz: row.get(5)?,
            experimental_mz: row.get(6)?,
            run_id: row.get(7)?,
            retention_time: row.get(8)?,
            scores: Vec::new(),
        });
    }
    attach_scores(conn, &mut input)?;
    input.peptides = conn
        .prepare(
            "SELECT pe.id, pe.sequence,
                    (SELECT COUNT(*) FROM modifications m WHERE m.peptide_id = pe.id),
                    d.is_decoy
             FROM peptides pe JOIN peptide_decoys d ON d.peptide_id = pe.id
             WHERE pe.id IN (SELECT peptide_id FROM psms WHERE rank = 1)
             ORDER BY pe.id",
        )?
        .query_map([], |row| {
            Ok(PeptideRow {
                id: row.get(0)?,
                sequence: row.get(1)?,
                modifications: row.get(2)?,
                is_decoy: row.get(3)?,
            })
        })?
        .collect::<Result<_, _>>()?;
    input.proteins = proteins(conn)?;
    input.peptide_proteins = peptide_proteins(conn, "SELECT peptide_id FROM psms WHERE rank = 1")?;
    Ok(input)
}

/// Rank-1 scores, merged into the PSMs they belong to (both ordered by id).
fn attach_scores(conn: &Connection, input: &mut Input) -> Result<(), StoreError> {
    let position: HashMap<i64, usize> = input
        .psms
        .iter()
        .enumerate()
        .map(|(index, psm)| (psm.id, index))
        .collect();
    let mut names: HashMap<String, usize> = HashMap::new();
    let mut statement = conn.prepare(
        "SELECT s.psm_id, s.name, s.value
         FROM psm_scores s JOIN psms p ON p.id = s.psm_id
         WHERE p.rank = 1 ORDER BY s.psm_id",
    )?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let psm_id: i64 = row.get(0)?;
        let name = text(row, 1)?;
        let index = match names.get(name) {
            Some(&index) => index,
            None => {
                names.insert(name.to_string(), input.score_names.len());
                input.score_names.push(name.to_string());
                input.score_names.len() - 1
            }
        };
        if let Some(&psm) = position.get(&psm_id) {
            input.psms[psm].scores.push((index, row.get(2)?));
        }
    }
    Ok(())
}

/// A borrowed text column, without allocating.
fn text<'r>(row: &'r rusqlite::Row<'_>, index: usize) -> Result<&'r str, rusqlite::Error> {
    row.get_ref(index)?.as_str().map_err(|error| {
        rusqlite::Error::FromSqlConversionFailure(
            index,
            rusqlite::types::Type::Text,
            Box::new(error),
        )
    })
}

fn proteins(conn: &Connection) -> Result<Vec<ProteinRow>, StoreError> {
    Ok(conn
        .prepare("SELECT id, accession, is_decoy FROM proteins ORDER BY id")?
        .query_map([], |row| {
            Ok(ProteinRow {
                id: row.get(0)?,
                accession: row.get(1)?,
                is_decoy: row.get(2)?,
            })
        })?
        .collect::<Result<_, _>>()?)
}

/// Distinct evidence pairs (from PSMs of any rank) for the peptides that
/// `peptides` (a SQL subquery of ids) selects.
fn peptide_proteins(conn: &Connection, peptides: &str) -> Result<Vec<(i64, i64)>, StoreError> {
    Ok(conn
        .prepare(&format!(
            "SELECT DISTINCT p.peptide_id, pp.protein_id
             FROM psm_proteins pp JOIN psms p ON p.id = pp.psm_id
             WHERE p.peptide_id IN ({peptides})
             ORDER BY 1, 2"
        ))?
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .collect::<Result<_, _>>()?)
}

pub(crate) fn inference_rows(conn: &Connection) -> Result<InferenceRows, StoreError> {
    let peptides = conn
        .prepare(
            "SELECT r.peptide_id, pe.sequence, r.is_decoy, COALESCE(r.group_score, r.score),
                    r.q_value, r.pep
             FROM peptide_results r JOIN peptides pe ON pe.id = r.peptide_id
             ORDER BY r.peptide_id",
        )?
        .query_map([], |row| {
            Ok(StoredPeptide {
                id: row.get(0)?,
                sequence: row.get(1)?,
                is_decoy: row.get(2)?,
                score: row.get(3)?,
                q_value: row.get(4)?,
                pep: row.get(5)?,
            })
        })?
        .collect::<Result<_, _>>()?;
    Ok(InferenceRows {
        peptides,
        peptide_proteins: peptide_proteins(conn, "SELECT peptide_id FROM peptide_results")?,
        proteins: proteins(conn)?,
    })
}

/// Clears every result and marks the settings stale (after an import).
pub(crate) fn invalidate(tx: &Transaction<'_>) -> Result<(), StoreError> {
    for table in PROTEIN_TABLES.iter().chain(PSM_TABLES) {
        tx.execute(&format!("DELETE FROM {table}"), [])?;
    }
    // Transfers rest on which precursors passed and where; new data
    // changes both.
    tx.execute("DELETE FROM transfers", [])?;
    tx.execute("DELETE FROM transfer_runs", [])?;
    tx.execute(
        "UPDATE validation
         SET validated_at = NULL, peptide_method = NULL, peptide_detail = NULL,
             peptide_features = NULL",
        [],
    )?;
    Ok(())
}

pub(crate) fn save(
    tx: &Transaction<'_>,
    settings: &Settings,
    validation: &Validation,
    validated_at: i64,
) -> Result<(), StoreError> {
    invalidate(tx)?;
    save_settings(tx, settings)?;
    tx.execute(
        "UPDATE validation SET validated_at = ?1, validated_with = ?2",
        params![validated_at, concat!("tuffold ", env!("CARGO_PKG_VERSION"))],
    )?;
    if let Some(method) = &validation.peptide_method {
        let scoring = describe_peptides(method);
        tx.execute(
            "UPDATE validation
             SET peptide_method = ?1, peptide_detail = ?2, peptide_features = ?3",
            params![scoring.method, scoring.detail, scoring.features.join("\n")],
        )?;
    }
    for (ordinal, engine) in (0_i64..).zip(&validation.engines) {
        let (method, detail, features) = describe(&engine.method);
        tx.execute(
            "INSERT INTO validation_engines
                 (ordinal, name, psms, decoys, method, detail, features, accepted)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
            params![
                ordinal,
                engine.name,
                count(engine.psms),
                count(engine.decoys),
                method,
                detail,
                features.join("\n"),
                count(engine.accepted)
            ],
        )?;
    }
    save_notes(tx, "psms", &validation.notes)?;
    let mut statement = tx
        .prepare("INSERT INTO psm_results (psm_id, score, q_value, pep) VALUES (?1, ?2, ?3, ?4)")?;
    for result in &validation.psms {
        statement.execute(params![
            result.psm_id,
            result.score,
            result.q_value,
            result.pep
        ])?;
    }
    let mut statement = tx.prepare(
        "INSERT INTO spectrum_results
             (spectrum_id, peptide_id, is_decoy, score, q_value, pep, engines)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
    )?;
    for found in &validation.spectra {
        statement.execute(params![
            found.spectrum_id,
            found.peptide_id,
            found.is_decoy,
            found.score,
            found.q_value,
            found.pep,
            found.engines
        ])?;
    }
    let mut statement = tx.prepare(
        "INSERT INTO peptide_results
             (peptide_id, is_decoy, score, q_value, pep, best_spectrum_id, group_score)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
    )?;
    for peptide in &validation.peptides {
        statement.execute(params![
            peptide.peptide_id,
            peptide.is_decoy,
            peptide.score,
            peptide.q_value,
            peptide.pep,
            peptide.best_spectrum_id,
            peptide.group_score
        ])?;
    }
    save_localizations(tx, &validation.localizations)?;
    save_proteins(tx, &validation.proteins)
}

fn save_localizations(
    tx: &Transaction<'_>,
    localizations: &[(i64, SpectrumLocalization)],
) -> Result<(), StoreError> {
    let mut summary_stmt = tx.prepare(
        "INSERT INTO ptm_spectrum_summary
             (spectrum_id, status, isoform_count, localized_peptide, best_ascore,
              min_probability, false_target_sites)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
    )?;
    let mut site_stmt = tx.prepare(
        "INSERT INTO ptm_site_scores
             (spectrum_id, is_decoy, position, residue, mass_delta, modification_name,
              probability, ascore, determining_matched, determining_total, assigned)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
    )?;
    for (spectrum_id, outcome) in localizations {
        match outcome {
            SpectrumLocalization::Unambiguous => {}
            SpectrumLocalization::MissingPeaks { isoform_count } => {
                summary_stmt.execute(params![
                    spectrum_id,
                    "missing_peaks",
                    count(*isoform_count),
                    Option::<String>::None,
                    Option::<f64>::None,
                    Option::<f64>::None,
                    0i64,
                ])?;
            }
            SpectrumLocalization::Unscorable { isoform_count } => {
                summary_stmt.execute(params![
                    spectrum_id,
                    "unscorable",
                    count(*isoform_count),
                    Option::<String>::None,
                    Option::<f64>::None,
                    Option::<f64>::None,
                    0i64,
                ])?;
            }
            SpectrumLocalization::Capped { isoform_count } => {
                summary_stmt.execute(params![
                    spectrum_id,
                    "capped",
                    count(*isoform_count),
                    Option::<String>::None,
                    Option::<f64>::None,
                    Option::<f64>::None,
                    0i64,
                ])?;
            }
            SpectrumLocalization::Scored(scored) => {
                summary_stmt.execute(params![
                    spectrum_id,
                    "scored",
                    count(scored.isoform_count),
                    &scored.localized_peptide,
                    scored.best_ascore,
                    scored.min_probability,
                    i64::from(scored.false_target_sites),
                ])?;
                for site in &scored.sites {
                    site_stmt.execute(params![
                        spectrum_id,
                        i64::from(site.is_decoy),
                        i64::from(site.position),
                        site.residue.to_string(),
                        site.mass_delta,
                        &site.modification_name,
                        site.probability,
                        site.ascore,
                        i64::from(site.determining_ions_matched),
                        i64::from(site.determining_ions_total),
                        i64::from(site.assigned),
                    ])?;
                }
            }
        }
    }
    Ok(())
}

/// Replaces protein groups (and their notes) after a threshold change.
pub(crate) fn save_inference(
    tx: &Transaction<'_>,
    thresholds: &Thresholds,
    inference: &ProteinInference,
) -> Result<(), StoreError> {
    for table in PROTEIN_TABLES {
        tx.execute(&format!("DELETE FROM {table}"), [])?;
    }
    tx.execute("DELETE FROM validation_notes WHERE stage = 'proteins'", [])?;
    save_thresholds(tx, thresholds)?;
    save_proteins(tx, inference)
}

pub(crate) fn save_thresholds(
    tx: &Transaction<'_>,
    thresholds: &Thresholds,
) -> Result<(), StoreError> {
    tx.execute(
        "UPDATE validation SET psm_fdr = ?1, peptide_fdr = ?2, protein_fdr = ?3,
                               min_peptides = ?4, threshold_mode = ?5",
        params![
            thresholds.psm,
            thresholds.peptide,
            thresholds.protein,
            thresholds.min_peptides,
            thresholds.mode.as_str()
        ],
    )?;
    Ok(())
}

fn save_settings(tx: &Transaction<'_>, settings: &Settings) -> Result<(), StoreError> {
    let thresholds = &settings.thresholds;
    tx.execute(
        "INSERT INTO validation
             (id, rescore, psm_fdr, peptide_fdr, protein_fdr, min_peptides, threshold_mode)
         VALUES (1, ?1, ?2, ?3, ?4, ?5, ?6)
         ON CONFLICT (id) DO UPDATE SET rescore = excluded.rescore,
             psm_fdr = excluded.psm_fdr, peptide_fdr = excluded.peptide_fdr,
             protein_fdr = excluded.protein_fdr, min_peptides = excluded.min_peptides,
             threshold_mode = excluded.threshold_mode",
        params![
            settings.rescore,
            thresholds.psm,
            thresholds.peptide,
            thresholds.protein,
            thresholds.min_peptides,
            thresholds.mode.as_str()
        ],
    )?;
    Ok(())
}

fn save_notes(tx: &Transaction<'_>, stage: &str, notes: &[Note]) -> Result<(), StoreError> {
    for note in notes {
        tx.execute(
            "INSERT INTO validation_notes (stage, code, message) VALUES (?1, ?2, ?3)",
            params![stage, note.code.as_str(), note.message],
        )?;
    }
    Ok(())
}

fn save_proteins(tx: &Transaction<'_>, inference: &ProteinInference) -> Result<(), StoreError> {
    save_notes(tx, "proteins", &inference.notes)?;
    let mut group_statement = tx.prepare(
        "INSERT INTO protein_groups (leader_protein_id, is_decoy, score, q_value, peptides, pep)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
    )?;
    let mut member_statement = tx.prepare(
        "INSERT INTO protein_group_members (group_id, protein_id, role) VALUES (?1, ?2, ?3)
         ON CONFLICT DO NOTHING",
    )?;
    let mut peptide_statement = tx.prepare(
        "INSERT INTO protein_group_peptides (group_id, peptide_id, role) VALUES (?1, ?2, ?3)",
    )?;
    for group in &inference.groups {
        group_statement.execute(params![
            group.leader,
            group.is_decoy,
            group.score,
            group.q_value,
            group.distinct_sequences,
            group.pep
        ])?;
        let id = tx.last_insert_rowid();
        for (protein, role) in &group.members {
            let role = match role {
                MemberRole::Leader => "leader",
                MemberRole::Equivalent => "equivalent",
                MemberRole::Subset => "subset",
            };
            member_statement.execute(params![id, protein, role])?;
        }
        for (peptide, role) in &group.peptides {
            let role = match role {
                PeptideRole::Unique => "unique",
                PeptideRole::Razor => "razor",
                PeptideRole::Shared => "shared",
            };
            peptide_statement.execute(params![id, peptide, role])?;
        }
    }
    Ok(())
}

/// `(method, detail, features)` for storage.
fn describe(method: &Method) -> (&'static str, String, Vec<String>) {
    match method {
        Method::Rescored { features, model } => (
            "rescored",
            format!(
                "Rescored by cross-validated {} over {} features.",
                model.describe(),
                features.len()
            ),
            features.clone(),
        ),
        Method::BestScore {
            score,
            higher_is_better,
            why,
        } => {
            let direction = if *higher_is_better { "higher" } else { "lower" };
            (
                "best_score",
                format!(
                    "Ranked by {score} ({direction} is better): {}.",
                    why.describe()
                ),
                vec![score.clone()],
            )
        }
        Method::Unscored { why } => (
            "unscored",
            format!("Not validated: {}.", why.describe()),
            Vec::new(),
        ),
    }
}

/// How peptides were scored, for storage. Peptide scoring never leaves a
/// peptide unscored: without pooling, its best spectrum's score stands.
fn describe_peptides(method: &Method) -> PeptideScoring {
    match method {
        Method::Rescored { model, features } => {
            let siblings = if features.iter().any(|feature| feature == SIBLINGS) {
                " and their proteins' other identified peptides"
            } else if features.iter().any(|feature| feature == SIBLINGS_IMPUTED) {
                " and their proteins' other identified peptides (decoys' counts imputed from \
                 low-scoring targets')"
            } else {
                ""
            };
            PeptideScoring {
                method: "pooled".to_string(),
                detail: format!(
                    "Pooled over each peptide's spectra (best and second-best spectrum scores, \
                     other spectra, other charges){siblings} by cross-validated {}.",
                    model.describe()
                ),
                features: features.clone(),
            }
        }
        Method::BestScore { why, .. } => PeptideScoring {
            method: "best_spectrum".to_string(),
            detail: format!(
                "Scored by their best spectrum, not pooled: {}.",
                why.describe()
            ),
            features: Vec::new(),
        },
        Method::Unscored { why } => PeptideScoring {
            method: "best_spectrum".to_string(),
            detail: format!("Scored by their best spectrum: {}.", why.describe()),
            features: Vec::new(),
        },
    }
}

fn count(value: usize) -> i64 {
    i64::try_from(value).unwrap_or(i64::MAX)
}

pub(crate) fn settings(conn: &Connection) -> Result<Option<(Settings, Option<i64>)>, StoreError> {
    Ok(conn
        .query_row(
            "SELECT rescore, psm_fdr, peptide_fdr, protein_fdr, min_peptides, validated_at,
                    threshold_mode
             FROM validation WHERE id = 1",
            [],
            |row| {
                let mode: String = row.get(6)?;
                let settings = Settings {
                    rescore: row.get(0)?,
                    thresholds: Thresholds {
                        // The column's CHECK allows only the two modes.
                        mode: ThresholdMode::parse(&mode).unwrap_or_default(),
                        psm: row.get(1)?,
                        peptide: row.get(2)?,
                        protein: row.get(3)?,
                        min_peptides: row.get(4)?,
                    },
                };
                Ok((settings, row.get(5)?))
            },
        )
        .optional()?)
}

pub(crate) fn summary(conn: &Connection) -> Result<Option<ValidationSummary>, StoreError> {
    let Some((settings, validated_at)) = settings(conn)? else {
        return Ok(None);
    };
    let thresholds = settings.thresholds;
    let engines = conn
        .prepare(
            "SELECT name, psms, decoys, method, detail, features, accepted
             FROM validation_engines ORDER BY ordinal",
        )?
        .query_map([], |row| {
            let features: String = row.get(5)?;
            Ok(EngineRecord {
                name: row.get(0)?,
                psms: usize_of(row.get(1)?),
                decoys: usize_of(row.get(2)?),
                method: row.get(3)?,
                detail: row.get(4)?,
                features: features
                    .lines()
                    .filter(|line| !line.is_empty())
                    .map(String::from)
                    .collect(),
                accepted: usize_of(row.get(6)?),
            })
        })?
        .collect::<Result<_, _>>()?;
    let notes = conn
        .prepare("SELECT code, message FROM validation_notes ORDER BY id")?
        .query_map([], |row| {
            Ok((row.get::<_, String>(0)?, row.get::<_, String>(1)?))
        })?
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter_map(|(code, message)| Some(Note::new(NoteCode::parse(&code)?, message)))
        .collect();
    let peptide_scoring = conn
        .query_row(
            "SELECT peptide_method, peptide_detail, peptide_features
             FROM validation WHERE id = 1",
            [],
            |row| {
                Ok((
                    row.get::<_, Option<String>>(0)?,
                    row.get::<_, Option<String>>(1)?,
                    row.get::<_, Option<String>>(2)?,
                ))
            },
        )
        .optional()?
        .and_then(|(method, detail, features)| {
            Some(PeptideScoring {
                method: method?,
                detail: detail.unwrap_or_default(),
                features: features
                    .unwrap_or_default()
                    .lines()
                    .filter(|line| !line.is_empty())
                    .map(String::from)
                    .collect(),
            })
        });
    // q-values in fdr mode, PEPs in probability mode (see Thresholds::bound).
    let measure = measure_column(thresholds.mode);
    let level = |table: &str, extra: &str, values: &[&dyn rusqlite::ToSql]| {
        let sql = format!(
            "SELECT COUNT(*) FILTER (WHERE is_decoy = 0 AND {measure} <= ?1 {extra}),
                    COUNT(*) FILTER (WHERE is_decoy = 1 AND {measure} <= ?1 {extra}),
                    MIN(q_value) FILTER (WHERE is_decoy = 0 {extra})
             FROM {table}"
        );
        conn.query_row(&sql, values, |row| {
            Ok(LevelCounts {
                targets: usize_of(row.get(0)?),
                decoys: usize_of(row.get(1)?),
                lowest_q: row.get(2)?,
            })
        })
    };
    let localization = localization_summary(conn, &thresholds, DEFAULT_SITE_PROBABILITY)?;
    Ok(Some(ValidationSummary {
        settings,
        validated_at,
        engines,
        peptide_scoring,
        notes,
        spectra: level("spectrum_results", "", &[&thresholds.bound(Level::Psm)])?,
        peptides: level("peptide_results", "", &[&thresholds.bound(Level::Peptide)])?,
        protein_groups: level(
            "protein_groups",
            "AND peptides >= ?2",
            &[&thresholds.bound(Level::Protein), &thresholds.min_peptides],
        )?,
        localization,
    }))
}

pub(crate) fn localization_summary(
    conn: &Connection,
    thresholds: &Thresholds,
    probability_threshold: f64,
) -> Result<LocalizationSummary, StoreError> {
    let measure = measure_column(thresholds.mode);
    let psm_bound = thresholds.bound(Level::Psm);

    let (
        scored_psms,
        localized_psms,
        missing_peaks_psms,
        capped_psms,
        unscorable_psms,
        false_target_sites,
    ) = conn.query_row(
        &format!(
            "SELECT COUNT(*) FILTER (WHERE ps.status = 'scored'),
                        COUNT(*) FILTER (WHERE ps.status = 'scored' AND ps.min_probability >= ?2),
                        COUNT(*) FILTER (WHERE ps.status = 'missing_peaks'),
                        COUNT(*) FILTER (WHERE ps.status = 'capped'),
                        COUNT(*) FILTER (WHERE ps.status = 'unscorable'),
                        COALESCE(SUM(ps.false_target_sites) FILTER (WHERE ps.status = 'scored'), 0)
                 FROM ptm_spectrum_summary ps
                 JOIN spectrum_results sr ON sr.spectrum_id = ps.spectrum_id
                 WHERE sr.is_decoy = 0 AND sr.{measure} <= ?1"
        ),
        params![psm_bound, probability_threshold],
        |row| {
            Ok((
                usize_of(row.get(0)?),
                usize_of(row.get(1)?),
                usize_of(row.get(2)?),
                usize_of(row.get(3)?),
                usize_of(row.get(4)?),
                usize_of(row.get(5)?),
            ))
        },
    )?;

    let (total_target_sites, localized_sites, decoy_sites_tested, passing_decoy_sites) = conn
        .query_row(
            &format!(
                "SELECT COUNT(*) FILTER (WHERE ss.is_decoy = 0),
                        COUNT(*) FILTER (WHERE ss.is_decoy = 0 AND ss.probability >= ?2),
                        COUNT(*) FILTER (WHERE ss.is_decoy = 1),
                        COUNT(*) FILTER (WHERE ss.is_decoy = 1 AND ss.probability >= ?2)
                 FROM ptm_site_scores ss
                 JOIN ptm_spectrum_summary ps ON ps.spectrum_id = ss.spectrum_id
                 JOIN spectrum_results sr ON sr.spectrum_id = ss.spectrum_id
                 WHERE sr.is_decoy = 0 AND sr.{measure} <= ?1 AND ps.status = 'scored'"
            ),
            params![psm_bound, probability_threshold],
            |row| {
                Ok((
                    usize_of(row.get(0)?),
                    usize_of(row.get(1)?),
                    usize_of(row.get(2)?),
                    usize_of(row.get(3)?),
                ))
            },
        )?;

    let estimated_flr = estimate_flr(
        false_target_sites,
        decoy_sites_tested,
        localized_sites,
        passing_decoy_sites,
    );

    Ok(LocalizationSummary {
        scored_psms,
        localized_psms,
        missing_peaks_psms,
        capped_psms,
        unscorable_psms,
        localized_sites,
        total_target_sites,
        probability_threshold,
        estimated_flr,
    })
}

/// The result column thresholds filter on in `mode`.
pub(crate) fn measure_column(mode: ThresholdMode) -> &'static str {
    match mode {
        ThresholdMode::Fdr => "q_value",
        ThresholdMode::Probability => "pep",
    }
}

fn usize_of(value: i64) -> usize {
    usize::try_from(value).unwrap_or_default()
}
