//! Writing an [`ImportBatch`] into the store.
//!
//! Merge rules, per natural key:
//! - biosamples, MS runs, spectra, peptides, proteins: inserted once; a later
//!   import never overwrites a recorded value, it only fills fields that were
//!   still unknown (e.g. a precursor m/z one engine omitted). A run keeps the
//!   biosample it was first imported under.
//! - protein decoy status: an explicit source flag (mzIdentML `isDecoy`) wins
//!   over one inferred from the accession; among equals, decoy wins.
//! - PSMs, scores, evidence: always inserted, attributed to their import. A
//!   PSM's decoy status is not stored: it follows from its proteins at query
//!   time, so later decoy information applies to earlier imports too.
//! - Repeated scores (same name) or evidence (same protein and start) on one
//!   PSM keep the first and are counted in the [`ImportReport`], never
//!   dropped silently.

use std::collections::HashMap;

use rusqlite::{OptionalExtension, Transaction, params};
use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{ImportBatch, ModSite, PeptideForm};

use crate::{ImportSource, StoreError};

/// What one import added.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ImportReport {
    pub import_id: i64,
    pub psms: usize,
    pub new_spectra: usize,
    pub new_peptides: usize,
    pub new_proteins: usize,
    pub new_ms_runs: usize,
    /// Scores whose name repeated on the same PSM (the first value is kept).
    pub duplicate_scores: usize,
    /// Evidence repeating a protein and start on the same PSM (kept once).
    pub duplicate_evidence: usize,
}

impl ImportReport {
    /// Add another import's counts (several engines from one source).
    pub fn absorb(&mut self, other: &ImportReport) {
        self.psms += other.psms;
        self.new_spectra += other.new_spectra;
        self.new_peptides += other.new_peptides;
        self.new_proteins += other.new_proteins;
        self.new_ms_runs += other.new_ms_runs;
        self.duplicate_scores += other.duplicate_scores;
        self.duplicate_evidence += other.duplicate_evidence;
    }
}

/// Ids resolved so far in this import, so each natural key costs one query.
#[derive(Default)]
struct Ids {
    runs: HashMap<String, i64>,
    spectra: HashMap<(i64, String), i64>,
    proteins: HashMap<String, i64>,
    peptides: HashMap<String, i64>,
}

pub(crate) fn write(
    tx: &Transaction<'_>,
    source: &ImportSource,
    batch: &ImportBatch,
    imported_at: i64,
) -> Result<ImportReport, StoreError> {
    tx.execute(
        "INSERT INTO imports (source_path, format, software, software_version, imported_at)
         VALUES (?1, ?2, ?3, ?4, ?5)",
        params![
            source.path,
            source.format.as_str(),
            batch.software.name,
            batch.software.version,
            imported_at
        ],
    )?;
    let mut report = ImportReport {
        import_id: tx.last_insert_rowid(),
        ..ImportReport::default()
    };
    let mut ids = Ids::default();

    for biosample in &batch.biosamples {
        tx.prepare_cached(
            "INSERT INTO biosamples (name, category) VALUES (?1, ?2)
             ON CONFLICT (name) DO UPDATE SET category = COALESCE(category, excluded.category)",
        )?
        .execute((&biosample.name, &biosample.category))?;
    }
    for run in &batch.ms_runs {
        let id = ensure_run(tx, &mut ids, &mut report, &run.name, &run.biosample)?;
        tx.prepare_cached(
            "UPDATE ms_runs SET spectra_file = COALESCE(spectra_file, ?2) WHERE id = ?1",
        )?
        .execute(params![id, run.spectra_file])?;
    }
    for spectrum in &batch.spectra {
        let run_id = ensure_run(
            tx,
            &mut ids,
            &mut report,
            &spectrum.ms_run,
            &spectrum.ms_run,
        )?;
        let spectrum_id =
            ensure_spectrum(tx, &mut ids, &mut report, run_id, &spectrum.key.canonical())?;
        tx.prepare_cached(
            "UPDATE spectra SET
                source_id = COALESCE(source_id, ?2),
                precursor_mz = COALESCE(precursor_mz, ?3),
                precursor_charge = COALESCE(precursor_charge, ?4),
                rt_seconds = COALESCE(rt_seconds, ?5)
             WHERE id = ?1",
        )?
        .execute(params![
            spectrum_id,
            spectrum.source_id,
            spectrum.precursor_mz,
            spectrum.precursor_charge,
            spectrum.rt_seconds
        ])?;
    }
    for protein in &batch.proteins {
        let id = ensure_protein(tx, &mut ids, &mut report, &protein.accession)?;
        tx.prepare_cached(
            "UPDATE proteins SET
                description = COALESCE(description, ?2),
                sequence = COALESCE(sequence, ?3),
                is_decoy = CASE
                    WHEN ?5 = 1 AND decoy_explicit = 0 THEN ?4
                    WHEN ?5 = 0 AND decoy_explicit = 1 THEN is_decoy
                    ELSE MAX(is_decoy, ?4)
                END,
                decoy_explicit = MAX(decoy_explicit, ?5)
             WHERE id = ?1",
        )?
        .execute(params![
            id,
            protein.description,
            protein.sequence,
            protein.is_decoy,
            protein.decoy_explicit
        ])?;
    }

    for psm in &batch.psms {
        let run_id = ensure_run(tx, &mut ids, &mut report, &psm.ms_run, &psm.ms_run)?;
        // A PSM without a declared spectrum keeps a bare spectrum record
        // rather than being dropped.
        let spectrum_id =
            ensure_spectrum(tx, &mut ids, &mut report, run_id, &psm.spectrum.canonical())?;
        let peptide_id = ensure_peptide(tx, &mut ids, &mut report, &psm.peptide)?;
        tx.prepare_cached(
            "INSERT INTO psms (import_id, spectrum_id, peptide_id, rank, charge, calculated_mz,
                               experimental_mz, pass_threshold)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)",
        )?
        .execute(params![
            report.import_id,
            spectrum_id,
            peptide_id,
            psm.rank,
            psm.charge,
            psm.calculated_mz,
            psm.experimental_mz,
            psm.pass_threshold,
        ])?;
        let psm_id = tx.last_insert_rowid();
        for score in &psm.scores {
            let inserted = tx
                .prepare_cached(
                    "INSERT INTO psm_scores (psm_id, name, accession, value) VALUES (?1, ?2, ?3, ?4)
                     ON CONFLICT (psm_id, name) DO NOTHING",
                )?
                .execute(params![psm_id, score.name, score.accession, score.value])?;
            report.duplicate_scores += usize::from(inserted == 0);
        }
        for evidence in &psm.evidence {
            let protein_id = ensure_protein(tx, &mut ids, &mut report, &evidence.accession)?;
            let inserted = tx
                .prepare_cached(
                    "INSERT INTO psm_proteins (psm_id, protein_id, start_pos, end_pos, pre, post)
                     VALUES (?1, ?2, ?3, ?4, ?5, ?6)
                     ON CONFLICT DO NOTHING",
                )?
                .execute(params![
                    psm_id,
                    protein_id,
                    evidence.start,
                    evidence.end,
                    evidence.pre.map(String::from),
                    evidence.post.map(String::from),
                ])?;
            report.duplicate_evidence += usize::from(inserted == 0);
        }
        report.psms += 1;
    }
    Ok(report)
}

/// The run's id, creating it (and its biosample) when new.
fn ensure_run(
    tx: &Transaction<'_>,
    ids: &mut Ids,
    report: &mut ImportReport,
    name: &str,
    biosample: &str,
) -> Result<i64, StoreError> {
    if let Some(&id) = ids.runs.get(name) {
        return Ok(id);
    }
    let existing: Option<i64> = tx
        .prepare_cached("SELECT id FROM ms_runs WHERE name = ?1")?
        .query_row([name], |row| row.get(0))
        .optional()?;
    let id = match existing {
        Some(id) => id,
        None => {
            tx.prepare_cached(
                "INSERT INTO biosamples (name) VALUES (?1) ON CONFLICT (name) DO NOTHING",
            )?
            .execute([biosample])?;
            tx.prepare_cached(
                "INSERT INTO ms_runs (biosample_id, name) SELECT id, ?2 FROM biosamples WHERE name = ?1",
            )?
            .execute([biosample, name])?;
            report.new_ms_runs += 1;
            tx.last_insert_rowid()
        }
    };
    ids.runs.insert(name.to_string(), id);
    Ok(id)
}

fn ensure_spectrum(
    tx: &Transaction<'_>,
    ids: &mut Ids,
    report: &mut ImportReport,
    run_id: i64,
    key: &str,
) -> Result<i64, StoreError> {
    let cache_key = (run_id, key.to_string());
    if let Some(&id) = ids.spectra.get(&cache_key) {
        return Ok(id);
    }
    let existing: Option<i64> = tx
        .prepare_cached("SELECT id FROM spectra WHERE ms_run_id = ?1 AND spectrum_key = ?2")?
        .query_row(params![run_id, key], |row| row.get(0))
        .optional()?;
    let id = match existing {
        Some(id) => id,
        None => {
            tx.prepare_cached("INSERT INTO spectra (ms_run_id, spectrum_key) VALUES (?1, ?2)")?
                .execute(params![run_id, key])?;
            report.new_spectra += 1;
            tx.last_insert_rowid()
        }
    };
    ids.spectra.insert(cache_key, id);
    Ok(id)
}

/// The protein's id; one first seen only through evidence gets its decoy
/// status from the accession convention.
fn ensure_protein(
    tx: &Transaction<'_>,
    ids: &mut Ids,
    report: &mut ImportReport,
    accession: &str,
) -> Result<i64, StoreError> {
    if let Some(&id) = ids.proteins.get(accession) {
        return Ok(id);
    }
    let existing: Option<i64> = tx
        .prepare_cached("SELECT id FROM proteins WHERE accession = ?1")?
        .query_row([accession], |row| row.get(0))
        .optional()?;
    let id = match existing {
        Some(id) => id,
        None => {
            tx.prepare_cached("INSERT INTO proteins (accession, is_decoy) VALUES (?1, ?2)")?
                .execute(params![accession, is_decoy_accession(accession)])?;
            report.new_proteins += 1;
            tx.last_insert_rowid()
        }
    };
    ids.proteins.insert(accession.to_string(), id);
    Ok(id)
}

fn ensure_peptide(
    tx: &Transaction<'_>,
    ids: &mut Ids,
    report: &mut ImportReport,
    form: &PeptideForm,
) -> Result<i64, StoreError> {
    let key = form.proforma();
    if let Some(&id) = ids.peptides.get(&key) {
        return Ok(id);
    }
    let existing: Option<i64> = tx
        .prepare_cached("SELECT id FROM peptides WHERE proforma = ?1")?
        .query_row([&key], |row| row.get(0))
        .optional()?;
    let id = match existing {
        Some(id) => id,
        None => {
            tx.prepare_cached("INSERT INTO peptides (proforma, sequence) VALUES (?1, ?2)")?
                .execute((&key, &form.sequence))?;
            let id = tx.last_insert_rowid();
            for (ordinal, modification) in (0_i64..).zip(&form.modifications) {
                let (kind, position) = match modification.site {
                    ModSite::NTerm => ("n_term", None),
                    ModSite::Residue(position) => ("residue", Some(position)),
                    ModSite::CTerm => ("c_term", None),
                };
                tx.prepare_cached(
                    "INSERT INTO modifications
                        (peptide_id, ordinal, site_kind, position, mass_delta, name, accession)
                     VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
                )?
                .execute(params![
                    id,
                    ordinal,
                    kind,
                    position,
                    modification.mass_delta,
                    modification.name,
                    modification.accession
                ])?;
            }
            report.new_peptides += 1;
            id
        }
    };
    ids.peptides.insert(key, id);
    Ok(id)
}
