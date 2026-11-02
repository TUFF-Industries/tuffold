//! MS1 precursor intensity: which precursors a run's identified spectra
//! point at, the features traced for them (see `tuffold_quant::ms1`), and
//! protein-group intensities at the saved thresholds.
//!
//! A precursor's intensity in a biosample (a peptide at one charge) is the
//! sum of its peak areas over every run of the biosample. Charges stay
//! apart: which charge states are identified varies from run to run, so
//! summing them would make a peptide look more abundant wherever more of
//! them happened to be identified. A group lists the precursors of its unique
//! and razor peptides; protein intensities are built from them in
//! `tuffold_quant::intensity`. A precursor counts in a run only where that
//! run has a passing target spectrum assigned its peptide at its charge. Without match-between-runs,
//! a peptide not identified in a run has no value there, the same as in
//! spectral counting.

use std::collections::HashMap;

use rusqlite::{Connection, Transaction, params};
use tuffold_core::chemistry::{mz, peptide_mass};
use tuffold_quant::ms1::Feature;
use tuffold_validate::Level;

use crate::browse::{Biosample, biosamples};
use crate::{StoreError, validation};

/// A precursor an identified spectrum points at: a rank-1 PSM's peptide at
/// its charge, with the m/z of its monoisotopic peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Precursor {
    pub peptide_id: i64,
    pub charge: u8,
    pub mz: f64,
}

/// When Tuffold's m/z for a PSM's peptide and the engine's own calculated
/// m/z disagree by more than this (ppm), the engine's is traced: the import
/// likely lacks a modification the engine applied (an unreported fixed one).
const ENGINE_MZ_PPM: f64 = 20.0;

/// The precursors a run's identified spectra name.
pub(crate) struct RunPrecursors {
    /// Spectrum key -> the distinct precursors its rank-1 PSMs name.
    pub by_key: HashMap<String, Vec<Precursor>>,
    /// Rank-1 PSMs that can't be traced: no mass for the peptide (a residue
    /// such as `X`) and no calculated m/z from the engine, or no charge.
    pub unplaceable: usize,
}

pub(crate) fn precursors(tx: &Transaction<'_>, run_id: i64) -> Result<RunPrecursors, StoreError> {
    let mut statement = tx.prepare(
        "SELECT s.spectrum_key, p.peptide_id, p.charge, pe.sequence,
                (SELECT COALESCE(SUM(m.mass_delta), 0) FROM modifications m
                  WHERE m.peptide_id = pe.id),
                p.calculated_mz
         FROM psms p
         JOIN spectra s ON s.id = p.spectrum_id
         JOIN peptides pe ON pe.id = p.peptide_id
         WHERE s.ms_run_id = ?1 AND p.rank = 1",
    )?;
    let mut found = RunPrecursors {
        by_key: HashMap::new(),
        unplaceable: 0,
    };
    let mut rows = statement.query([run_id])?;
    while let Some(row) = rows.next()? {
        let (key, peptide_id, charge, sequence, delta, engine_mz): (
            String,
            i64,
            i64,
            String,
            f64,
            Option<f64>,
        ) = (
            row.get(0)?,
            row.get(1)?,
            row.get(2)?,
            row.get(3)?,
            row.get(4)?,
            row.get(5)?,
        );
        let Some(charge) = u8::try_from(charge).ok().filter(|&c| c > 0) else {
            found.unplaceable += 1;
            continue;
        };
        let ours = peptide_mass(&sequence).map(|mass| mz(mass + delta, charge));
        let engine_mz = engine_mz.filter(|v| v.is_finite() && *v > 0.0);
        let traced = match (ours, engine_mz) {
            (Some(ours), Some(engine)) if (ours - engine).abs() > engine * ENGINE_MZ_PPM * 1e-6 => {
                engine
            }
            (Some(ours), _) => ours,
            (None, Some(engine)) => engine,
            (None, None) => {
                found.unplaceable += 1;
                continue;
            }
        };
        let precursors = found.by_key.entry(key).or_default();
        // Engines that agree on the peptide name the same precursor once.
        if !precursors
            .iter()
            .any(|p| (p.peptide_id, p.charge) == (peptide_id, charge))
        {
            precursors.push(Precursor {
                peptide_id,
                charge,
                mz: traced,
            });
        }
    }
    Ok(found)
}

/// Replace the run's features with a new tracing over `ms1_scans` scans.
/// Returns (precursors traced, precursors with a peak, features replaced
/// from an earlier tracing).
pub(crate) fn save_features(
    tx: &Transaction<'_>,
    run_id: i64,
    ms1_scans: usize,
    features: &HashMap<(i64, u8), Option<Feature>>,
) -> Result<(usize, usize, usize), StoreError> {
    let replaced: i64 = tx.query_row(
        "SELECT COUNT(*) FROM precursor_features WHERE ms_run_id = ?1",
        [run_id],
        |row| row.get(0),
    )?;
    tx.execute(
        "UPDATE ms_runs SET ms1_scans = ?2 WHERE id = ?1",
        params![run_id, i64::try_from(ms1_scans).unwrap_or(i64::MAX)],
    )?;
    tx.execute(
        "DELETE FROM precursor_features WHERE ms_run_id = ?1",
        [run_id],
    )?;
    let mut insert = tx.prepare(
        "INSERT INTO precursor_features
             (ms_run_id, peptide_id, charge, area, apex_rt, apex_intensity, start_rt, end_rt,
              points)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
    )?;
    let mut measured = 0;
    for (&(peptide_id, charge), feature) in features {
        measured += usize::from(feature.is_some());
        insert.execute(params![
            run_id,
            peptide_id,
            charge,
            feature.map(|f| f.area),
            feature.map(|f| f.apex_rt),
            feature.map(|f| f.apex_intensity),
            feature.map(|f| f.start_rt),
            feature.map(|f| f.end_rt),
            feature.map(|f| i64::try_from(f.points).unwrap_or(i64::MAX)),
        ])?;
    }
    Ok((
        features.len(),
        measured,
        usize::try_from(replaced).unwrap_or_default(),
    ))
}

/// How much MS1 data a run has.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ms1Run {
    pub run: String,
    pub biosample: String,
    /// MS1 scans read; `None` when no MS1 data was ever attached.
    pub ms1_scans: Option<usize>,
    /// Precursors traced, and of those, ones with a peak.
    pub traced: usize,
    pub measured: usize,
    /// Precursors identified in the run but never traced (imported after
    /// its spectra file was attached).
    pub untraced: usize,
}

/// One precursor (a peptide at one charge) of a group.
#[derive(Debug, Clone, PartialEq)]
pub struct PeptideIntensity {
    pub peptide_id: i64,
    pub proforma: String,
    pub charge: u8,
    /// Summed peak areas, parallel to [`IntensityList::biosamples`].
    pub intensities: Vec<Option<f64>>,
    /// Whether a biosample's value includes a match-between-runs transfer
    /// (parallel to `intensities`).
    pub transferred: Vec<bool>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IntensityListRow {
    pub group_id: i64,
    pub accession: String,
    pub description: Option<String>,
    /// Precursors of unique and razor peptides with a peak in some biosample.
    pub peptides: Vec<PeptideIntensity>,
}

/// Protein groups passing the saved thresholds with their peptides' MS1
/// intensities (in accession order).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct IntensityList {
    pub biosamples: Vec<Biosample>,
    pub rows: Vec<IntensityListRow>,
    pub runs: Vec<Ms1Run>,
    /// Run measurements that came from match-between-runs transfers
    /// (summed into the biosamples' values).
    pub transferred: usize,
}

impl IntensityList {
    /// Whether any run has MS1 data at all.
    pub fn has_ms1(&self) -> bool {
        self.runs
            .iter()
            .any(|run| run.ms1_scans.is_some_and(|n| n > 0))
    }
}

pub(crate) fn intensity_list(conn: &Connection) -> Result<IntensityList, StoreError> {
    let thresholds = validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default();
    let measure = validation::measure_column(thresholds.mode);
    let biosamples = biosamples(conn)?;
    let column: HashMap<i64, usize> = biosamples
        .iter()
        .enumerate()
        .map(|(i, b)| (b.id, i))
        .collect();
    let mut rows: Vec<IntensityListRow> = conn
        .prepare(&format!(
            "SELECT g.id, p.accession, p.description
             FROM protein_groups g JOIN proteins p ON p.id = g.leader_protein_id
             WHERE g.is_decoy = 0 AND g.{measure} <= ?1 AND g.peptides >= ?2"
        ))?
        .query_map(
            params![thresholds.bound(Level::Protein), thresholds.min_peptides],
            |row| {
                Ok(IntensityListRow {
                    group_id: row.get(0)?,
                    accession: row.get(1)?,
                    description: row.get(2)?,
                    peptides: Vec::new(),
                })
            },
        )?
        .collect::<Result<_, _>>()?;
    let index: HashMap<i64, usize> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (row.group_id, i))
        .collect();
    // A precursor's measurement in a run: its own peak where a passing
    // spectrum there names it at its charge, or else an accepted transfer
    // (match-between-runs) of a precursor passing in some run.
    let passing = |run: &str| {
        format!(
            "EXISTS (
               SELECT 1 FROM spectrum_results sr
               JOIN spectra s ON s.id = sr.spectrum_id
               JOIN psms p ON p.spectrum_id = sr.spectrum_id
                   AND p.peptide_id = sr.peptide_id AND p.rank = 1
               WHERE sr.peptide_id = m.peptide_id {run}
                 AND p.charge = m.charge AND sr.is_decoy = 0 AND sr.{measure} <= ?1)"
        )
    };
    let own = passing("AND s.ms_run_id = m.ms_run_id");
    let anywhere = passing("");
    let mut statement = conn.prepare(&format!(
        "WITH measured AS (
             SELECT m.ms_run_id, m.peptide_id, m.charge, m.area, 0 AS transferred
             FROM precursor_features m WHERE m.area IS NOT NULL AND {own}
             UNION ALL
             SELECT m.ms_run_id, m.peptide_id, m.charge, m.area, 1
             FROM transfers m
             WHERE m.area IS NOT NULL AND m.q_value <= ?2 AND {anywhere}
               AND NOT ({own}))
         SELECT gp.group_id, m.peptide_id, pe.proforma, m.charge, r.biosample_id,
                SUM(m.area), SUM(m.transferred)
         FROM protein_group_peptides gp
         JOIN measured m ON m.peptide_id = gp.peptide_id
         JOIN peptides pe ON pe.id = m.peptide_id
         JOIN ms_runs r ON r.id = m.ms_run_id
         WHERE gp.role IN ('unique', 'razor')
         GROUP BY gp.group_id, m.peptide_id, m.charge, r.biosample_id
         ORDER BY gp.group_id, pe.proforma, m.charge"
    ))?;
    let mut found = statement.query(params![
        thresholds.bound(Level::Psm),
        crate::mbr::TRANSFER_FDR
    ])?;
    let mut transferred = 0;
    while let Some(row) = found.next()? {
        let (group, peptide_id, charge, biosample): (i64, i64, u8, i64) =
            (row.get(0)?, row.get(1)?, row.get(3)?, row.get(4)?);
        let (Some(&i), Some(&s)) = (index.get(&group), column.get(&biosample)) else {
            continue;
        };
        let peptides = &mut rows[i].peptides;
        if peptides
            .last()
            .is_none_or(|p| (p.peptide_id, p.charge) != (peptide_id, charge))
        {
            peptides.push(PeptideIntensity {
                peptide_id,
                proforma: row.get(2)?,
                charge,
                intensities: vec![None; biosamples.len()],
                transferred: vec![false; biosamples.len()],
            });
        }
        let from_transfers = usize::try_from(row.get::<_, i64>(6)?).unwrap_or_default();
        if let Some(peptide) = peptides.last_mut() {
            peptide.intensities[s] = row.get(5)?;
            peptide.transferred[s] = from_transfers > 0;
        }
        transferred += from_transfers;
    }
    rows.sort_by(|a, b| a.accession.cmp(&b.accession));
    Ok(IntensityList {
        biosamples,
        rows,
        runs: ms1_runs(conn)?,
        transferred,
    })
}

fn ms1_runs(conn: &Connection) -> Result<Vec<Ms1Run>, StoreError> {
    let count = |value: i64| usize::try_from(value).unwrap_or_default();
    conn.prepare(
        "SELECT r.name, b.name, r.ms1_scans,
                (SELECT COUNT(*) FROM precursor_features f WHERE f.ms_run_id = r.id),
                (SELECT COUNT(*) FROM precursor_features f
                  WHERE f.ms_run_id = r.id AND f.area IS NOT NULL),
                (SELECT COUNT(*) FROM (
                    SELECT DISTINCT p.peptide_id, p.charge FROM psms p
                    JOIN spectra s ON s.id = p.spectrum_id
                    WHERE s.ms_run_id = r.id AND p.rank = 1
                      AND NOT EXISTS (
                          SELECT 1 FROM precursor_features f
                          WHERE f.ms_run_id = r.id AND f.peptide_id = p.peptide_id
                            AND f.charge = p.charge)))
         FROM ms_runs r JOIN biosamples b ON b.id = r.biosample_id
         ORDER BY r.name",
    )?
    .query_map([], |row| {
        Ok(Ms1Run {
            run: row.get(0)?,
            biosample: row.get(1)?,
            ms1_scans: row.get::<_, Option<i64>>(2)?.map(count),
            traced: count(row.get(3)?),
            measured: count(row.get(4)?),
            untraced: count(row.get(5)?),
        })
    })?
    .collect::<Result<_, _>>()
    .map_err(Into::into)
}
