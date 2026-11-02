//! Match-between-runs (schema v15): the inputs a transfer needs, and the
//! transfers found. A transfer is a precursor identified in other runs but
//! not in this one, measured where it should elute here; it never replaces
//! anything measured, and imports clear it (as they clear validation).

use std::collections::HashMap;

use rusqlite::{Connection, params};
use tuffold_validate::Level;

use crate::{StoreError, ms1, validation};

/// Transfers are used where their q-value is at most this.
pub const TRANSFER_FDR: f64 = 0.01;

/// One identified precursor of a run.
#[derive(Debug, Clone, PartialEq)]
pub struct MbrPrecursor {
    pub peptide_id: i64,
    pub charge: u8,
    /// The m/z MS1 tracing uses ([`ms1`]'s choice).
    pub mz: f64,
    /// Median retention time of its rank-1 spectra (seconds).
    pub rt: f64,
    /// Median ion mobility (1/K0) of those spectra that have one (timsTOF).
    pub mobility: Option<f64>,
    /// A spectrum of it passes the saved thresholds.
    pub passing: bool,
    /// Apex intensity of its traced MS1 peak in this run, if one was found.
    pub apex_intensity: Option<f64>,
}

/// A run's identified precursors with a retention time.
#[derive(Debug, Clone, PartialEq)]
pub struct MbrRun {
    pub id: i64,
    pub name: String,
    pub precursors: Vec<MbrPrecursor>,
}

/// A transfer into a run, as found.
#[derive(Debug, Clone, PartialEq)]
pub struct Transfer {
    pub peptide_id: i64,
    pub charge: u8,
    pub predicted_rt: f64,
    pub area: Option<f64>,
    pub apex_rt: Option<f64>,
    /// `None`: no peak, beaten by its decoy, or a peak too thin to use
    /// (see the engine's `Support`).
    pub q_value: Option<f64>,
}

/// What matching did in one run.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TransferSummary {
    pub candidates: usize,
    pub known: usize,
    pub trained: bool,
    pub targets_won: usize,
    pub decoys_won: usize,
    pub accepted: usize,
}

/// A spectrum: its key, time, the peptide passing on it (if any), and its
/// precursor's ion mobility.
type SpectrumRow = (String, Option<f64>, Option<i64>, Option<f64>);

/// A precursor's identifications in one run.
struct Seen {
    mz: f64,
    rts: Vec<f64>,
    passing: bool,
    mobilities: Vec<f64>,
}

pub(crate) fn inputs(conn: &mut Connection) -> Result<Vec<MbrRun>, StoreError> {
    let thresholds = validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default();
    let measure = validation::measure_column(thresholds.mode);
    let runs: Vec<(i64, String)> = conn
        .prepare("SELECT id, name FROM ms_runs ORDER BY name")?
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .collect::<Result<_, _>>()?;
    let tx = conn.transaction()?;
    let mut out = Vec::new();
    for (id, name) in runs {
        let by_key = ms1::precursors(&tx, id)?.by_key;
        let apexes: HashMap<(i64, u8), f64> = tx
            .prepare(
                "SELECT peptide_id, charge, apex_intensity FROM precursor_features
                 WHERE ms_run_id = ?1 AND area IS NOT NULL AND apex_intensity > 0",
            )?
            .query_map(params![id], |row| {
                Ok(((row.get(0)?, row.get(1)?), row.get(2)?))
            })?
            .collect::<Result<_, _>>()?;
        // Each spectrum's time and the peptide passing on it, if any.
        let spectra: Vec<SpectrumRow> = tx
            .prepare(&format!(
                "SELECT s.spectrum_key, s.rt_seconds,
                        (SELECT sr.peptide_id FROM spectrum_results sr
                          WHERE sr.spectrum_id = s.id AND sr.is_decoy = 0
                            AND sr.{measure} <= ?2),
                        s.ion_mobility
                 FROM spectra s WHERE s.ms_run_id = ?1"
            ))?
            .query_map(params![id, thresholds.bound(Level::Psm)], |row| {
                Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?))
            })?
            .collect::<Result<_, _>>()?;
        let mut found: HashMap<(i64, u8), Seen> = HashMap::new();
        for (key, rt, passing_peptide, mobility) in spectra {
            let (Some(rt), Some(precursors)) = (rt, by_key.get(&key)) else {
                continue;
            };
            for p in precursors {
                let seen = found.entry((p.peptide_id, p.charge)).or_insert(Seen {
                    mz: p.mz,
                    rts: Vec::new(),
                    passing: false,
                    mobilities: Vec::new(),
                });
                seen.rts.push(rt);
                seen.passing |= passing_peptide == Some(p.peptide_id);
                seen.mobilities.extend(mobility.filter(|m| m.is_finite()));
            }
        }
        let mut precursors: Vec<MbrPrecursor> = found
            .into_iter()
            .map(|((peptide_id, charge), seen)| MbrPrecursor {
                peptide_id,
                charge,
                mz: seen.mz,
                rt: median(seen.rts).unwrap_or(f64::NAN),
                mobility: median(seen.mobilities),
                passing: seen.passing,
                apex_intensity: apexes.get(&(peptide_id, charge)).copied(),
            })
            .collect();
        precursors.sort_by_key(|p| (p.peptide_id, p.charge));
        out.push(MbrRun {
            id,
            name,
            precursors,
        });
    }
    Ok(out)
}

fn median(mut values: Vec<f64>) -> Option<f64> {
    values.sort_by(f64::total_cmp);
    let n = values.len();
    match n {
        0 => None,
        _ if n % 2 == 1 => Some(values[n / 2]),
        _ => Some((values[n / 2 - 1] + values[n / 2]) / 2.0),
    }
}

/// Replace a run's transfers and summary.
pub(crate) fn save(
    conn: &mut Connection,
    run_id: i64,
    transfers: &[Transfer],
    summary: &TransferSummary,
) -> Result<(), StoreError> {
    let tx = conn.transaction()?;
    tx.execute("DELETE FROM transfers WHERE ms_run_id = ?1", [run_id])?;
    tx.execute("DELETE FROM transfer_runs WHERE ms_run_id = ?1", [run_id])?;
    let mut insert = tx.prepare(
        "INSERT INTO transfers
             (ms_run_id, peptide_id, charge, predicted_rt, area, apex_rt, q_value)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
    )?;
    for t in transfers {
        insert.execute(params![
            run_id,
            t.peptide_id,
            t.charge,
            t.predicted_rt,
            t.area,
            t.apex_rt,
            t.q_value
        ])?;
    }
    drop(insert);
    let count = |n: usize| i64::try_from(n).unwrap_or(i64::MAX);
    tx.execute(
        "INSERT INTO transfer_runs
             (ms_run_id, candidates, known, trained, targets_won, decoys_won, accepted)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)",
        params![
            run_id,
            count(summary.candidates),
            count(summary.known),
            summary.trained,
            count(summary.targets_won),
            count(summary.decoys_won),
            count(summary.accepted)
        ],
    )?;
    tx.commit()?;
    Ok(())
}
