//! Attaching measured spectra (peaks) to identified spectra.
//!
//! Only spectra that some PSM already references are stored, so a project
//! carries the evidence it needs without duplicating whole raw files. A
//! spectrum attaches when its run and key equal an identified spectrum's;
//! retention time and precursor values fill in only where still unknown, and
//! peaks once stored are never replaced, except partial ones: the peaks a
//! project file kept (a Scaffold project's at most 300) give way to a
//! spectra file's. Keys are matched in memory before
//! any peak data is compressed, so unidentified spectra cost almost nothing.
//!
//! An MS3 scan naming an identified spectrum as its parent (SPS-MS3, whose
//! TMT reporters are in the MS3) is stored on that spectrum
//! (`spectra.ms3_peaks`), once, whatever peaks the spectrum has.
//!
//! The attacher also knows which precursors each identified spectrum names
//! ([`Attacher::precursors`]), so the caller can trace them through the
//! file's MS1 scans and record the features in the same transaction
//! ([`Attacher::record_ms1`]).

use std::collections::HashMap;

use rusqlite::{OptionalExtension, Transaction, params};
use tuffold_core::PeakSpectrum;

use crate::StoreError;
use crate::ms1::{self, Precursor};
use crate::peaks::{self, Peaks};
use tuffold_core::SpectrumKey;
use tuffold_quant::ms1::Feature;

/// What one spectra file contributed.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AttachReport {
    pub run: String,
    /// MS/MS spectra read from the file.
    pub read: usize,
    /// Of those, spectra matching an identified spectrum (peaks stored).
    pub attached: usize,
    /// Spectra skipped because their identified spectrum already had peaks
    /// (from an earlier file, or an earlier spectrum with the same key).
    pub already_attached: usize,
    /// Identified spectra of this run that still have no peaks.
    pub identified_without_peaks: usize,
    /// MS3 scans read, and of those, ones stored on the identified MS/MS
    /// spectrum they name as parent (SPS-MS3). An MS3 scan whose parent
    /// isn't identified, or already has one, is skipped.
    pub ms3_scans: usize,
    pub ms3_attached: usize,
    /// Identified spectra whose precursor purity was measured (MS1 scans
    /// around them and an isolation window in the file).
    pub purity_measured: usize,
    /// MS1 scans read (0 when the file has none: MGF, or MS/MS-only mzML).
    pub ms1_scans: usize,
    /// Precursors traced through those scans, and of those, ones with a
    /// peak. Both 0 without MS1 scans.
    pub precursors_traced: usize,
    pub precursors_measured: usize,
    /// Features of an earlier tracing of the run that this one replaced.
    pub precursors_replaced: usize,
    /// Rank-1 PSMs of the run whose precursor can't be placed in m/z (see
    /// [`Attacher::precursors`]).
    pub precursors_unplaceable: usize,
    /// MS1 scans without a retention time, which can't be placed in time
    /// and aren't traced (counted by whoever read the scans).
    pub ms1_scans_without_rt: usize,
    /// MS1 scans that can't be trusted (a timsTOF frame that won't decode,
    /// an mzML scan holding impossible peaks), left out of tracing and
    /// purity (counted by whoever read the scans).
    pub ms1_scans_damaged: usize,
}

/// An identified spectrum of the run: its row id and whether it has peaks
/// and MS3 peaks.
#[derive(Clone, Copy)]
struct Target {
    id: i64,
    has_peaks: bool,
    has_ms3: bool,
}

/// An open attach operation on one run; commits on [`Attacher::finish`].
pub struct Attacher<'p> {
    tx: Transaction<'p>,
    run_id: i64,
    /// These peaks are partial (from a project file), filling only spectra
    /// without any.
    partial: bool,
    targets: HashMap<String, Target>,
    precursors: HashMap<String, Vec<Precursor>>,
    report: AttachReport,
}

impl<'p> Attacher<'p> {
    pub(crate) fn begin(tx: Transaction<'p>, run: &str, partial: bool) -> Result<Self, StoreError> {
        let run_id: Option<i64> = tx
            .query_row("SELECT id FROM ms_runs WHERE name = ?1", [run], |row| {
                row.get(0)
            })
            .optional()?;
        let Some(run_id) = run_id else {
            let mut statement = tx.prepare("SELECT name FROM ms_runs ORDER BY name")?;
            let known = statement
                .query_map([], |row| row.get::<_, String>(0))?
                .collect::<Result<Vec<_>, _>>()?;
            drop(statement);
            return Err(StoreError::UnknownRun {
                run: run.to_string(),
                known,
            });
        };
        let targets = {
            // Partial peaks count as present only to more partial peaks.
            let mut statement = tx.prepare(
                "SELECT spectrum_key, id, peaks IS NOT NULL AND (peaks_partial = 0 OR ?2),
                        ms3_peaks IS NOT NULL
                 FROM spectra WHERE ms_run_id = ?1",
            )?;
            statement
                .query_map(params![run_id, partial], |row| {
                    Ok((
                        row.get::<_, String>(0)?,
                        Target {
                            id: row.get(1)?,
                            has_peaks: row.get(2)?,
                            has_ms3: row.get(3)?,
                        },
                    ))
                })?
                .collect::<Result<HashMap<_, _>, _>>()?
        };
        let found = ms1::precursors(&tx, run_id)?;
        Ok(Self {
            tx,
            run_id,
            partial,
            targets,
            precursors: found.by_key,
            report: AttachReport {
                run: run.to_string(),
                precursors_unplaceable: found.unplaceable,
                ..AttachReport::default()
            },
        })
    }

    pub fn add(&mut self, spectrum: &PeakSpectrum) -> Result<(), StoreError> {
        if spectrum.ms_level >= 3 {
            return self.add_ms3(spectrum);
        }
        self.report.read += 1;
        let Some(target) = self.targets.get_mut(&spectrum.key.canonical()) else {
            return Ok(()); // not identified: nothing to attach to
        };
        if target.has_peaks {
            self.report.already_attached += 1;
            // Peaks attached before mobility was recorded get it now.
            if let Some(mobility) = spectrum.precursor_mobility {
                self.tx
                    .prepare_cached(
                        "UPDATE spectra SET ion_mobility = COALESCE(ion_mobility, ?2)
                         WHERE id = ?1",
                    )?
                    .execute(params![target.id, mobility])?;
            }
            return Ok(());
        }
        target.has_peaks = true;
        self.tx
            .prepare_cached(
                "UPDATE spectra SET
                    peaks = ?2,
                    peaks_partial = ?6,
                    rt_seconds = COALESCE(rt_seconds, ?3),
                    precursor_mz = COALESCE(precursor_mz, ?4),
                    precursor_charge = COALESCE(precursor_charge, ?5),
                    ion_mobility = COALESCE(ion_mobility, ?7)
                 WHERE id = ?1 AND (peaks IS NULL OR peaks_partial = 1)",
            )?
            .execute(params![
                target.id,
                peaks::encode(&spectrum.mz, &spectrum.intensity),
                spectrum.rt_seconds,
                spectrum.precursor_mz,
                spectrum.precursor_charge,
                self.partial,
                spectrum.precursor_mobility,
            ])?;
        self.report.attached += 1;
        Ok(())
    }

    /// Store an MS3 scan's peaks on the identified spectrum it names as
    /// parent.
    fn add_ms3(&mut self, spectrum: &PeakSpectrum) -> Result<(), StoreError> {
        self.report.ms3_scans += 1;
        let Some(target) = spectrum
            .parent
            .as_deref()
            .map(|parent| SpectrumKey::from_native_id(parent).canonical())
            .and_then(|key| self.targets.get_mut(&key))
        else {
            return Ok(());
        };
        if target.has_ms3 {
            return Ok(());
        }
        target.has_ms3 = true;
        self.tx
            .prepare_cached("UPDATE spectra SET ms3_peaks = ?2 WHERE id = ?1")?
            .execute(params![
                target.id,
                peaks::encode(&spectrum.mz, &spectrum.intensity)
            ])?;
        self.report.ms3_attached += 1;
        Ok(())
    }

    /// The precursors that rank-1 PSMs of the spectrum `key` name, at the m/z
    /// of their monoisotopic peak: Tuffold's own for the peptide, or the
    /// engine's calculated m/z where the two disagree by more than 20 ppm
    /// (the import likely lacks a modification the engine applied).
    pub fn precursors(&self, key: &SpectrumKey) -> &[Precursor] {
        self.precursors
            .get(&key.canonical())
            .map_or(&[], Vec::as_slice)
    }

    /// Record an identified spectrum's precursor purity (see
    /// `tuffold_quant::purity`); `None` leaves it unknown.
    pub fn record_purity(&mut self, key: &str, purity: Option<f64>) -> Result<(), StoreError> {
        let (Some(target), Some(purity)) = (self.targets.get(key), purity) else {
            return Ok(());
        };
        self.tx
            .prepare_cached("UPDATE spectra SET precursor_purity = ?2 WHERE id = ?1")?
            .execute(params![target.id, purity])?;
        self.report.purity_measured += 1;
        Ok(())
    }

    /// Replace the run's precursor features with ones traced through
    /// `ms1_scans` MS1 scans. Without MS1 scans nothing changes: earlier
    /// features stay.
    pub fn record_ms1(
        &mut self,
        ms1_scans: usize,
        features: &HashMap<(i64, u8), Option<Feature>>,
    ) -> Result<(), StoreError> {
        self.report.ms1_scans = ms1_scans;
        if ms1_scans == 0 {
            return Ok(());
        }
        let (traced, measured, replaced) =
            ms1::save_features(&self.tx, self.run_id, ms1_scans, features)?;
        self.report.precursors_traced = traced;
        self.report.precursors_measured = measured;
        self.report.precursors_replaced = replaced;
        Ok(())
    }

    pub fn finish(mut self) -> Result<AttachReport, StoreError> {
        let missing: i64 = self.tx.query_row(
            "SELECT COUNT(*) FROM spectra s
             WHERE s.ms_run_id = ?1 AND s.peaks IS NULL
               AND EXISTS (SELECT 1 FROM psms p WHERE p.spectrum_id = s.id)",
            [self.run_id],
            |row| row.get(0),
        )?;
        self.report.identified_without_peaks = usize::try_from(missing).unwrap_or_default();
        self.tx.commit()?;
        Ok(self.report)
    }
}

/// Peaks stored for one spectrum, if any.
pub(crate) fn stored_peaks(
    conn: &rusqlite::Connection,
    run: &str,
    key: &str,
) -> Result<Option<Peaks>, StoreError> {
    let blob: Option<Option<Vec<u8>>> = conn
        .query_row(
            "SELECT s.peaks FROM spectra s JOIN ms_runs r ON r.id = s.ms_run_id
             WHERE r.name = ?1 AND s.spectrum_key = ?2",
            [run, key],
            |row| row.get(0),
        )
        .optional()?;
    Ok(blob.flatten().and_then(|blob| peaks::decode(&blob)))
}
