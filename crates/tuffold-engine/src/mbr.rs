//! Match-between-runs: a precursor identified in other runs but not in this
//! one is measured where it should elute here, and the found peaks are
//! judged by target-decoy competition (`tuffold_validate::transfer`).
//!
//! It runs after validation, on the runs' spectra files (read again for
//! their MS1 scans): donors are passing precursors, and each other run's
//! retention times are mapped onto this run's through the precursors both
//! identified (`tuffold_quant::align`). A precursor's predicted time is the
//! median of its donors' mapped times. Every looked-for peak gets a decoy
//! at a time borrowed from another identified precursor of this run (at
//! least [`DECOY_DISTANCE`] alignment spreads away), and the run's own
//! passing precursors, predicted the same way from the other runs, are the
//! known transfers the model trains on. Accepted transfers (q ≤ 1%) fill a
//! run's value only where the run identified nothing of that precursor, and
//! only if the peak is as well supported as nearly all of the run's own
//! identifications ([`SUPPORT_QUANTILE`]). Normalization never uses them
//! (`tuffold_quant::intensity`).
//!
//! On timsTOF runs a transfer is also measured at its ion mobility: the
//! median of its donors' mobilities, each donor's shifted by the median
//! difference between this run's and the donor's mobilities over the
//! precursors both passed (calibration drifts a little between runs). A
//! decoy keeps its target's m/z and mobility and only borrows a time.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use tuffold_core::centroid::centroid;
use tuffold_io::{SpectraError, read_scans, spectra_run_name};
use tuffold_quant::align::{Alignment, align};
use tuffold_quant::ms1::{Feature, TraceSettings, Tracer};
use tuffold_store::{MbrRun, Project, StoreError, TRANSFER_FDR, Transfer, TransferSummary};
use tuffold_validate::transfer::{MIN_KNOWN, TransferKind, TransferRow, judge};

/// A decoy's borrowed time is at least this many alignment spreads (and
/// [`MIN_SPREAD`] seconds' worth) from the predicted one.
pub const DECOY_DISTANCE: f64 = 4.0;
/// Alignment spreads below this (seconds) are scored as this: no alignment
/// is that good, and a near-zero spread would make every error huge.
pub const MIN_SPREAD: f64 = 10.0;
/// A transferred peak must be as well supported (apex intensity and scans)
/// as all but this share of the run's own identified precursors.
pub const SUPPORT_QUANTILE: f64 = 0.05;
const PROTON: f64 = 1.007_276_467;
/// Expected M+1 / M per dalton of peptide mass, from averagine (Senko et
/// al. 1995: C 4.9384, H 7.7583, N 1.3577, O 1.4773, S 0.0417 per
/// 111.1254 Da) and each element's heavy-isotope ratio (¹³C/¹²C 0.010816,
/// ²H/¹H 0.000115, ¹⁵N/¹⁴N 0.003673, ¹⁷O/¹⁶O 0.000381, ³³S/³²S 0.0079):
/// (0.053414 + 0.000892 + 0.004987 + 0.000563 + 0.000329) / 111.1254.
const M1_PER_DALTON: f64 = 0.060_185 / 111.125_4;

#[derive(Debug, thiserror::Error)]
pub enum MbrError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error("this project has no current validation results; validate it first")]
    NotValidated,
    #[error("{file}: no run named '{run}' in this project (runs: {})", known.join(", "))]
    UnknownRun {
        file: String,
        run: String,
        known: Vec<String>,
    },
    #[error("{0}")]
    Read(String),
}

/// What matching did in one run.
#[derive(Debug, Clone, PartialEq)]
pub struct MbrRunReport {
    pub run: String,
    pub summary: TransferSummary,
    /// Other runs aligned to this one, of those with identifications.
    pub aligned: usize,
    pub donors: usize,
    /// MS1 scans of the run's spectra file that were left out as damaged
    /// (a timsTOF frame that won't decode, an mzML scan with impossible
    /// peaks), as when the file was attached.
    pub ms1_scans_damaged: usize,
}

/// A precursor passing in other runs: its m/z, and per donor its time
/// mapped onto this run (with the alignment's spread), its mobility shifted
/// onto this run, and its log10 apex intensity put on this run's scale.
#[derive(Debug, Default)]
struct Donors {
    mz: f64,
    mapped: Vec<(f64, f64)>,
    mobilities: Vec<f64>,
    intensities: Vec<f64>,
}

type Predicted = HashMap<(i64, u8), Donors>;

/// One looked-for peak.
#[derive(Debug, Clone, Copy)]
struct Looked {
    pair: usize,
    kind: TransferKind,
    peptide_id: i64,
    charge: u8,
    mz: f64,
    predicted: f64,
    /// Where it is looked for (the predicted time, or a decoy's borrowed one).
    at: f64,
    spread: f64,
    /// Ion mobility (1/K0) it is looked for at, when the donors had one.
    mobility: Option<f64>,
    /// The donors' median log10 apex intensity on this run's scale.
    intensity: Option<f64>,
}

/// Match between the runs whose spectra `files` are given.
pub fn match_between_runs(
    project: &mut Project,
    files: &[PathBuf],
) -> Result<Vec<MbrRunReport>, MbrError> {
    let validated = project
        .validation_summary()?
        .is_some_and(|summary| summary.validated_at.is_some());
    if !validated {
        return Err(MbrError::NotValidated);
    }
    let runs = project.mbr_inputs()?;
    let mut reports = Vec::new();
    for file in files {
        let name = spectra_run_name(file);
        let Some(acceptor) = runs.iter().find(|r| r.name == name) else {
            return Err(MbrError::UnknownRun {
                file: file.display().to_string(),
                run: name,
                known: runs.iter().map(|r| r.name.clone()).collect(),
            });
        };
        let scans = read_scans(file).map_err(|e| MbrError::Read(e.to_string()))?;
        reports.push(match_run(project, &runs, acceptor, scans)?);
    }
    Ok(reports)
}

/// Match into one run, given its scans in file order.
pub(crate) fn match_run(
    project: &mut Project,
    runs: &[MbrRun],
    acceptor: &MbrRun,
    scans: impl Iterator<Item = Result<tuffold_core::PeakSpectrum, SpectraError>>,
) -> Result<MbrRunReport, MbrError> {
    let (looked, aligned, donors) = plan(acceptor, runs);
    let (features, ms1_scans_damaged) = trace(scans, &looked)?;
    let (transfers, summary) = judged(&looked, &features);
    project.save_transfers(acceptor.id, &transfers, &summary)?;
    Ok(MbrRunReport {
        run: acceptor.name.clone(),
        summary,
        aligned,
        donors,
        ms1_scans_damaged,
    })
}

/// What to look for in `acceptor`, and how many other runs were aligned
/// (of those with identifications).
fn plan(acceptor: &MbrRun, runs: &[MbrRun]) -> (Vec<Looked>, usize, usize) {
    let own: HashMap<(i64, u8), f64> = acceptor
        .precursors
        .iter()
        .map(|p| ((p.peptide_id, p.charge), p.rt))
        .collect();
    let passing_here: HashMap<(i64, u8), f64> = acceptor
        .precursors
        .iter()
        .filter(|p| p.passing)
        .map(|p| ((p.peptide_id, p.charge), p.rt))
        .collect();
    let mobility_here: HashMap<(i64, u8), f64> = acceptor
        .precursors
        .iter()
        .filter(|p| p.passing)
        .filter_map(|p| p.mobility.map(|m| ((p.peptide_id, p.charge), m)))
        .collect();
    let intensity_here: HashMap<(i64, u8), f64> = acceptor
        .precursors
        .iter()
        .filter(|p| p.passing)
        .filter_map(|p| Some(((p.peptide_id, p.charge), log_intensity(p.apex_intensity?)?)))
        .collect();
    // Each passing precursor of the other runs: its mapped times and spreads.
    let mut predicted: Predicted = HashMap::new();
    let (mut aligned, mut donors) = (0, 0);
    for donor in runs.iter().filter(|r| r.id != acceptor.id) {
        if donor.precursors.is_empty() {
            continue;
        }
        donors += 1;
        let pairs: Vec<(f64, f64)> = donor
            .precursors
            .iter()
            .filter(|p| p.passing)
            .filter_map(|p| {
                passing_here
                    .get(&(p.peptide_id, p.charge))
                    .map(|&r| (p.rt, r))
            })
            .collect();
        let Some(alignment) = align(&pairs) else {
            continue;
        };
        aligned += 1;
        let shifts: Vec<f64> = donor
            .precursors
            .iter()
            .filter(|p| p.passing)
            .filter_map(|p| Some(mobility_here.get(&(p.peptide_id, p.charge))? - p.mobility?))
            .collect();
        let shift = if shifts.is_empty() {
            0.0
        } else {
            median(&shifts)
        };
        // How much brighter this run reads than the donor: the median log10
        // apex ratio over precursors both passed with a peak in both (at
        // least MIN_KNOWN of them, or the donor's intensities aren't used).
        let ratios: Vec<f64> = donor
            .precursors
            .iter()
            .filter(|p| p.passing)
            .filter_map(|p| {
                Some(
                    intensity_here.get(&(p.peptide_id, p.charge))?
                        - log_intensity(p.apex_intensity?)?,
                )
            })
            .collect();
        let scale = (ratios.len() >= MIN_KNOWN).then(|| median(&ratios));
        add_donor(&mut predicted, donor, &alignment, shift, scale);
    }
    // The run's identified times, to borrow decoy times from.
    let mut pool: Vec<f64> = own.values().copied().collect();
    pool.sort_by(f64::total_cmp);
    let mut keys: Vec<&(i64, u8)> = predicted.keys().collect();
    keys.sort();
    let mut looked = Vec::new();
    for (pair, key) in keys.into_iter().enumerate() {
        let Donors {
            mz,
            mapped,
            mobilities,
            intensities,
        } = &predicted[key];
        let kind = if passing_here.contains_key(key) {
            TransferKind::Known
        } else if own.contains_key(key) {
            continue; // identified here, but not passing: neither
        } else {
            TransferKind::Target
        };
        let times: Vec<f64> = mapped.iter().map(|m| m.0).collect();
        let spreads: Vec<f64> = mapped.iter().map(|m| m.1).collect();
        let at = median(&times);
        let spread = median(&spreads).max(MIN_SPREAD);
        let Some(decoy_at) = borrowed(&pool, at, spread, key) else {
            continue; // no time far enough: neither looked for
        };
        let row = Looked {
            pair,
            kind,
            peptide_id: key.0,
            charge: key.1,
            mz: *mz,
            predicted: at,
            at,
            spread,
            mobility: (!mobilities.is_empty()).then(|| median(mobilities)),
            intensity: (!intensities.is_empty()).then(|| median(intensities)),
        };
        let decoy_kind = match kind {
            TransferKind::Known => TransferKind::KnownDecoy,
            _ => TransferKind::Decoy,
        };
        looked.push(row);
        // The decoy is a transfer predicted at the borrowed time: judged by
        // its error from there, as a wrong peak at a plausible time would be.
        // (Judged from the target's time, every decoy lost on retention time
        // alone, and on iPRG 2016 2.7% of transfers at "1%" were peptides of
        // the other mixture's PrEST pool.)
        looked.push(Looked {
            kind: decoy_kind,
            predicted: decoy_at,
            at: decoy_at,
            ..row
        });
    }
    (looked, aligned, donors)
}

/// A donor's passing precursors, their times mapped onto this run, their
/// mobilities shifted by `mobility_shift`, and their log10 apex intensities
/// by `intensity_scale` (none without one).
fn add_donor(
    predicted: &mut Predicted,
    donor: &MbrRun,
    alignment: &Alignment,
    mobility_shift: f64,
    intensity_scale: Option<f64>,
) {
    for p in donor.precursors.iter().filter(|p| p.passing) {
        let entry = predicted
            .entry((p.peptide_id, p.charge))
            .or_insert_with(|| Donors {
                mz: p.mz,
                ..Donors::default()
            });
        entry.mapped.push((alignment.map(p.rt), alignment.spread));
        entry
            .mobilities
            .extend(p.mobility.map(|m| m + mobility_shift));
        if let (Some(apex), Some(scale)) =
            (p.apex_intensity.and_then(log_intensity), intensity_scale)
        {
            entry.intensities.push(apex + scale);
        }
    }
}

/// A decoy time: an identified time of the run at least
/// [`DECOY_DISTANCE`] spreads from `at`, picked by a hash of the precursor
/// so the choice is fixed.
fn borrowed(pool: &[f64], at: f64, spread: f64, key: &(i64, u8)) -> Option<f64> {
    if pool.is_empty() {
        return None;
    }
    let far = DECOY_DISTANCE * spread;
    let start = usize::try_from(mix(key.0, key.1) % pool.len() as u64).unwrap_or(0);
    (0..pool.len())
        .map(|step| pool[(start + step) % pool.len()])
        .find(|&t| (t - at).abs() > far)
}

/// splitmix64 of the precursor: a fixed, well-spread choice.
fn mix(peptide_id: i64, charge: u8) -> u64 {
    let mut z = (peptide_id as u64)
        .wrapping_mul(31)
        .wrapping_add(u64::from(charge))
        .wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Trace every looked-for peak through the run's MS1 scans, leaving out (and
/// counting) damaged ones.
fn trace(
    scans: impl Iterator<Item = Result<tuffold_core::PeakSpectrum, SpectraError>>,
    looked: &[Looked],
) -> Result<(HashMap<usize, Option<Feature>>, usize), MbrError> {
    let mut order: Vec<usize> = (0..looked.len()).collect();
    order.sort_by(|&a, &b| looked[a].at.total_cmp(&looked[b].at));
    let mut next = 0;
    let mut tracer = Tracer::new(TraceSettings::default());
    let mut damaged = 0;
    for scan in scans {
        let scan = match scan {
            Ok(scan) => scan,
            Err(SpectraError::DamagedMs1 { .. }) => {
                damaged += 1;
                continue;
            }
            Err(error) => return Err(MbrError::Read(error.to_string())),
        };
        let Some(rt) = scan.rt_seconds.filter(|_| scan.ms_level == 1) else {
            continue;
        };
        while let Some(&i) = order.get(next).filter(|&&i| looked[i].at <= rt) {
            let l = &looked[i];
            tracer.add_identification(i, l.mz, l.charge, l.at, l.mobility);
            next += 1;
        }
        if scan.profile {
            let (mz, intensity) = centroid(&scan.mz, &scan.intensity);
            tracer.add_ms1(rt, &mz, &intensity, &[]);
        } else {
            tracer.add_ms1(rt, &scan.mz, &scan.intensity, &scan.mobility);
        }
    }
    Ok((tracer.finish(), damaged))
}

/// How thin a transferred peak may be and still be used.
///
/// Target-decoy competition says a peak is the right *precursor*, not that
/// its area is a good measurement: transfers that win are weaker than the
/// run's identifications (a precursor not picked for MS/MS usually is), and
/// the weakest, near the noise floor, are cut short (few scans above it) and
/// read far too low. So a winner is used only if its peak reaches the apex
/// intensity and has the scans that [`SUPPORT_QUANTILE`] of the run's own
/// identified precursors do (found at the same predicted times). The floors
/// come from the run's data, apply to winners after the 1% judgement (so
/// they cannot change what competes), and are zero when too few known peaks
/// exist to set them.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Support {
    apex_intensity: f64,
    points: usize,
}

impl Support {
    /// The floors from the run's known peaks (identified precursors looked
    /// for where the other runs predict them).
    fn of_known(peaks: &[Feature]) -> Self {
        if peaks.len() < MIN_KNOWN {
            return Self::default();
        }
        let mut apexes: Vec<f64> = peaks.iter().map(|f| f.apex_intensity).collect();
        apexes.sort_by(f64::total_cmp);
        let mut points: Vec<usize> = peaks.iter().map(|f| f.points).collect();
        points.sort_unstable();
        // The value with ⌊n · q⌋ smaller ones.
        let below = (peaks.len() as f64 * SUPPORT_QUANTILE) as usize;
        Self {
            apex_intensity: apexes[below],
            points: points[below],
        }
    }

    fn holds(&self, peak: &Feature) -> bool {
        peak.apex_intensity >= self.apex_intensity && peak.points >= self.points
    }
}

/// Judge the found peaks; the targets as transfers, with what happened.
fn judged(
    looked: &[Looked],
    features: &HashMap<usize, Option<Feature>>,
) -> (Vec<Transfer>, TransferSummary) {
    let found = |i: usize| features.get(&i).copied().flatten();
    // The run's m/z offset: the known peaks' median error.
    let known_ppm: Vec<f64> = (0..looked.len())
        .filter(|&i| looked[i].kind == TransferKind::Known)
        .filter_map(|i| found(i).map(|f| f.apex_ppm))
        .collect();
    let center = if known_ppm.is_empty() {
        0.0
    } else {
        median(&known_ppm)
    };
    let log_apex = |f: &Feature| log_intensity(f.apex_intensity).unwrap_or(0.0);
    // A peak whose donors had no measured peak gets the median intensity
    // error of every found peak that has one (the same for a target and its
    // decoy, which share donors).
    let errors: Vec<f64> = (0..looked.len())
        .filter_map(|i| Some((log_apex(&found(i)?) - looked[i].intensity?).abs()))
        .collect();
    let typical = if errors.is_empty() {
        0.0
    } else {
        median(&errors)
    };
    let rows: Vec<TransferRow> = looked
        .iter()
        .enumerate()
        .map(|(i, l)| TransferRow {
            pair: l.pair,
            kind: l.kind,
            features: found(i).map(|f| {
                let mass = (l.mz - PROTON) * f64::from(l.charge);
                let expected = (mass * M1_PER_DALTON).max(1e-3);
                [
                    (f.apex_rt - l.predicted).abs() / l.spread,
                    (f.apex_ppm - center).abs(),
                    (f.apex_ratio / expected).ln().abs(),
                    log_apex(&f),
                    f.points as f64,
                    l.intensity.map_or(typical, |e| (log_apex(&f) - e).abs()),
                ]
            }),
        })
        .collect();
    let judgement = judge(&rows);
    let support = Support::of_known(
        &(0..looked.len())
            .filter(|&i| looked[i].kind == TransferKind::Known)
            .filter_map(found)
            .collect::<Vec<_>>(),
    );
    let mut transfers = Vec::new();
    let mut summary = TransferSummary {
        known: judgement.known.0,
        trained: judgement.trained,
        targets_won: judgement.won.0,
        decoys_won: judgement.won.1,
        ..TransferSummary::default()
    };
    let mut seen = HashSet::new();
    for (i, l) in looked.iter().enumerate() {
        if l.kind != TransferKind::Target || !seen.insert(l.pair) {
            continue;
        }
        summary.candidates += 1;
        let feature = found(i);
        // A winner whose peak is too thin to measure well is not used.
        let q = judgement.q[i].filter(|_| feature.is_some_and(|f| support.holds(&f)));
        if q.is_some_and(|q| q <= TRANSFER_FDR) {
            summary.accepted += 1;
        }
        transfers.push(Transfer {
            peptide_id: l.peptide_id,
            charge: l.charge,
            predicted_rt: l.predicted,
            area: feature.map(|f| f.area),
            apex_rt: feature.map(|f| f.apex_rt),
            q_value: q,
        });
    }
    (transfers, summary)
}

/// log10 of an apex intensity, floored at 1; none for a value no peak has.
fn log_intensity(apex: f64) -> Option<f64> {
    (apex.is_finite() && apex > 0.0).then(|| apex.max(1.0).log10())
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decoy_times_are_borrowed_far_from_the_prediction() {
        // Identified times 100..=1000 s every 100; a prediction at 500 with
        // spread 20 needs a time more than 4 × 20 = 80 s away: anything but
        // 500 itself (the nearest others are 100 s away).
        let pool: Vec<f64> = (1..=10).map(|i| f64::from(i) * 100.0).collect();
        for key in [(1, 2), (7, 3), (42, 2)] {
            let t = borrowed(&pool, 500.0, 20.0, &key).unwrap();
            assert!((t - 500.0).abs() > 80.0 && pool.contains(&t), "{t}");
            // The same precursor always borrows the same time.
            assert_eq!(borrowed(&pool, 500.0, 20.0, &key), Some(t));
        }
        // Nothing far enough: no decoy (and so no transfer is looked for).
        assert_eq!(borrowed(&[480.0, 520.0], 500.0, 20.0, &(1, 2)), None);
    }

    #[test]
    fn damaged_ms1_scans_are_skipped_and_counted_and_other_errors_stop() {
        let damaged = |scan: &str| SpectraError::DamagedMs1 {
            file: "run.mzML".into(),
            scan: scan.into(),
            problem: "impossible peaks".into(),
        };
        let ms1 = tuffold_core::PeakSpectrum {
            key: tuffold_core::SpectrumKey::Scan(2),
            native_id: "scan=2".into(),
            ms_level: 1,
            rt_seconds: Some(60.0),
            precursor_mz: None,
            precursor_charge: None,
            mz: vec![500.0],
            intensity: vec![10.0],
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        };

        let (_, count) = trace(
            [Err(damaged("scan=1")), Ok(ms1), Err(damaged("scan=3"))].into_iter(),
            &[],
        )
        .unwrap();
        assert_eq!(count, 2);

        let invalid = SpectraError::Invalid {
            file: "run.mzML".into(),
            id: "scan=4".into(),
            problem: "bad".into(),
        };
        assert!(matches!(
            trace([Err(invalid)].into_iter(), &[]),
            Err(MbrError::Read(_))
        ));
    }

    /// The demo's MS1 scans (with chromatography on) hold elution peaks only
    /// for each run's identified precursors. With control-A's identifications of every
    /// fifth precursor dropped (their peaks stay), those are the only true
    /// transfers into it: every other candidate is absent from its scans.
    #[test]
    fn transfers_recover_dropped_identifications_and_not_absent_ones() {
        use tuffold_core::SourceFormat;
        use tuffold_store::ImportSource;
        use tuffold_validate::Settings;

        use crate::demo::{DEMO_SEED, generate_with};
        use crate::import::attach_scans;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mbr.tuff");
        let demo = generate_with(DEMO_SEED, true);
        let fifth = |sequence: &str| {
            mix(
                sequence.len() as i64 * 7919 + i64::from(sequence.as_bytes()[0]),
                0,
            )
            .is_multiple_of(5)
        };
        let mut batch = demo.batch.clone();
        let dropped: HashSet<(String, u8)> = batch
            .psms
            .iter()
            .filter(|p| p.ms_run == "control-A" && p.rank == 1 && fifth(&p.peptide.sequence))
            .map(|p| {
                (
                    p.peptide.sequence.to_string(),
                    u8::try_from(p.charge).unwrap(),
                )
            })
            .collect();
        batch
            .psms
            .retain(|p| !(p.ms_run == "control-A" && fifth(&p.peptide.sequence)));
        let mut project = Project::create(&path).unwrap();
        let source = ImportSource {
            path: "MBR test (synthetic)".into(),
            format: SourceFormat::Synthetic,
        };
        project.import(&source, &batch).unwrap();
        for (run, scans) in demo.peaks.clone() {
            attach_scans(&mut project, &run, scans.into_iter().map(Ok)).unwrap();
        }
        crate::validate_project(&mut project, &Settings::default()).unwrap();

        let runs = project.mbr_inputs().unwrap();
        let acceptor = runs.iter().find(|r| r.name == "control-A").unwrap().clone();
        let scans = demo
            .peaks
            .iter()
            .find(|(run, _)| run == "control-A")
            .map(|(_, scans)| scans.clone())
            .unwrap();
        let report = match_run(&mut project, &runs, &acceptor, scans.into_iter().map(Ok)).unwrap();
        assert_eq!((report.aligned, report.donors), (3, 3));
        assert_eq!(report.ms1_scans_damaged, 0);
        // The run's own precursors are found where the others predict them.
        // (Its decoys, looked for in clean synthetic scans, mostly find
        // nothing, too few to train on: ranked by retention-time error.)
        assert!(report.summary.known > 100, "{:?}", report.summary);
        drop(project);

        // Every target looked for, and whether it won its pair (q given).
        let looked: Vec<(String, u8, Option<f64>)> = rusqlite::Connection::open(&path)
            .unwrap()
            .prepare(
                "SELECT pe.sequence, t.charge, t.q_value FROM transfers t
                 JOIN peptides pe ON pe.id = t.peptide_id
                 JOIN ms_runs r ON r.id = t.ms_run_id WHERE r.name = 'control-A'",
            )
            .unwrap()
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        let is_dropped = |l: &&(String, u8, Option<f64>)| dropped.contains(&(l.0.clone(), l.1));
        let looked_dropped = looked.iter().filter(is_dropped).count();
        let won: Vec<_> = looked.iter().filter(|l| l.2.is_some()).collect();
        let right = won.iter().filter(|l| is_dropped(l)).count();
        let wrong = won.len() - right;
        // 34 of 42 without the support floors, 33 with them (they set the
        // thinnest few aside).
        assert!(
            right * 4 >= looked_dropped * 3,
            "{right} of {looked_dropped} recovered"
        );
        assert!(
            wrong <= report.summary.decoys_won + 1,
            "{wrong} absent precursors won"
        );
    }

    #[test]
    fn donor_intensities_are_put_on_this_runs_scale() {
        // Both runs pass precursors 1..=60 (10 s apart, the same times); each
        // reads 1e6 here and 1e7 in the donor, so the donor is 10x brighter
        // (log10 ratio -1). The donor also passes precursor 1000 at 5e7, so
        // here it is expected at 5e6 (log10 6.69897); its decoy shares that.
        // Precursor 1001 has no peak in the donor: no expected intensity.
        // A second donor shares peaks with this run for only 29 precursors
        // (one short of MIN_KNOWN), so it has no scale and adds nothing,
        // though it reads precursor 1000 at 1e9.
        let precursor = |id: i64, apex: Option<f64>| tuffold_store::MbrPrecursor {
            peptide_id: id,
            charge: 2,
            mz: 500.0 + id as f64,
            rt: 10.0 * (id % 100) as f64,
            mobility: None,
            passing: true,
            apex_intensity: apex,
        };
        let here = MbrRun {
            id: 1,
            name: "here".into(),
            precursors: (1..=60).map(|i| precursor(i, Some(1e6))).collect(),
        };
        let mut bright: Vec<_> = (1..=60).map(|i| precursor(i, Some(1e7))).collect();
        bright.push(precursor(1000, Some(5e7)));
        bright.push(precursor(1001, None));
        let dark: Vec<_> = (1..=60)
            .map(|i| precursor(i, (i < 30).then_some(1e8)))
            .chain([precursor(1000, Some(1e9))])
            .collect();
        let donors = [
            here.clone(),
            MbrRun {
                id: 2,
                name: "bright".into(),
                precursors: bright,
            },
            MbrRun {
                id: 3,
                name: "dark".into(),
                precursors: dark,
            },
        ];
        let (looked, aligned, _) = plan(&here, &donors);
        assert_eq!(aligned, 2);
        let of = |id: i64| looked.iter().filter(move |l| l.peptide_id == id);
        assert_eq!(of(1000).count(), 2);
        for l in of(1000) {
            assert!((l.intensity.unwrap() - 5e6_f64.log10()).abs() < 1e-9);
        }
        assert_eq!(of(1001).count(), 2);
        assert!(of(1001).all(|l| l.intensity.is_none()));
        // A known transfer comes back to its own intensity: 1e7 / 10.
        let known = of(30).find(|l| l.kind == TransferKind::Known).unwrap();
        assert!((known.intensity.unwrap() - 6.0).abs() < 1e-9);
    }

    #[test]
    fn transfers_are_looked_for_at_the_donors_mobility_shifted_onto_this_run() {
        // Both runs pass precursors 1..=60 (10 s apart, the same times in
        // both); the donor's mobilities read 0.004 lower than this run's,
        // and it also passes precursor 1000 at 1/K0 0.950.
        let precursor = |id: i64, rt: f64, mobility: f64| tuffold_store::MbrPrecursor {
            peptide_id: id,
            charge: 2,
            mz: 500.0 + id as f64,
            rt,
            mobility: Some(mobility),
            passing: true,
            apex_intensity: None,
        };
        let shared = |shift: f64| {
            (1..=60).map(move |i| precursor(i, 10.0 * i as f64, 0.9 + 0.005 * i as f64 + shift))
        };
        let here = MbrRun {
            id: 1,
            name: "here".into(),
            precursors: shared(0.0).collect(),
        };
        let mut donor_precursors: Vec<_> = shared(-0.004).collect();
        donor_precursors.push(precursor(1000, 305.0, 0.950));
        let donor = MbrRun {
            id: 2,
            name: "donor".into(),
            precursors: donor_precursors,
        };
        let (looked, aligned, _) = plan(&here, &[here.clone(), donor]);
        assert_eq!(aligned, 1);
        // The transfer and its decoy: 0.950 + 0.004.
        let transfer: Vec<&Looked> = looked.iter().filter(|l| l.peptide_id == 1000).collect();
        assert_eq!(transfer.len(), 2);
        for l in transfer {
            assert!(
                (l.mobility.unwrap() - 0.954).abs() < 1e-9,
                "{:?}",
                l.mobility
            );
        }
        // A known transfer comes back to its own mobility: 0.9 + 0.005 × 30.
        let known = looked
            .iter()
            .find(|l| l.peptide_id == 30 && l.kind == TransferKind::Known)
            .unwrap();
        assert!((known.mobility.unwrap() - 1.05).abs() < 1e-9);
    }

    fn peak(apex_rt: f64, apex_intensity: f64, points: usize) -> Feature {
        Feature {
            area: apex_intensity * 10.0,
            apex_rt,
            apex_intensity,
            start_rt: apex_rt - 5.0,
            end_rt: apex_rt + 5.0,
            points,
            apex_ppm: 0.0,
            // A 998 Da precursor's expected M+1 / M (see `M1_PER_DALTON`).
            apex_ratio: (500.0 - PROTON) * 2.0 * M1_PER_DALTON,
        }
    }

    #[test]
    fn support_floors_are_the_known_peaks_low_end() {
        // 100 known peaks with apex intensities 100, 200, … 10,000 and 1, 2,
        // … 100 scans. At 5% the floor is the 6th smallest of each (five
        // peaks lie below it): 600 and 6.
        let known: Vec<Feature> = (1..=100)
            .map(|i| peak(100.0, 100.0 * i as f64, i))
            .collect();
        let support = Support::of_known(&known);
        assert_eq!(
            support,
            Support {
                apex_intensity: 600.0,
                points: 6
            }
        );
        assert!(support.holds(&peak(100.0, 600.0, 6)));
        assert!(!support.holds(&peak(100.0, 599.0, 6)));
        assert!(!support.holds(&peak(100.0, 600.0, 5)));
        // Too few known peaks to set a floor: everything holds.
        let few = Support::of_known(&known[..MIN_KNOWN - 1]);
        assert_eq!(few, Support::default());
        assert!(few.holds(&peak(100.0, 1.0, 3)));
    }

    #[test]
    fn a_thin_winner_is_not_used_but_a_strong_one_is() {
        // 40 known precursors found where predicted (apex 1,000,000 + 1,000 k
        // and 20 + k % 5 scans), each with a known decoy found 300 s away
        // (so the model learns that retention-time error means a wrong peak).
        // 150 targets are found at their predicted time with strong peaks
        // (1.5 million, 22 scans) and no decoy finds anything, so each wins
        // its pair and q = 1/150 ≤ 1%. One of them is thin (5,000, 2 scans):
        // the floors from the known peaks (1,002,000 and 20: ⌊40 · 0.05⌋ = 2
        // smaller ones) put it out, so 149 are accepted, not 150.
        let looked_row = |pair: usize, kind: TransferKind, at: f64| Looked {
            pair,
            kind,
            peptide_id: pair as i64,
            charge: 2,
            mz: 500.0,
            predicted: at,
            at,
            spread: 10.0,
            mobility: None,
            intensity: None,
        };
        let mut looked = Vec::new();
        let mut features = HashMap::new();
        let mut add = |row: Looked, found: Option<Feature>| {
            features.insert(looked.len(), found);
            looked.push(row);
        };
        for k in 0..40 {
            let (intensity, points) = (1_000_000.0 + 1_000.0 * k as f64, 20 + k % 5);
            let at = 100.0 + k as f64;
            add(
                looked_row(1000 + k, TransferKind::Known, at),
                Some(peak(at + 1.0, intensity, points)),
            );
            add(
                looked_row(1000 + k, TransferKind::KnownDecoy, at),
                Some(peak(at + 300.0, intensity, points)),
            );
        }
        for t in 0..150 {
            let at = 500.0 + t as f64;
            let strong = t != 7;
            let found = peak(
                at + 1.0,
                if strong { 1.5e6 } else { 5_000.0 },
                if strong { 22 } else { 2 },
            );
            add(looked_row(t, TransferKind::Target, at), Some(found));
            add(looked_row(t, TransferKind::Decoy, at + 300.0), None);
        }
        let (transfers, summary) = judged(&looked, &features);
        assert_eq!((summary.candidates, summary.accepted), (150, 149));
        let thin = transfers.iter().find(|t| t.peptide_id == 7).unwrap();
        // Measured, won its pair, and set aside.
        assert!(thin.area.is_some() && thin.q_value.is_none());
        assert!(
            transfers
                .iter()
                .filter(|t| t.peptide_id != 7)
                .all(|t| t.q_value.is_some_and(|q| q <= TRANSFER_FDR))
        );
    }

    #[test]
    fn expected_isotope_ratio_follows_averagine() {
        // A 1,500 Da peptide: 1500 × 0.060185 / 111.1254 = 0.8124.
        assert!((1500.0 * M1_PER_DALTON - 0.8124).abs() < 1e-4);
    }
}
