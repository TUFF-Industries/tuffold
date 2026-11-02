//! How well a PSM's measured fragments match the model's prediction for its
//! peptide: rescoring features that need no training on the project, so
//! they are decoy-safe by construction (a decoy's prediction is as good as
//! a target's; only the spectrum decides how well it fits).
//!
//! Collision energy and instrument aren't in the project, so they are
//! calibrated per project: the pair whose predictions fit a sample of rank-1
//! PSMs best (90th percentile of the Pearson correlation, so the wrong
//! matches in it don't decide it). The caller picks the sample: the engines'
//! best-scoring PSMs, targets and decoys alike (`tuffold-engine`'s `deep`).
//! When even the best pair fits poorly, there are no features.
//!
//! Fragments are matched within 20 ppm (Orbitrap and TOF MS/MS) or 0.5 Da
//! (ion-trap MS/MS), whichever fits the sample better. A wide window also
//! catches unrelated peaks, which lifts every fit a little, so 0.5 Da is
//! used only when it fits clearly better: on high-resolution projects 20 ppm
//! fits better by 0.004–0.02, on ion-trap ones 0.5 Da by 0.29–0.50.

use std::collections::{BTreeMap, HashMap, HashSet};

use tuffold_core::PeptideForm;
use tuffold_core::fragments::{IonSeries, fragment_ions};

use crate::encode::Encoded;
use crate::ms2::{CHANNELS, Instrument, Ms2Model};

pub const FEATURE_NAMES: [&str; 5] = [
    "predicted spectrum: spectral angle",
    "predicted spectrum: Pearson (√intensity)",
    "predicted spectrum: share of top-6 fragments seen",
    "predicted spectrum: b-ion Pearson",
    "predicted spectrum: y-ion Pearson",
];

/// How far a fragment's peak may be from its m/z.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    Ppm(f64),
    Da(f64),
}

impl std::fmt::Display for Tolerance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Tolerance::Ppm(ppm) => write!(f, "{ppm} ppm"),
            Tolerance::Da(da) => write!(f, "{da} Da"),
        }
    }
}

const HIGH_RESOLUTION: Tolerance = Tolerance::Ppm(20.0);
const LOW_RESOLUTION: Tolerance = Tolerance::Da(0.5);
/// How much better 0.5 Da must fit than 20 ppm to be used.
const LOW_RESOLUTION_MARGIN: f64 = 0.05;
const BATCH: usize = 256;
const NCE_GRID: [f64; 9] = [20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0];
const INSTRUMENTS: [Instrument; 4] = [
    Instrument::QExactive,
    Instrument::Lumos,
    Instrument::TimsTof,
    Instrument::SciexTof,
];
/// Below this (90th-percentile Pearson of the calibration sample) the model
/// doesn't describe the spectra.
pub const MIN_FIT: f64 = 0.6;
const FIT_QUANTILE: f64 = 0.9;
/// Calibration first tries every setting on this many PSMs, then the best
/// [`SHORTLIST`] per tolerance on the whole sample.
const SCREEN_SAMPLE: usize = 500;
const SHORTLIST: usize = 3;

/// Predicted rows by (peptide, precursor charge): a prediction depends on
/// nothing else, so each is made once however many spectra match it.
pub type Predictions = HashMap<(i64, u8), Vec<[f32; CHANNELS]>>;

/// One rank-1 PSM with its spectrum's peaks (sorted by m/z).
pub struct PsmSpectrum<'a> {
    pub psm: i64,
    pub peptide: i64,
    pub form: &'a PeptideForm,
    pub encoded: &'a Encoded,
    pub charge: u8,
    pub mz: &'a [f64],
    pub intensity: &'a [f32],
}

/// The MS2 model at a project's calibrated settings.
#[derive(Clone, Copy)]
pub struct Ms2Scorer<'a> {
    model: &'a Ms2Model,
    pub nce: f64,
    pub instrument: Instrument,
    pub tolerance: Tolerance,
    /// The calibration sample's 90th-percentile Pearson at these settings.
    pub fit: f64,
}

/// The best settings for `sample`, or `None` when none fits.
pub(crate) fn calibrate<'a>(
    model: &'a Ms2Model,
    sample: &[PsmSpectrum<'_>],
) -> Option<Ms2Scorer<'a>> {
    let settings: Vec<Ms2Scorer<'a>> = INSTRUMENTS
        .iter()
        .flat_map(|&instrument| {
            NCE_GRID.iter().map(move |&nce| Ms2Scorer {
                model,
                nce,
                instrument,
                tolerance: HIGH_RESOLUTION,
                fit: 0.0,
            })
        })
        .collect();
    let stride = (sample.len() / SCREEN_SAMPLE).max(1);
    let screen: Vec<&PsmSpectrum<'_>> = sample.iter().step_by(stride).collect();
    let mut screened = Vec::with_capacity(settings.len());
    for scorer in &settings {
        screened.push((scorer, scorer.fits(&screen)?));
    }
    let mut shortlist: Vec<&Ms2Scorer<'a>> = Vec::new();
    for tolerance in 0..2 {
        screened.sort_by(|a, b| b.1[tolerance].total_cmp(&a.1[tolerance]));
        for (scorer, _) in screened.iter().take(SHORTLIST) {
            if !shortlist.iter().any(|known| std::ptr::eq(*known, *scorer)) {
                shortlist.push(scorer);
            }
        }
    }
    let whole: Vec<&PsmSpectrum<'_>> = sample.iter().collect();
    let mut best: [Option<Ms2Scorer<'a>>; 2] = [None, None];
    for scorer in shortlist {
        let fits = scorer.fits(&whole)?;
        for ((&fit, &tolerance), best) in fits.iter().zip(&TOLERANCES).zip(&mut best) {
            if best.is_none_or(|known| fit > known.fit) {
                *best = Some(Ms2Scorer {
                    tolerance,
                    fit,
                    ..*scorer
                });
            }
        }
    }
    let [high, low] = best;
    let chosen = match (high, low) {
        (Some(high), Some(low)) if wide_fits_better(high.fit, low.fit) => low,
        (high, _) => high?,
    };
    (chosen.fit >= MIN_FIT).then_some(chosen)
}

const TOLERANCES: [Tolerance; 2] = [HIGH_RESOLUTION, LOW_RESOLUTION];

fn wide_fits_better(high_resolution: f64, low_resolution: f64) -> bool {
    low_resolution > high_resolution + LOW_RESOLUTION_MARGIN
}

impl Ms2Scorer<'_> {
    /// [`FEATURE_NAMES`] for each PSM the model can predict, reusing and
    /// adding to `predictions`.
    pub fn score(
        &self,
        psms: &[PsmSpectrum<'_>],
        predictions: &mut Predictions,
    ) -> Vec<(i64, Vec<f64>)> {
        let psms: Vec<&PsmSpectrum<'_>> = psms.iter().collect();
        self.predict(&psms, predictions);
        psms.iter()
            .filter_map(|psm| {
                let rows = predictions.get(&(psm.peptide, psm.charge))?;
                Some((psm.psm, features(psm, rows, self.tolerance)?))
            })
            .collect()
    }

    /// The 90th-percentile Pearson at each of [`TOLERANCES`]; `None` when
    /// the model can predict none of `psms`.
    fn fits(&self, psms: &[&PsmSpectrum<'_>]) -> Option<[f64; 2]> {
        let mut predictions = Predictions::new();
        self.predict(psms, &mut predictions);
        let mut fits = [0.0; 2];
        for (fit, &tolerance) in fits.iter_mut().zip(&TOLERANCES) {
            let mut pearsons: Vec<f64> = psms
                .iter()
                .filter_map(|psm| {
                    let rows = predictions.get(&(psm.peptide, psm.charge))?;
                    Some(features(psm, rows, tolerance)?[1])
                })
                .collect();
            if pearsons.is_empty() {
                return None;
            }
            pearsons.sort_by(f64::total_cmp);
            *fit = pearsons[((pearsons.len() - 1) as f64 * FIT_QUANTILE) as usize];
        }
        Some(fits)
    }

    /// Predict every (peptide, charge) of `psms` not yet in `predictions`.
    fn predict(&self, psms: &[&PsmSpectrum<'_>], predictions: &mut Predictions) {
        let mut queued = HashSet::new();
        let mut groups: BTreeMap<(usize, u8), Vec<&PsmSpectrum<'_>>> = BTreeMap::new();
        for psm in psms.iter().filter(|psm| psm.encoded.len() >= 4) {
            let key = (psm.peptide, psm.charge);
            if !predictions.contains_key(&key) && queued.insert(key) {
                groups
                    .entry((psm.encoded.len(), psm.charge))
                    .or_default()
                    .push(psm);
            }
        }
        for ((_, charge), group) in groups {
            for batch in group.chunks(BATCH) {
                let encoded: Vec<&Encoded> = batch.iter().map(|psm| psm.encoded).collect();
                let Ok(predicted) = self
                    .model
                    .predict(&encoded, charge, self.nce, self.instrument)
                else {
                    continue;
                };
                for (psm, rows) in batch.iter().zip(predicted) {
                    predictions.insert((psm.peptide, charge), rows);
                }
            }
        }
    }
}

/// Predicted against observed fragments (the fragment charges
/// `fragment_ions` allows for the precursor).
fn features(
    psm: &PsmSpectrum<'_>,
    rows: &[[f32; CHANNELS]],
    tolerance: Tolerance,
) -> Option<Vec<f64>> {
    let ions = fragment_ions(psm.form, Some(psm.charge))?;
    let length = psm.form.sequence.len();
    let (mut predicted, mut observed, mut is_b) = (Vec::new(), Vec::new(), Vec::new());
    for ion in &ions {
        let ordinal = ion.ordinal as usize;
        let (row, channel) = match ion.series {
            IonSeries::B => (ordinal - 1, usize::from(ion.charge) - 1),
            IonSeries::Y => (length - 1 - ordinal, usize::from(ion.charge) + 1),
        };
        predicted.push(f64::from(rows.get(row)?[channel].max(0.0)));
        observed.push(peak_at(psm.mz, psm.intensity, ion.mz, tolerance));
        is_b.push(ion.series == IonSeries::B);
    }
    let strongest = observed.iter().copied().fold(0.0, f64::max);
    if strongest > 0.0 {
        observed.iter_mut().for_each(|value| *value /= strongest);
    }
    let dot: f64 = predicted.iter().zip(&observed).map(|(p, o)| p * o).sum();
    let norm = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let cosine = (dot / (norm(&predicted) * norm(&observed) + 1e-12)).clamp(-1.0, 1.0);
    let angle = 1.0 - 2.0 * cosine.acos() / std::f64::consts::PI;
    let root = |v: &[f64]| v.iter().map(|x| x.sqrt()).collect::<Vec<f64>>();
    let pick = |v: &[f64], b: bool| -> Vec<f64> {
        v.iter()
            .zip(&is_b)
            .filter(|(_, x)| **x == b)
            .map(|(v, _)| v.sqrt())
            .collect()
    };
    let mut top: Vec<usize> = (0..predicted.len()).collect();
    top.sort_by(|&a, &b| predicted[b].total_cmp(&predicted[a]));
    let top = &top[..top.len().min(6)];
    let seen = top.iter().filter(|&&i| observed[i] > 0.0).count() as f64 / top.len().max(1) as f64;
    Some(vec![
        angle,
        pearson(&root(&predicted), &root(&observed)),
        seen,
        pearson(&pick(&predicted, true), &pick(&observed, true)),
        pearson(&pick(&predicted, false), &pick(&observed, false)),
    ])
}

/// The most intense peak within the tolerance, or 0.
fn peak_at(mz: &[f64], intensity: &[f32], target: f64, tolerance: Tolerance) -> f64 {
    let window = match tolerance {
        Tolerance::Ppm(ppm) => target * ppm * 1e-6,
        Tolerance::Da(da) => da,
    };
    let start = mz.partition_point(|&m| m < target - window);
    mz[start..]
        .iter()
        .zip(&intensity[start..])
        .take_while(|(m, _)| **m <= target + window)
        .map(|(_, i)| f64::from(*i))
        .fold(0.0, f64::max)
}

/// 0 when either side is constant (nothing to correlate).
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 2 {
        return 0.0;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peaks_are_found_within_the_tolerance() {
        let mz = [100.0, 200.0, 200.002, 300.0];
        let intensity = [1.0, 5.0, 7.0, 2.0];
        assert_eq!(peak_at(&mz, &intensity, 200.001, HIGH_RESOLUTION), 7.0);
        assert_eq!(peak_at(&mz, &intensity, 200.3, HIGH_RESOLUTION), 0.0);
        assert_eq!(peak_at(&mz, &intensity, 200.3, LOW_RESOLUTION), 7.0);
        assert_eq!(peak_at(&mz, &intensity, 250.0, LOW_RESOLUTION), 0.0);
        assert!((pearson(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]) - 1.0).abs() < 1e-12);
        assert_eq!(pearson(&[1.0, 1.0], &[2.0, 3.0]), 0.0);
    }

    #[test]
    fn the_wide_window_needs_a_clearly_better_fit() {
        // Fits measured on high-resolution (first two) and ion-trap projects.
        assert!(!wide_fits_better(0.968, 0.948));
        assert!(!wide_fits_better(0.984, 0.980));
        assert!(!wide_fits_better(0.90, 0.94));
        assert!(wide_fits_better(0.646, 0.980));
        assert!(wide_fits_better(0.266, 0.764));
    }
}
