//! Retention-time error as a rescoring feature.
//!
//! A peptide elutes at a time its sequence largely sets; a wrong match was
//! made to a spectrum taken at a time unrelated to its peptide. Each run gets
//! a model of retention time from sequence, trained on the engine's
//! confident targets, and each PSM's error against it, in units of the run's
//! typical error, becomes a feature. It lifts wrong targets and decoys alike
//! (the prediction depends on the sequence only, and neither kind of wrong
//! match elutes where its sequence says), so decoys keep modelling the
//! errors; unlike pooling, it also penalizes a wrong spectrum on a peptide
//! that is present, which elutes elsewhere.
//!
//! The model is classical: ridge regression of the time on amino-acid
//! counts, length and modification count, then a running median of the
//! observed times along that prediction, which follows the gradient's
//! curve. One point per peptide (its confident spectra's median time), so
//! an abundant peptide smeared over many spectra doesn't dominate. It is
//! cross-fitted on the rescoring folds: no PSM's error comes from a model
//! that saw its spectrum. On PXD006447 (136-min gradient) the prediction is
//! off by a median 4.5–5 min.
//!
//! A [`RetentionPredictor`] (a deep-learned model the caller supplies)
//! replaces the regression when given: it is trained the same way, on the
//! fold's confident peptides of the run, and its predictions go through the
//! same running median. Peptides it can't predict (an unsupported
//! modification) keep the classical model's error, so every row has one.
//!
//! Without a retention time on every PSM, or with too few confident
//! peptides in a run to train on, there is no feature.

use std::collections::{BTreeMap, HashMap};

use crate::input::{PeptideRow, PsmRow};
use crate::lda::cholesky_solve;
use crate::rescore::{FOLDS, Rows, confident_by_best_score, fold_of};

pub const RT_ERROR: &str = "|RT error| (run's typical errors)";
pub(crate) const LOG_RT_ERROR: &str = "log(1 + |RT error|)";
pub const DEEP_RT_ERROR: &str = "|RT error, deep model| (run's typical errors)";
pub(crate) const LOG_DEEP_RT_ERROR: &str = "log(1 + |RT error, deep model|)";

/// Learns one run's retention times from peptides and predicts others'.
pub trait RetentionPredictor: Sync {
    /// Trained on (peptide id, retention time) pairs of one run; `None` when
    /// it can't be.
    fn train(&self, training: &[(i64, f64)]) -> Option<Box<dyn TrainedRetention + '_>>;
}

pub trait TrainedRetention {
    /// A raw prediction per peptide id (any monotone scale); `None` for one
    /// it can't predict.
    fn predict(&self, peptides: &[i64]) -> Vec<Option<f64>>;
}

/// Confident peptides a run needs in each fold's training set: well above
/// the model's 23 coefficients.
const MIN_TRAINING_PEPTIDES: usize = 100;
const RIDGE: f64 = 1.0;
const AMINO_ACIDS: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// Counts per amino acid, any other letter, length, modifications, and an
/// intercept.
const WIDTH: usize = AMINO_ACIDS.len() + 4;

/// Each row's retention-time error in its run's typical errors (parallel to
/// the engine's rows), and whether the deep model made it; `None` when it
/// can't be had.
pub(crate) fn rt_errors(
    rows: &Rows<'_>,
    psms: &[&PsmRow],
    peptides: &[&PeptideRow],
    predictor: Option<&dyn RetentionPredictor>,
) -> Option<(Vec<f64>, bool)> {
    let times: Vec<f64> = psms
        .iter()
        .map(|psm| psm.retention_time.filter(|t| t.is_finite()))
        .collect::<Option<_>>()?;
    let mut confident = vec![false; psms.len()];
    for row in confident_by_best_score(rows) {
        confident[row] = true;
    }
    let mut errors = vec![f64::NAN; psms.len()];
    // The deep model's errors are used only if it trains in every fold and
    // run (a mix of the two scales would rank runs unevenly).
    let mut deep = vec![f64::NAN; psms.len()];
    let mut deep_ok = predictor.is_some();
    let runs: BTreeMap<i64, Vec<usize>> =
        psms.iter()
            .enumerate()
            .fold(BTreeMap::new(), |mut runs, (row, psm)| {
                runs.entry(psm.run_id).or_default().push(row);
                runs
            });
    for fold in 0..FOLDS {
        for members in runs.values() {
            let in_fold = |row: &usize| fold_of(rows.folds[*row]) == fold;
            // One training point per peptide: the median time of its
            // confident spectra outside this fold.
            let mut times_of: BTreeMap<i64, Vec<f64>> = BTreeMap::new();
            for &row in members
                .iter()
                .filter(|row| confident[**row] && !in_fold(row))
            {
                times_of
                    .entry(peptides[row].id)
                    .or_default()
                    .push(times[row]);
            }
            if times_of.len() < MIN_TRAINING_PEPTIDES {
                return None;
            }
            let examples: HashMap<i64, &PeptideRow> = members
                .iter()
                .map(|&row| (peptides[row].id, peptides[row]))
                .collect();
            let peptide_times: Vec<(i64, f64)> = times_of
                .into_iter()
                .map(|(id, mut times)| (id, median(&mut times)))
                .collect();
            let training: Vec<(Vec<f64>, f64)> = peptide_times
                .iter()
                .map(|&(id, time)| (encode(examples[&id]), time))
                .collect();
            let model = Model::fit(&training)?;
            let fold_rows: Vec<usize> =
                members.iter().copied().filter(|row| in_fold(row)).collect();
            for &row in &fold_rows {
                let predicted = model.predict(&encode(peptides[row]));
                errors[row] = (times[row] - predicted).abs() / model.typical_error;
            }
            if deep_ok {
                deep_ok = predictor.is_some_and(|predictor| {
                    deep_errors(
                        predictor,
                        &peptide_times,
                        &fold_rows,
                        peptides,
                        &times,
                        &mut deep,
                    )
                });
            }
        }
    }
    if deep_ok {
        for (deep, classical) in deep.iter_mut().zip(&errors) {
            if deep.is_nan() {
                *deep = *classical;
            }
        }
        errors = deep;
    }
    errors
        .iter()
        .all(|e| e.is_finite())
        .then_some((errors, deep_ok))
}

/// Overwrite the fold's errors with the deep model's where it predicts;
/// `false` when it can't be trained (the classical errors stand).
fn deep_errors(
    predictor: &dyn RetentionPredictor,
    peptide_times: &[(i64, f64)],
    fold_rows: &[usize],
    peptides: &[&PeptideRow],
    times: &[f64],
    errors: &mut [f64],
) -> bool {
    let Some(trained) = predictor.train(peptide_times) else {
        return false;
    };
    let ids: Vec<i64> = peptide_times.iter().map(|(id, _)| *id).collect();
    let points: Vec<(f64, f64)> = trained
        .predict(&ids)
        .into_iter()
        .zip(peptide_times)
        .filter_map(|(raw, &(_, time))| Some((raw?, time)))
        .collect();
    if points.len() < MIN_TRAINING_PEPTIDES {
        return false;
    }
    let Some(curve) = Curve::fit(points) else {
        return false;
    };
    let wanted: Vec<i64> = fold_rows.iter().map(|&row| peptides[row].id).collect();
    for (&row, raw) in fold_rows.iter().zip(trained.predict(&wanted)) {
        if let Some(raw) = raw {
            errors[row] = (times[row] - curve.at(raw)).abs() / curve.typical_error;
        }
    }
    true
}

/// A run's retention-time model: a ridge regression on [`encode`]d
/// sequences, then a running median of the observed times along it.
pub(crate) struct Model {
    weights: Vec<f64>,
    curve: Curve,
    /// Median absolute error on the training peptides.
    pub typical_error: f64,
}

/// Observed times as a function of a raw prediction: the running median of
/// the training peptides' times along it, which follows the gradient.
pub(crate) struct Curve {
    /// Raw predictions, ascending, and the running median of their times.
    points: Vec<(f64, f64)>,
    /// Median absolute error on the training peptides.
    pub typical_error: f64,
}

impl Curve {
    /// `None` without an error scale (a perfect fit).
    pub fn fit(mut points: Vec<(f64, f64)>) -> Option<Self> {
        points.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        let half = (points.len() / 50).max(10);
        let running: Vec<(f64, f64)> = (0..points.len())
            .map(|i| {
                let window = &points[i.saturating_sub(half)..(i + half + 1).min(points.len())];
                let mut times: Vec<f64> = window.iter().map(|(_, y)| *y).collect();
                (points[i].0, median(&mut times))
            })
            .collect();
        let mut curve = Self {
            points: running,
            typical_error: 0.0,
        };
        let mut residuals: Vec<f64> = points
            .iter()
            .map(|&(x, y)| (y - curve.at(x)).abs())
            .collect();
        curve.typical_error = median(&mut residuals);
        (curve.typical_error > 0.0).then_some(curve)
    }

    /// The curve at `raw`, interpolated, and flat beyond its ends.
    pub fn at(&self, raw: f64) -> f64 {
        let after = self.points.partition_point(|(at, _)| *at < raw);
        match (
            after.checked_sub(1).map(|i| self.points[i]),
            self.points.get(after),
        ) {
            (Some((x0, y0)), Some(&(x1, y1))) if x1 > x0 => y0 + (y1 - y0) * (raw - x0) / (x1 - x0),
            (_, Some(&(_, y))) | (Some((_, y)), None) => y,
            (None, None) => f64::NAN,
        }
    }
}

impl Model {
    /// `None` when the regression can't be solved or the fit is perfect
    /// (no error scale).
    pub fn fit(training: &[(Vec<f64>, f64)]) -> Option<Self> {
        let mut gram = vec![0.0; WIDTH * WIDTH];
        let mut moment = vec![0.0; WIDTH];
        for (x, y) in training {
            for i in 0..WIDTH {
                moment[i] += x[i] * y;
                for j in 0..WIDTH {
                    gram[i * WIDTH + j] += x[i] * x[j];
                }
            }
        }
        for i in 0..WIDTH {
            gram[i * WIDTH + i] += RIDGE;
        }
        let weights = cholesky_solve(gram, moment, WIDTH)?;
        let linear = |x: &[f64]| x.iter().zip(&weights).map(|(a, b)| a * b).sum::<f64>();
        let curve = Curve::fit(training.iter().map(|(x, y)| (linear(x), *y)).collect())?;
        Some(Self {
            typical_error: curve.typical_error,
            weights,
            curve,
        })
    }

    /// Predicted retention time: the curve at the linear prediction.
    pub fn predict(&self, x: &[f64]) -> f64 {
        let linear: f64 = x.iter().zip(&self.weights).map(|(a, b)| a * b).sum();
        self.curve.at(linear)
    }
}

/// Amino-acid counts, other letters, length, modifications, intercept.
pub(crate) fn encode(peptide: &PeptideRow) -> Vec<f64> {
    let mut x = vec![0.0; WIDTH];
    for residue in peptide.sequence.bytes() {
        let slot = AMINO_ACIDS
            .iter()
            .position(|&aa| aa == residue)
            .unwrap_or(AMINO_ACIDS.len());
        x[slot] += 1.0;
    }
    x[AMINO_ACIDS.len() + 1] = peptide.sequence.len() as f64;
    x[AMINO_ACIDS.len() + 2] = f64::from(peptide.modifications);
    x[AMINO_ACIDS.len() + 3] = 1.0;
    x
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    let middle = values.len() / 2;
    if values.len() % 2 == 1 {
        values[middle]
    } else {
        (values[middle - 1] + values[middle]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peptide(sequence: &str) -> PeptideRow {
        PeptideRow {
            id: 0,
            sequence: sequence.to_string(),
            modifications: 0,
            is_decoy: false,
        }
    }

    #[test]
    fn sequences_are_encoded_as_counts() {
        let x = encode(&PeptideRow {
            modifications: 2,
            ..peptide("AAKXR")
        });
        assert_eq!(x[0], 2.0, "A");
        assert_eq!(x[8], 1.0, "K");
        assert_eq!(x[14], 1.0, "R");
        assert_eq!(x[20], 1.0, "X is another letter");
        assert_eq!(x[21..], [5.0, 2.0, 1.0]);
    }

    #[test]
    fn the_model_follows_a_curved_gradient() {
        // Time is a curved function of hydrophobicity (counts of L), with
        // the gradient flattening at the end: t = 100 · (1 − e^(−n/4)).
        let training: Vec<(Vec<f64>, f64)> = (0..200)
            .map(|i| {
                let leucines = i % 20;
                let sequence = format!("{}{}K", "L".repeat(leucines), "G".repeat(5));
                let time = 100.0 * (1.0 - (-(leucines as f64) / 4.0).exp());
                // Each sequence appears 10 times, half at +1 and half at
                // −1, so the fit has an error scale of 1.
                let noise = if (i / 20) % 2 == 0 { 1.0 } else { -1.0 };
                (encode(&peptide(&sequence)), time + noise)
            })
            .collect();
        let model = Model::fit(&training).expect("fits");
        for leucines in [2, 8, 15] {
            let sequence = format!("{}{}K", "L".repeat(leucines), "G".repeat(5));
            let truth = 100.0 * (1.0 - (-(leucines as f64) / 4.0).exp());
            let predicted = model.predict(&encode(&peptide(&sequence)));
            assert!(
                (predicted - truth).abs() < 2.0,
                "{leucines}: {predicted} vs {truth}"
            );
        }
        assert!(
            (0.5..=1.5).contains(&model.typical_error),
            "{}",
            model.typical_error
        );
    }
}
