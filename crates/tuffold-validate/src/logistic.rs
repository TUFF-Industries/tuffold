//! L2-regularized logistic regression, the model the rescorer trains.
//!
//! Features are standardized on the training rows; the weights minimize
//! `Σ log(1 + e^(−y·(w·x + b))) + ridge/2 · |w|²` (the intercept `b` is not
//! penalized), found by Newton's method with a backtracking line search, so
//! the fit is deterministic. The ridge keeps weights finite when the classes
//! separate completely and correlated engine scores (`deltacn`,
//! `deltacnstar`) from making the Hessian singular.
//!
//! The rescorer trains this and Fisher's discriminant (`lda.rs`) and keeps
//! whichever accepts more: trained on only the confident targets against
//! the decoys, the discriminant's positives are cut off along the very
//! score being learned, which distorts the within-class covariance it
//! relies on (on real Sage results it barely beat the best single score,
//! where this matched mokapot's linear SVM), but on Gaussian classes with
//! equal covariance, its ideal case, it does better.

use crate::features::Features;
use crate::lda::cholesky_solve;

/// Ridge weight, in units of one PSM's loss (scikit-learn's `C = 1`).
const RIDGE: f64 = 1.0;
/// Keeps the Hessian positive definite along the unpenalized intercept.
const INTERCEPT_RIDGE: f64 = 1e-9;
/// Columns whose training standard deviation is below this are ignored.
const MIN_SPREAD: f64 = 1e-12;
const MAX_NEWTON_STEPS: usize = 100;
const MAX_HALVINGS: usize = 40;
/// Stop when no weight moves more than this (standardized units).
const STEP_TOLERANCE: f64 = 1e-9;
/// Armijo constant for the line search.
const SUFFICIENT_DECREASE: f64 = 1e-4;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Logistic {
    mean: Vec<f64>,
    /// `1 / standard deviation`, or 0 for ignored columns.
    inverse_spread: Vec<f64>,
    weights: Vec<f64>,
    intercept: f64,
}

impl Logistic {
    /// `None` when a class is empty or Newton's method breaks down.
    pub fn fit(features: &Features, positives: &[usize], negatives: &[usize]) -> Option<Self> {
        if positives.is_empty() || negatives.is_empty() {
            return None;
        }
        let rows: Vec<usize> = positives.iter().chain(negatives).copied().collect();
        let (mean, inverse_spread) = standardization(features, &rows);
        let active: Vec<usize> = (0..features.columns.len())
            .filter(|&column| inverse_spread[column] > 0.0)
            .collect();
        // Row-major standardized design, with a trailing 1 for the intercept.
        let width = active.len() + 1;
        let mut design = Vec::with_capacity(rows.len() * width);
        for &row in &rows {
            design.extend(active.iter().map(|&column| {
                (features.value(column, row) - mean[column]) * inverse_spread[column]
            }));
            design.push(1.0);
        }
        let labels: Vec<f64> = (0..rows.len())
            .map(|i| if i < positives.len() { 1.0 } else { 0.0 })
            .collect();
        let solved = newton(&design, &labels, width)?;
        let mut weights = vec![0.0; features.columns.len()];
        for (&column, weight) in active.iter().zip(&solved) {
            weights[column] = *weight;
        }
        Some(Self {
            mean,
            inverse_spread,
            weights,
            intercept: solved[width - 1],
        })
    }

    pub fn score(&self, features: &Features, row: usize) -> f64 {
        self.intercept
            + (0..self.weights.len())
                .map(|column| {
                    self.weights[column]
                        * (features.value(column, row) - self.mean[column])
                        * self.inverse_spread[column]
                })
                .sum::<f64>()
    }
}

/// Minimizes the penalized loss over `design` (row-major, `width` columns,
/// the last one the intercept). `None` if the Hessian is not positive
/// definite or the loss stops being finite.
fn newton(design: &[f64], labels: &[f64], width: usize) -> Option<Vec<f64>> {
    let mut beta = vec![0.0; width];
    let mut current = loss(design, labels, &beta, width);
    for _ in 0..MAX_NEWTON_STEPS {
        let (gradient, hessian) = derivatives(design, labels, &beta, width);
        let step = cholesky_solve(hessian, gradient.clone(), width)?;
        let slope: f64 = gradient.iter().zip(&step).map(|(g, s)| g * s).sum();
        let mut scale = 1.0;
        let mut accepted = None;
        for _ in 0..MAX_HALVINGS {
            let trial: Vec<f64> = beta.iter().zip(&step).map(|(b, s)| b - scale * s).collect();
            let value = loss(design, labels, &trial, width);
            if value.is_finite() && value <= current - SUFFICIENT_DECREASE * scale * slope {
                accepted = Some((trial, value));
                break;
            }
            scale /= 2.0;
        }
        // No decrease along a Newton direction: already at the minimum.
        let Some((trial, value)) = accepted else {
            break;
        };
        let moved = beta
            .iter()
            .zip(&trial)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        beta = trial;
        current = value;
        if moved < STEP_TOLERANCE {
            break;
        }
    }
    current.is_finite().then_some(beta)
}

fn penalty(beta: &[f64], width: usize) -> f64 {
    beta[..width - 1].iter().map(|w| w * w).sum::<f64>() * RIDGE / 2.0
}

/// Penalized negative log-likelihood.
fn loss(design: &[f64], labels: &[f64], beta: &[f64], width: usize) -> f64 {
    let data: f64 = design
        .chunks_exact(width)
        .zip(labels)
        .map(|(x, &label)| {
            let z: f64 = x.iter().zip(beta).map(|(a, b)| a * b).sum();
            // log(1 + e^(−z)) for positives, log(1 + e^z) for negatives,
            // computed without overflow.
            let margin = if label > 0.5 { z } else { -z };
            softplus(-margin)
        })
        .sum();
    data + penalty(beta, width)
}

/// `log(1 + e^x)` without overflow.
fn softplus(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

/// Gradient and Hessian (row-major) of [`loss`].
fn derivatives(design: &[f64], labels: &[f64], beta: &[f64], width: usize) -> (Vec<f64>, Vec<f64>) {
    let mut gradient = vec![0.0; width];
    let mut hessian = vec![0.0; width * width];
    for (x, &label) in design.chunks_exact(width).zip(labels) {
        let z: f64 = x.iter().zip(beta).map(|(a, b)| a * b).sum();
        let p = 1.0 / (1.0 + (-z).exp());
        let residual = p - label;
        let curvature = p * (1.0 - p);
        for i in 0..width {
            gradient[i] += residual * x[i];
            for j in 0..=i {
                hessian[i * width + j] += curvature * x[i] * x[j];
            }
        }
    }
    for i in 0..width {
        for j in 0..i {
            hessian[j * width + i] = hessian[i * width + j];
        }
        let ridge = if i + 1 == width {
            INTERCEPT_RIDGE
        } else {
            gradient[i] += RIDGE * beta[i];
            RIDGE
        };
        hessian[i * width + i] += ridge;
    }
    (gradient, hessian)
}

/// Per-column mean and inverse standard deviation over `rows`.
fn standardization(features: &Features, rows: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let count = rows.len() as f64;
    features
        .columns
        .iter()
        .map(|column| {
            let mean = rows.iter().map(|&row| column[row]).sum::<f64>() / count;
            let variance = rows
                .iter()
                .map(|&row| (column[row] - mean).powi(2))
                .sum::<f64>()
                / count;
            let spread = variance.sqrt();
            (
                mean,
                if spread > MIN_SPREAD {
                    1.0 / spread
                } else {
                    0.0
                },
            )
        })
        .unzip()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn features(columns: Vec<Vec<f64>>) -> Features {
        let rows = columns[0].len();
        Features {
            names: (0..columns.len()).map(|i| format!("f{i}")).collect(),
            engine_scores: columns.len(),
            columns,
            rows,
        }
    }

    #[test]
    fn the_model_weights_the_separating_feature() {
        // Column 0 separates the classes completely; column 1 is the same
        // pattern in both classes (no information); column 2 is constant.
        let f = features(vec![
            vec![5.0, 6.0, 7.0, 0.0, 1.0, 2.0],
            vec![1.0, -2.0, 1.0, 1.0, -2.0, 1.0],
            vec![3.0; 6],
        ]);
        let model = Logistic::fit(&f, &[0, 1, 2], &[3, 4, 5]).unwrap();
        assert!(model.weights[0] > 0.0);
        // By symmetry the uninformative column's weight is exactly 0 at the
        // optimum; Newton stops within its tolerance.
        assert!(model.weights[1].abs() < 1e-6, "{:?}", model.weights);
        assert_eq!(model.weights[2], 0.0);
        // Separable classes: the ridge keeps the weight finite.
        assert!(model.weights[0].is_finite() && model.weights[0] < 100.0);
        let lowest_positive = (0..3).map(|r| model.score(&f, r)).fold(f64::MAX, f64::min);
        let highest_negative = (3..6).map(|r| model.score(&f, r)).fold(f64::MIN, f64::max);
        assert!(lowest_positive > highest_negative);
        assert!(Logistic::fit(&f, &[], &[3]).is_none());
    }

    #[test]
    fn overlapping_classes_reach_the_hand_derived_optimum() {
        // One feature, x = 0 three times (1 positive) and x = 1 three times
        // (2 positives). Standardized, x is −1 or +1. Swapping x and the
        // labels maps the data onto itself, so the intercept is 0 and the
        // score is w·x. The gradient in w is
        //   [3σ(w) − 2] − [3σ(−w) − 1] + w = 6σ(w) − 4 + w = 0
        // (σ(−w) = 1 − σ(w); the ridge adds w). Bisection by hand:
        // w = 0.4032 gives 3.99989, w = 0.4033 gives 4.00018, so w ≈ 0.40324
        // and the score gap between x = 1 and x = 0 is 2w ≈ 0.8065.
        let f = features(vec![vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]);
        let model = Logistic::fit(&f, &[2, 4, 5], &[0, 1, 3]).unwrap();
        let gap = model.score(&f, 3) - model.score(&f, 0);
        assert!((gap - 0.8065).abs() < 1e-3, "{gap}");
        assert!(model.intercept.abs() < 1e-9, "{}", model.intercept);
    }
}
