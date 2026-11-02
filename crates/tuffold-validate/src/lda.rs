//! Fisher's linear discriminant, the model the rescorer trains.
//!
//! Features are standardized on the training rows; the weights solve
//! `(S + ridge * I) w = mean(positives) - mean(negatives)`, where `S` is the
//! pooled within-class covariance. The ridge keeps correlated engine scores
//! (e.g. `deltacn` and `deltacnstar`) from making `S` singular.

use crate::features::Features;

const RIDGE: f64 = 0.01;
/// Columns whose training standard deviation is below this are ignored.
const MIN_SPREAD: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Lda {
    mean: Vec<f64>,
    /// `1 / standard deviation`, or 0 for ignored columns.
    inverse_spread: Vec<f64>,
    weights: Vec<f64>,
}

impl Lda {
    /// `None` when a class is empty or the system cannot be solved.
    pub fn fit(features: &Features, positives: &[usize], negatives: &[usize]) -> Option<Self> {
        if positives.is_empty() || negatives.is_empty() {
            return None;
        }
        let rows: Vec<usize> = positives.iter().chain(negatives).copied().collect();
        let (mean, inverse_spread) = standardization(features, &rows);
        let active: Vec<usize> = (0..features.columns.len())
            .filter(|&column| inverse_spread[column] > 0.0)
            .collect();
        let standardized = |column: usize, row: usize| {
            (features.value(column, row) - mean[column]) * inverse_spread[column]
        };
        let class_mean = |members: &[usize]| -> Vec<f64> {
            active
                .iter()
                .map(|&column| {
                    let total: f64 = members.iter().map(|&row| standardized(column, row)).sum();
                    total / members.len() as f64
                })
                .collect()
        };
        let (positive_mean, negative_mean) = (class_mean(positives), class_mean(negatives));
        let size = active.len();
        let mut scatter = vec![0.0; size * size];
        for (members, center) in [(positives, &positive_mean), (negatives, &negative_mean)] {
            for &row in members {
                let centered: Vec<f64> = active
                    .iter()
                    .zip(center)
                    .map(|(&column, mean)| standardized(column, row) - mean)
                    .collect();
                for i in 0..size {
                    for j in 0..=i {
                        scatter[i * size + j] += centered[i] * centered[j];
                    }
                }
            }
        }
        let degrees = (rows.len().saturating_sub(2)).max(1) as f64;
        for i in 0..size {
            for j in 0..=i {
                scatter[i * size + j] /= degrees;
                scatter[j * size + i] = scatter[i * size + j];
            }
            scatter[i * size + i] += RIDGE;
        }
        let difference: Vec<f64> = positive_mean
            .iter()
            .zip(&negative_mean)
            .map(|(p, n)| p - n)
            .collect();
        let solved = cholesky_solve(scatter, difference, size)?;
        let mut weights = vec![0.0; features.columns.len()];
        for (&column, weight) in active.iter().zip(solved) {
            weights[column] = weight;
        }
        Some(Self {
            mean,
            inverse_spread,
            weights,
        })
    }

    pub fn score(&self, features: &Features, row: usize) -> f64 {
        (0..self.weights.len())
            .map(|column| {
                self.weights[column]
                    * (features.value(column, row) - self.mean[column])
                    * self.inverse_spread[column]
            })
            .sum()
    }
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

/// Solves `a x = b` for a symmetric positive definite `a` (row-major,
/// `size` x `size`); `None` when `a` is not positive definite.
pub(crate) fn cholesky_solve(mut a: Vec<f64>, b: Vec<f64>, size: usize) -> Option<Vec<f64>> {
    // Factor in place: the lower triangle becomes L with a = L L^T.
    for j in 0..size {
        let diagonal = a[j * size + j] - (0..j).map(|k| a[j * size + k].powi(2)).sum::<f64>();
        if diagonal <= 0.0 || !diagonal.is_finite() {
            return None;
        }
        let pivot = diagonal.sqrt();
        a[j * size + j] = pivot;
        for i in j + 1..size {
            let dot: f64 = (0..j).map(|k| a[i * size + k] * a[j * size + k]).sum();
            a[i * size + j] = (a[i * size + j] - dot) / pivot;
        }
    }
    let mut y = vec![0.0; size];
    for i in 0..size {
        let dot: f64 = (0..i).map(|k| a[i * size + k] * y[k]).sum();
        y[i] = (b[i] - dot) / a[i * size + i];
    }
    let mut x = vec![0.0; size];
    for i in (0..size).rev() {
        let dot: f64 = (i + 1..size).map(|k| a[k * size + i] * x[k]).sum();
        x[i] = (y[i] - dot) / a[i * size + i];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cholesky_solves_a_small_system() {
        // [[4, 2], [2, 3]] x = [2, 1]: det 8, inverse [[3, -2], [-2, 4]] / 8,
        // so x = [(6 - 2) / 8, (-4 + 4) / 8] = [0.5, 0].
        let x = cholesky_solve(vec![4.0, 2.0, 2.0, 3.0], vec![2.0, 1.0], 2).unwrap();
        assert!((x[0] - 0.5).abs() < 1e-12 && x[1].abs() < 1e-12, "{x:?}");
        assert!(cholesky_solve(vec![1.0, 2.0, 2.0, 1.0], vec![1.0, 1.0], 2).is_none());
    }

    #[test]
    fn the_discriminant_weights_the_separating_feature() {
        // Column 0 separates the classes; column 1 has the same mean in both
        // and no within-class covariance with column 0 (centered column 0 is
        // [-1, 0, 1] per class, and [-1, 0, 1] . [1, -2, 1] = 0), so it gets no
        // weight; column 2 is constant and ignored.
        let features = Features {
            names: vec!["a".into(), "b".into(), "c".into()],
            engine_scores: 3,
            columns: vec![
                vec![5.0, 6.0, 7.0, 0.0, 1.0, 2.0],
                vec![1.0, -2.0, 1.0, 1.0, -2.0, 1.0],
                vec![3.0; 6],
            ],
            rows: 6,
        };
        let lda = Lda::fit(&features, &[0, 1, 2], &[3, 4, 5]).unwrap();
        assert!(lda.weights[0] > 0.0);
        assert!(lda.weights[1].abs() < 1e-9, "{:?}", lda.weights);
        assert_eq!(lda.weights[2], 0.0);
        let lowest_positive = (0..3)
            .map(|r| lda.score(&features, r))
            .fold(f64::MAX, f64::min);
        let highest_negative = (3..6)
            .map(|r| lda.score(&features, r))
            .fold(f64::MIN, f64::max);
        assert!(lowest_positive > highest_negative);
        assert!(Lda::fit(&features, &[], &[3]).is_none());
    }
}
