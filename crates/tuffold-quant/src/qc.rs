//! Sample QC: do the samples look like what they're labelled as? Three views
//! of the per-sample log2 values a comparison would use:
//!
//! - **Per sample**: groups measured and their median, for a run that failed
//!   or a channel that was hardly loaded.
//! - **Correlation**: Pearson's r for every pair of samples, over the groups
//!   both measured. Replicates should sit highest. On spectral counts, zeros
//!   and abundance raise every r, so read it relatively.
//! - **PCA**: over the groups measured in every sample, centred per group.
//!   Nothing is imputed: a guessed value would pull samples together or
//!   apart, and PCA is what people judge batch effects by. Samples are few and
//!   groups many, so it comes from the samples' Gram matrix.
//!
//! Samples without any values are left out, as comparisons leave them out;
//! samples without a category stay in, since a mislabelled or unlabelled
//! sample is what this is for.

use crate::Sample;
use crate::anova::GroupValues;

/// Fewer shared groups than this say nothing about a correlation or a PCA.
pub const MIN_GROUPS: usize = 3;
/// Two samples have one principal component and nothing to plot.
pub const MIN_PCA_SAMPLES: usize = 3;

const JACOBI_SWEEPS: usize = 60;

#[derive(Debug, Clone, PartialEq)]
pub struct SampleStats {
    pub name: String,
    pub category: Option<String>,
    /// Groups with a value in this sample.
    pub measured: usize,
    pub median_log2: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Pca {
    /// Groups measured in every sample, which the components are built from.
    pub groups: usize,
    /// Share of the variance on PC1 and PC2.
    pub explained: [f64; 2],
    /// (PC1, PC2) per sample, parallel to [`SampleQc::samples`]; each
    /// component's largest score is positive, so results are repeatable.
    pub scores: Vec<(f64, f64)>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SampleQc {
    pub samples: Vec<SampleStats>,
    /// Pearson's r, `correlation[a][b]`; `None` with fewer than
    /// [`MIN_GROUPS`] shared groups or no variation in one sample.
    pub correlation: Vec<Vec<Option<f64>>>,
    /// `None` with fewer than [`MIN_PCA_SAMPLES`] samples, fewer than
    /// [`MIN_GROUPS`] groups measured in all, or no variance among them.
    pub pca: Option<Pca>,
    /// Groups looked at in all.
    pub groups: usize,
}

/// QC of the samples with data. `with_data` and each group's `log2` are
/// parallel to `samples`.
pub fn sample_qc(samples: &[Sample], with_data: &[bool], groups: &[GroupValues]) -> SampleQc {
    let shown: Vec<usize> = (0..samples.len())
        .filter(|&s| with_data.get(s).copied().unwrap_or(false))
        .collect();
    // Per sample, in `shown` order; rows parallel to `groups`.
    let values: Vec<Vec<Option<f64>>> = groups
        .iter()
        .map(|group| {
            shown
                .iter()
                .map(|&s| group.log2.get(s).copied().flatten())
                .collect()
        })
        .collect();
    let stats = shown
        .iter()
        .enumerate()
        .map(|(column, &s)| {
            let mut measured: Vec<f64> = values.iter().filter_map(|row| row[column]).collect();
            measured.sort_by(f64::total_cmp);
            SampleStats {
                name: samples[s].name.clone(),
                category: samples[s].category.clone(),
                measured: measured.len(),
                median_log2: median(&measured),
            }
        })
        .collect();
    let n = shown.len();
    let mut correlation = vec![vec![None; n]; n];
    for a in 0..n {
        for b in a..n {
            let r = pearson(values.iter().filter_map(|row| Some((row[a]?, row[b]?))));
            correlation[a][b] = r;
            correlation[b][a] = r;
        }
    }
    SampleQc {
        samples: stats,
        correlation,
        pca: pca(&values, n),
        groups: groups.len(),
    }
}

fn median(sorted: &[f64]) -> Option<f64> {
    let middle = sorted.len() / 2;
    match sorted.len() {
        0 => None,
        n if n % 2 == 1 => Some(sorted[middle]),
        _ => Some((sorted[middle - 1] + sorted[middle]) / 2.0),
    }
}

fn pearson(pairs: impl Iterator<Item = (f64, f64)>) -> Option<f64> {
    let (mut n, mut mean_x, mut mean_y) = (0.0, 0.0, 0.0);
    let (mut xx, mut yy, mut xy) = (0.0, 0.0, 0.0);
    // Welford's, so long columns don't lose precision.
    for (x, y) in pairs {
        n += 1.0;
        let (dx, dy) = (x - mean_x, y - mean_y);
        mean_x += dx / n;
        mean_y += dy / n;
        xx += dx * (x - mean_x);
        yy += dy * (y - mean_y);
        xy += dx * (y - mean_y);
    }
    (n >= MIN_GROUPS as f64 && xx > 0.0 && yy > 0.0)
        .then(|| (xy / (xx * yy).sqrt()).clamp(-1.0, 1.0))
}

fn pca(values: &[Vec<Option<f64>>], n: usize) -> Option<Pca> {
    if n < MIN_PCA_SAMPLES {
        return None;
    }
    let complete: Vec<Vec<f64>> = values
        .iter()
        .filter_map(|row| row.iter().copied().collect::<Option<Vec<f64>>>())
        .collect();
    if complete.len() < MIN_GROUPS {
        return None;
    }
    let mut gram = vec![vec![0.0; n]; n];
    for row in &complete {
        let mean = row.iter().sum::<f64>() / n as f64;
        for (a, line) in gram.iter_mut().enumerate() {
            for (b, cell) in line.iter_mut().enumerate() {
                *cell += (row[a] - mean) * (row[b] - mean);
            }
        }
    }
    let (eigenvalues, vectors) = jacobi(gram);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| eigenvalues[b].total_cmp(&eigenvalues[a]));
    let total: f64 = eigenvalues.iter().map(|v| v.max(0.0)).sum();
    if total <= f64::EPSILON {
        return None;
    }
    let component = |rank: usize| -> (f64, Vec<f64>) {
        let at = order[rank];
        let value = eigenvalues[at].max(0.0);
        let mut scores: Vec<f64> = (0..n).map(|s| vectors[s][at] * value.sqrt()).collect();
        let largest = scores
            .iter()
            .copied()
            .max_by(|a, b| a.abs().total_cmp(&b.abs()))
            .unwrap_or(0.0);
        if largest < 0.0 {
            scores.iter_mut().for_each(|v| *v = -*v);
        }
        (value / total, scores)
    };
    let ((first, xs), (second, ys)) = (component(0), component(1));
    Some(Pca {
        groups: complete.len(),
        explained: [first, second],
        scores: xs.into_iter().zip(ys).collect(),
    })
}

/// The eigenvalues of a symmetric matrix, and its eigenvectors as columns
/// (cyclic Jacobi rotations).
fn jacobi(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut vectors = vec![vec![0.0; n]; n];
    for (i, row) in vectors.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let scale: f64 = (0..n)
        .map(|i| a[i][i].abs())
        .sum::<f64>()
        .max(f64::MIN_POSITIVE);
    for _ in 0..JACOBI_SWEEPS {
        let off: f64 = (0..n)
            .flat_map(|p| (p + 1..n).map(move |q| (p, q)))
            .map(|(p, q)| a[p][q] * a[p][q])
            .sum();
        if off.sqrt() <= 1e-13 * scale {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q].abs() <= f64::MIN_POSITIVE {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                let rotate = |row: &mut Vec<f64>| {
                    let (vp, vq) = (row[p], row[q]);
                    row[p] = c * vp - s * vq;
                    row[q] = s * vp + c * vq;
                };
                a.iter_mut().for_each(rotate);
                vectors.iter_mut().for_each(rotate);
                let (top, bottom) = a.split_at_mut(q);
                for (x, y) in top[p].iter_mut().zip(&mut bottom[0]) {
                    let (apk, aqk) = (*x, *y);
                    *x = c * apk - s * aqk;
                    *y = s * apk + c * aqk;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), vectors)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn samples(names: &[&str]) -> Vec<Sample> {
        names
            .iter()
            .map(|name| Sample {
                name: (*name).to_string(),
                category: name.get(..1).map(str::to_string),
            })
            .collect()
    }

    fn group(id: i64, log2: &[Option<f64>]) -> GroupValues {
        GroupValues {
            group_id: id,
            accession: format!("P{id}"),
            log2: log2.to_vec(),
        }
    }

    #[test]
    fn correlation_uses_shared_groups_and_says_when_it_cannot() {
        let s = samples(&["a", "b", "c"]);
        let groups = vec![
            group(1, &[Some(1.0), Some(-1.0), Some(5.0)]),
            group(2, &[Some(2.0), Some(-2.0), None]),
            group(3, &[Some(3.0), Some(-3.0), None]),
            group(4, &[None, None, Some(7.0)]),
        ];
        let qc = sample_qc(&s, &[true; 3], &groups);
        let r = &qc.correlation;
        assert!((r[0][1].unwrap() + 1.0).abs() < 1e-12, "exactly opposed");
        assert!((r[0][0].unwrap() - 1.0).abs() < 1e-12);
        assert_eq!(r[0][1], r[1][0]);
        assert_eq!(r[0][2], None, "one shared group");
        assert_eq!(qc.samples[0].measured, 3);
        assert_eq!(qc.samples[2].measured, 2);
        assert_eq!(qc.samples[2].median_log2, Some(6.0));
        assert_eq!(qc.pca, None, "one group measured in all");
    }

    #[test]
    fn samples_without_data_are_left_out_uncategorized_ones_stay() {
        let mut s = samples(&["a", "b", "dead", "c"]);
        s[3].category = None;
        let groups: Vec<GroupValues> = (0..5)
            .map(|i| {
                let i = f64::from(i);
                group(
                    i as i64,
                    &[Some(i), Some(i * i), None, Some((i * 2.0).sin())],
                )
            })
            .collect();
        let qc = sample_qc(&s, &[true, true, false, true], &groups);
        let names: Vec<&str> = qc.samples.iter().map(|v| v.name.as_str()).collect();
        assert_eq!(names, ["a", "b", "c"]);
        assert_eq!(qc.samples[2].category, None);
        assert_eq!(
            qc.pca.unwrap().groups,
            5,
            "the dead sample doesn't exclude groups"
        );
    }

    /// Two conditions of three, one group in three changing by 4 log2: PC1 is
    /// the condition, and the scores keep all the variance.
    #[test]
    fn pca_separates_conditions_and_keeps_the_variance() {
        let s = samples(&["a1", "a2", "a3", "b1", "b2", "b3"]);
        let groups: Vec<GroupValues> = (0..60)
            .map(|g| {
                let base = 10.0 + f64::from(g % 7);
                let row: Vec<Option<f64>> = (0..6)
                    .map(|c| {
                        let change = if c >= 3 && g % 3 == 0 { 4.0 } else { 0.0 };
                        let noise = 0.05 * (f64::from(g * 31 + c * 17)).sin();
                        Some(base + change + noise)
                    })
                    .collect();
                group(i64::from(g), &row)
            })
            .collect();
        let qc = sample_qc(&s, &[true; 6], &groups);
        let pca = qc.pca.unwrap();
        assert_eq!(pca.groups, 60);
        assert!(pca.explained[0] > 0.95, "{:?}", pca.explained);
        assert!(pca.explained[1] < pca.explained[0]);
        let (a, b) = pca.scores.split_at(3);
        let (a_max, b_min) = (
            a.iter().map(|p| p.0).fold(f64::MIN, f64::max),
            b.iter().map(|p| p.0).fold(f64::MAX, f64::min),
        );
        let (a_min, b_max) = (
            a.iter().map(|p| p.0).fold(f64::MAX, f64::min),
            b.iter().map(|p| p.0).fold(f64::MIN, f64::max),
        );
        assert!(a_max < b_min || b_max < a_min, "{:?}", pca.scores);

        // The first two components' squared scores are their share of the
        // total sum of squares about each group's mean.
        let total: f64 = groups
            .iter()
            .map(|g| {
                let row: Vec<f64> = g.log2.iter().flatten().copied().collect();
                let mean = row.iter().sum::<f64>() / 6.0;
                row.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
            })
            .sum();
        let first: f64 = pca.scores.iter().map(|p| p.0 * p.0).sum();
        let second: f64 = pca.scores.iter().map(|p| p.1 * p.1).sum();
        assert!((first / total - pca.explained[0]).abs() < 1e-9);
        assert!((second / total - pca.explained[1]).abs() < 1e-9);

        // Replicates correlate better than they do with the other condition.
        let r = &qc.correlation;
        assert!(r[0][1].unwrap() > r[0][3].unwrap());
    }

    #[test]
    fn pca_matches_a_decomposition_worked_by_hand() {
        // Centred rows (1, −1, 0) and (1, 1, −2): orthogonal, sums of squares
        // 2 and 6, so 75% and 25%, and PC1's scores are ±(1, 1, −2) exactly.
        let s = samples(&["a", "b", "c"]);
        let groups = vec![
            group(1, &[Some(6.0), Some(4.0), Some(5.0)]),
            group(2, &[Some(6.0), Some(6.0), Some(3.0)]),
            group(3, &[Some(9.0), Some(9.0), Some(9.0)]),
        ];
        let pca = sample_qc(&s, &[true; 3], &groups).pca.unwrap();
        assert!(
            (pca.explained[0] - 0.75).abs() < 1e-12,
            "{:?}",
            pca.explained
        );
        assert!((pca.explained[1] - 0.25).abs() < 1e-12);
        // The largest score, c's, is made positive.
        for (score, want) in pca.scores.iter().zip([-1.0, -1.0, 2.0]) {
            assert!((score.0 - want).abs() < 1e-9, "{:?}", pca.scores);
        }
    }

    #[test]
    fn nothing_to_decompose_is_said_not_invented() {
        let s = samples(&["a", "b", "c"]);
        let flat: Vec<GroupValues> = (0..4)
            .map(|g| group(g, &[Some(2.0), Some(2.0), Some(2.0)]))
            .collect();
        let qc = sample_qc(&s, &[true; 3], &flat);
        assert_eq!(qc.pca, None);
        assert_eq!(qc.correlation[0][1], None);
        let two = sample_qc(&s[..2], &[true; 2], &flat);
        assert_eq!(two.pca, None, "two samples have nothing to plot");
    }
}
