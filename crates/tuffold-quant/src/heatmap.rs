//! A heat map of protein groups across samples. Each row's log2 values are
//! centred and scaled (z-scores over the samples that measured it), so rows
//! show their pattern rather than their abundance. Rows are ordered by
//! average-linkage hierarchical clustering, so groups that rise and fall
//! together sit together. Columns are grouped by category.
//!
//! Missing values stay missing: they are drawn as gaps, and the distance
//! between two rows uses only the columns both measured, scaled up to the
//! full width. Two rows sharing fewer than two columns are as far apart as
//! the map allows (one column says nothing about a pattern), and rows
//! measured in fewer than two columns aren't clustered at all: they follow
//! the clustered rows in their given order.

use crate::Sample;
use crate::anova::GroupValues;

#[derive(Debug, Clone, PartialEq)]
pub struct HeatmapRow {
    pub group_id: i64,
    pub accession: String,
    /// z-scores, parallel to [`Heatmap::columns`].
    pub z: Vec<Option<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Heatmap {
    /// Sample names, grouped by category.
    pub columns: Vec<String>,
    pub column_categories: Vec<Option<String>>,
    /// In clustering order.
    pub rows: Vec<HeatmapRow>,
}

/// The samples a comparison of every category uses (with data and a
/// category), grouped by category in first-seen order, as the ANOVA lists
/// them. Uncategorized samples are left out, as the ANOVA leaves them out.
pub fn column_order(samples: &[Sample], with_data: &[bool]) -> Vec<usize> {
    let shown = |s: usize| with_data.get(s).copied().unwrap_or(false);
    let mut categories: Vec<&str> = Vec::new();
    for (s, sample) in samples.iter().enumerate() {
        if let Some(c) = sample.category.as_deref()
            && shown(s)
            && !categories.contains(&c)
        {
            categories.push(c);
        }
    }
    categories
        .iter()
        .flat_map(|c| {
            (0..samples.len())
                .filter(move |&s| shown(s) && samples[s].category.as_deref() == Some(*c))
        })
        .collect()
}

/// The heat map of `groups` over the samples `columns` (indices into
/// `samples`, in display order).
pub fn heatmap(samples: &[Sample], columns: &[usize], groups: &[GroupValues]) -> Heatmap {
    let rows: Vec<HeatmapRow> = groups
        .iter()
        .map(|group| HeatmapRow {
            group_id: group.group_id,
            accession: group.accession.clone(),
            z: z_scores(
                &columns
                    .iter()
                    .map(|&s| group.log2.get(s).copied().flatten())
                    .collect::<Vec<_>>(),
            ),
        })
        .collect();
    let order = cluster_order(&rows, columns.len());
    Heatmap {
        columns: columns.iter().map(|&s| samples[s].name.clone()).collect(),
        column_categories: columns
            .iter()
            .map(|&s| samples[s].category.clone())
            .collect(),
        rows: order.into_iter().map(|i| rows[i].clone()).collect(),
    }
}

/// (v − mean) / sd over the present values; a row with one value or no
/// spread is all 0 where measured.
fn z_scores(values: &[Option<f64>]) -> Vec<Option<f64>> {
    let present: Vec<f64> = values.iter().flatten().copied().collect();
    let n = present.len() as f64;
    let mean = present.iter().sum::<f64>() / n.max(1.0);
    let sd = if present.len() >= 2 {
        (present.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt()
    } else {
        0.0
    };
    values
        .iter()
        .map(|v| v.map(|v| if sd > 1e-12 { (v - mean) / sd } else { 0.0 }))
        .collect()
}

/// Euclidean distance over shared columns, scaled to `width` columns;
/// `None` when the rows share none.
fn distance(a: &[Option<f64>], b: &[Option<f64>], width: usize) -> Option<f64> {
    let shared: Vec<f64> = a
        .iter()
        .zip(b)
        .filter_map(|(x, y)| Some((x.as_ref()? - y.as_ref()?).powi(2)))
        .collect();
    (shared.len() >= 2)
        .then(|| (shared.iter().sum::<f64>() * width as f64 / shared.len() as f64).sqrt())
}

/// Leaf order of average-linkage clustering (UPGMA). Merging two clusters
/// puts the one holding the earlier row first, so the order is
/// deterministic.
fn cluster_order(rows: &[HeatmapRow], width: usize) -> Vec<usize> {
    let (patterned, sparse): (Vec<usize>, Vec<usize>) =
        (0..rows.len()).partition(|&i| rows[i].z.iter().flatten().count() >= 2);
    let picked: Vec<HeatmapRow> = patterned.iter().map(|&i| rows[i].clone()).collect();
    upgma(&picked, width)
        .into_iter()
        .map(|i| patterned[i])
        .chain(sparse)
        .collect()
}

fn upgma(rows: &[HeatmapRow], width: usize) -> Vec<usize> {
    let n = rows.len();
    if n < 3 {
        return (0..n).collect();
    }
    let raw: Vec<Vec<Option<f64>>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| distance(&rows[i].z, &rows[j].z, width))
                .collect()
        })
        .collect();
    // Rows sharing nothing: farther than any measured pair.
    let farthest = raw
        .iter()
        .flatten()
        .flatten()
        .copied()
        .fold(0.0_f64, f64::max);
    let mut d: Vec<Vec<f64>> = raw
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|v| v.unwrap_or(farthest * 2.0 + 1.0))
                .collect()
        })
        .collect();
    let mut clusters: Vec<Option<Vec<usize>>> = (0..n).map(|i| Some(vec![i])).collect();
    for _ in 1..n {
        let mut best: Option<(usize, usize, f64)> = None;
        for i in 0..n {
            if clusters[i].is_none() {
                continue;
            }
            for j in i + 1..n {
                if clusters[j].is_some() && best.is_none_or(|(_, _, b)| d[i][j] < b) {
                    best = Some((i, j, d[i][j]));
                }
            }
        }
        let Some((i, j, _)) = best else { break };
        let (Some(a), Some(b)) = (clusters[i].take(), clusters[j].take()) else {
            break;
        };
        let (na, nb) = (a.len() as f64, b.len() as f64);
        for k in 0..n {
            if clusters[k].is_some() {
                let merged = (na * d[i][k] + nb * d[j][k]) / (na + nb);
                d[i][k] = merged;
                d[k][i] = merged;
            }
        }
        let (first, second) = if a[0] <= b[0] { (a, b) } else { (b, a) };
        clusters[i] = Some(first.into_iter().chain(second).collect());
    }
    clusters.into_iter().flatten().flatten().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str, category: Option<&str>) -> Sample {
        Sample {
            name: name.into(),
            category: category.map(str::to_string),
        }
    }

    fn group(id: i64, log2: &[Option<f64>]) -> GroupValues {
        GroupValues {
            group_id: id,
            accession: format!("P{id}"),
            log2: log2.to_vec(),
        }
    }

    #[test]
    fn columns_group_by_category() {
        let samples = [
            sample("t1", Some("treated")),
            sample("c1", Some("control")),
            sample("x", None),
            sample("t2", Some("treated")),
            sample("c2", Some("control")),
        ];
        let order = column_order(&samples, &[true, true, true, true, false]);
        // Treated first (seen first), c2 has no data, x has no category.
        assert_eq!(order, [0, 3, 1]);
        // A category's first sample without data doesn't set the order.
        let order = column_order(&samples, &[false, true, true, true, true]);
        assert_eq!(order, [1, 4, 3]);
    }

    #[test]
    fn rows_are_z_scored_and_similar_rows_sit_together() {
        let samples: Vec<Sample> = ["a1", "a2", "b1", "b2"]
            .iter()
            .map(|n| sample(n, Some(&n[..1])))
            .collect();
        let up = |id, base: f64| {
            group(
                id,
                &[
                    Some(base),
                    Some(base + 0.1),
                    Some(base + 3.0),
                    Some(base + 3.1),
                ],
            )
        };
        let down = |id, base: f64| {
            group(
                id,
                &[
                    Some(base + 3.0),
                    Some(base + 3.1),
                    Some(base),
                    Some(base + 0.1),
                ],
            )
        };
        // Interleaved up and down rows come out grouped.
        let map = heatmap(
            &samples,
            &[0, 1, 2, 3],
            &[
                up(1, 10.0),
                down(2, 5.0),
                up(3, 20.0),
                down(4, 7.0),
                group(5, &[Some(1.0), None, None, Some(3.0)]),
            ],
        );
        let order: Vec<i64> = map.rows.iter().map(|r| r.group_id).collect();
        let at = |id| order.iter().position(|&g| g == id).unwrap();
        assert_eq!(at(1).abs_diff(at(3)), 1, "{order:?}");
        assert_eq!(at(2).abs_diff(at(4)), 1, "{order:?}");
        // A row measured once has no pattern: it comes last, unclustered.
        let lonely = heatmap(
            &samples,
            &[0, 1, 2, 3],
            &[
                group(9, &[Some(4.0), None, None, None]),
                up(1, 10.0),
                down(2, 5.0),
                up(3, 20.0),
            ],
        );
        assert_eq!(lonely.rows.last().unwrap().group_id, 9);
        assert_eq!(lonely.rows.last().unwrap().z[0], Some(0.0));
        // [1, 3]: mean 2, sd √2, so z = ∓1/√2; missing stays missing.
        let five = map.rows.iter().find(|r| r.group_id == 5).unwrap();
        let half = std::f64::consts::FRAC_1_SQRT_2;
        assert!(
            (five.z[0].unwrap() + half).abs() < 1e-12 && (five.z[3].unwrap() - half).abs() < 1e-12
        );
        assert_eq!((five.z[1], five.z[2]), (None, None));
        assert_eq!(map.columns, ["a1", "a2", "b1", "b2"]);
    }
}
