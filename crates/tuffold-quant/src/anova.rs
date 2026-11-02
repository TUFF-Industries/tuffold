//! Comparing every category at once: one-way ANOVA per protein group, on the
//! same per-sample log2 values the two-category tests use (log2 of
//! normalized counts + 1, or log2 median-polish intensities).
//!
//! For spectral counts, each category's spread is floored at its Poisson
//! counting noise ([`crate::stats::counting_floor`]); MS1 intensities have
//! no floor.
//!
//! Classic one-way ANOVA pools the within-category variance. With two or
//! three replicates per category, that pooled estimate is far steadier than
//! each category's own, which is what Welch's ANOVA would need. Only the
//! categories that measured a group take part in its test (missing values
//! stay missing). A group needs two such categories, and at least one degree
//! of freedom within them (more values than categories). Results get
//! Benjamini-Hochberg q-values across the groups tested. A group with no
//! spread at all within categories can't be tested (F is undefined) and is
//! left out, though its difference may be real; with measured intensities
//! that essentially never happens.

use crate::stats::{benjamini_hochberg, incomplete_beta};
use crate::{QuantError, Sample};

/// `P(F(d1, d2) ≥ f)`, the upper tail of the F distribution.
pub fn f_upper_tail(f: f64, d1: f64, d2: f64) -> f64 {
    if f <= 0.0 {
        return 1.0;
    }
    incomplete_beta(d2 / 2.0, d1 / 2.0, d2 / (d2 + d1 * f))
}

/// One-way ANOVA over groups of values: `(F, df between, df within, p)`.
/// Empty groups are ignored. `None` with fewer than two non-empty groups,
/// no within-group degree of freedom, or no within-group variance.
pub fn one_way_anova(groups: &[Vec<f64>]) -> Option<(f64, f64, f64, f64)> {
    one_way_anova_floored(groups, None)
}

/// One-way ANOVA where each group's within sum of squares is at least
/// `floor(group) × (n − 1)`, so no group can claim less spread than the
/// floor allows (counting noise for spectral counts; see
/// [`crate::stats::counting_floor`]).
pub fn one_way_anova_floored(
    groups: &[Vec<f64>],
    floor: Option<fn(&[f64]) -> f64>,
) -> Option<(f64, f64, f64, f64)> {
    let groups: Vec<&Vec<f64>> = groups.iter().filter(|g| !g.is_empty()).collect();
    let k = groups.len();
    let n: usize = groups.iter().map(|g| g.len()).sum();
    if k < 2 || n <= k {
        return None;
    }
    let mean = |values: &[f64]| values.iter().sum::<f64>() / values.len() as f64;
    let grand = groups.iter().flat_map(|g| g.iter()).sum::<f64>() / n as f64;
    let between: f64 = groups
        .iter()
        .map(|g| g.len() as f64 * (mean(g) - grand).powi(2))
        .sum();
    let within: f64 = groups
        .iter()
        .map(|g| {
            let m = mean(g);
            let observed = g.iter().map(|v| (v - m).powi(2)).sum::<f64>();
            let least = floor.map_or(0.0, |f| f(g) * (g.len() - 1) as f64);
            observed.max(least)
        })
        .sum();
    let (d1, d2) = ((k - 1) as f64, (n - k) as f64);
    // No spread within categories, up to rounding: nothing to test against.
    if within <= 1e-24 * grand.powi(2).max(1.0) * n as f64 {
        return None;
    }
    let f = (between / d1) / (within / d2);
    let p = f_upper_tail(f, d1, d2);
    (f.is_finite() && p.is_finite()).then_some((f, d1, d2, p))
}

/// One protein group's per-sample log2 values.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupValues {
    pub group_id: i64,
    pub accession: String,
    /// Parallel to the samples; `None` where not measured.
    pub log2: Vec<Option<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnovaRow {
    pub group_id: i64,
    pub accession: String,
    /// Mean log2 value per category (parallel to
    /// [`AnovaComparison::categories`]), over the samples measuring it.
    pub means: Vec<Option<f64>>,
    pub measured: Vec<usize>,
    /// Highest minus lowest category mean, in log2 units.
    pub range: Option<f64>,
    pub f: Option<f64>,
    pub p: Option<f64>,
    pub q: Option<f64>,
}

/// Every category compared at once, most significant first (untested groups
/// last, by range).
#[derive(Debug, Clone, PartialEq)]
pub struct AnovaComparison {
    pub categories: Vec<String>,
    /// The samples in each category that took part.
    pub members: Vec<Vec<String>>,
    pub rows: Vec<AnovaRow>,
}

/// Compare every category with data. `with_data[s]` says whether sample `s`
/// has any values at all (samples without are left out). Samples without a
/// category are not compared.
pub fn compare_all(
    samples: &[Sample],
    with_data: &[bool],
    groups: &[GroupValues],
    floor: Option<fn(&[f64]) -> f64>,
) -> Result<AnovaComparison, QuantError> {
    let mut categories: Vec<String> = Vec::new();
    let mut members: Vec<Vec<usize>> = Vec::new();
    for (s, sample) in samples.iter().enumerate() {
        let Some(category) = &sample.category else {
            continue;
        };
        if !with_data.get(s).copied().unwrap_or(false) {
            continue;
        }
        match categories.iter().position(|c| c == category) {
            Some(i) => members[i].push(s),
            None => {
                categories.push(category.clone());
                members.push(vec![s]);
            }
        }
    }
    if categories.len() < 2 {
        return Err(QuantError::TooFewCategories(categories));
    }
    let mut rows: Vec<AnovaRow> = groups
        .iter()
        .map(|group| {
            let values: Vec<Vec<f64>> = members
                .iter()
                .map(|m| {
                    m.iter()
                        .filter_map(|&s| group.log2.get(s).copied().flatten())
                        .collect()
                })
                .collect();
            let means: Vec<Option<f64>> = values
                .iter()
                .map(|v| (!v.is_empty()).then(|| v.iter().sum::<f64>() / v.len() as f64))
                .collect();
            let present: Vec<f64> = means.iter().flatten().copied().collect();
            let range = (present.len() >= 2).then(|| {
                present.iter().copied().fold(f64::MIN, f64::max)
                    - present.iter().copied().fold(f64::MAX, f64::min)
            });
            let test = one_way_anova_floored(&values, floor);
            AnovaRow {
                group_id: group.group_id,
                accession: group.accession.clone(),
                measured: values.iter().map(Vec::len).collect(),
                means,
                range,
                f: test.map(|t| t.0),
                p: test.map(|t| t.3),
                q: None,
            }
        })
        .collect();
    let q = benjamini_hochberg(&rows.iter().map(|r| r.p).collect::<Vec<_>>());
    for (row, q) in rows.iter_mut().zip(q) {
        row.q = q;
    }
    // q-values tie in runs (they're monotone), so p breaks ties first.
    rows.sort_by(|x, y| {
        x.q.unwrap_or(2.0)
            .total_cmp(&y.q.unwrap_or(2.0))
            .then(x.p.unwrap_or(2.0).total_cmp(&y.p.unwrap_or(2.0)))
            .then(y.range.unwrap_or(-1.0).total_cmp(&x.range.unwrap_or(-1.0)))
            .then_with(|| x.accession.cmp(&y.accession))
    });
    Ok(AnovaComparison {
        categories,
        members: members
            .iter()
            .map(|m| m.iter().map(|&s| samples[s].name.clone()).collect())
            .collect(),
        rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-9 * expected.abs().max(1e-12),
            "{actual} vs {expected}"
        );
    }

    #[test]
    fn anova_matches_scipy() {
        // scipy.stats.f_oneway([1,2,3], [2,4,6], [5,5.5,7]):
        // F 5.43835616438356, p 0.04493556513850041. By hand: means 2, 4,
        // 5.833; SSB 22.0556 on 2 df, SSW 12.1667 on 6 df.
        let (f, d1, d2, p) = one_way_anova(&[
            vec![1.0, 2.0, 3.0],
            vec![2.0, 4.0, 6.0],
            vec![5.0, 5.5, 7.0],
        ])
        .unwrap();
        close(f, 5.438_356_164_383_56);
        assert_eq!((d1, d2), (2.0, 6.0));
        close(p, 0.044_935_565_138_500_41);
        // Unequal sizes: f_oneway([10.1,10.4], [12,11.6,12.3], [10.2,9.9]).
        let (f, _, _, p) =
            one_way_anova(&[vec![10.1, 10.4], vec![12.0, 11.6, 12.3], vec![10.2, 9.9]]).unwrap();
        close(f, 33.847_241_867_043_735);
        close(p, 0.003_112_780_513_866_589_4);
        // scipy.stats.f.sf.
        close(f_upper_tail(3.5, 2.0, 9.0), 0.075_084_686_279_296_9);
        close(f_upper_tail(250.0, 1.0, 3.0), 0.000_549_973_797_350_436_2);
    }

    #[test]
    fn anova_needs_two_groups_and_a_spread() {
        assert_eq!(one_way_anova(&[vec![1.0, 2.0], vec![]]), None);
        // Three values in three groups: no within-group degree of freedom.
        assert_eq!(one_way_anova(&[vec![1.0], vec![2.0], vec![3.0]]), None);
        assert_eq!(one_way_anova(&[vec![4.0, 4.0], vec![5.0, 5.0]]), None);
    }

    fn sample(name: &str, category: Option<&str>) -> Sample {
        Sample {
            name: name.into(),
            category: category.map(str::to_string),
        }
    }

    #[test]
    fn every_category_is_compared_and_missing_stays_missing() {
        let samples = [
            sample("a1", Some("a")),
            sample("a2", Some("a")),
            sample("b1", Some("b")),
            sample("b2", Some("b")),
            sample("c1", Some("c")),
            sample("c2", Some("c")),
            sample("x", None),
            sample("empty", Some("d")),
        ];
        let with_data = [true, true, true, true, true, true, true, false];
        let group = |id, log2: [Option<f64>; 8]| GroupValues {
            group_id: id,
            accession: format!("P{id}"),
            log2: log2.to_vec(),
        };
        let groups = [
            // Up in c.
            group(
                1,
                [
                    Some(10.0),
                    Some(10.2),
                    Some(10.1),
                    Some(9.9),
                    Some(14.0),
                    Some(14.1),
                    Some(1.0),
                    None,
                ],
            ),
            // Flat.
            group(
                2,
                [
                    Some(8.0),
                    Some(8.3),
                    Some(8.1),
                    Some(8.2),
                    Some(8.2),
                    Some(8.0),
                    None,
                    None,
                ],
            ),
            // Only in a: measured but untestable.
            group(
                3,
                [Some(5.0), Some(5.5), None, None, None, None, None, None],
            ),
        ];
        let anova = compare_all(&samples, &with_data, &groups, None).unwrap();
        // "x" has no category and "d" no data: neither takes part.
        assert_eq!(anova.categories, ["a", "b", "c"]);
        assert_eq!(anova.members[2], ["c1", "c2"]);
        let row = |id| anova.rows.iter().find(|r| r.group_id == id).unwrap();
        assert_eq!(anova.rows[0].group_id, 1);
        assert!(row(1).q.unwrap() < 0.05);
        // Means 10.1, 10.0, 14.05: range 4.05.
        assert!((row(1).range.unwrap() - 4.05).abs() < 1e-9);
        assert!(row(2).p.unwrap() > 0.5);
        assert_eq!((row(3).measured.clone(), row(3).p), (vec![2, 0, 0], None));
        assert_eq!(row(3).means, [Some(5.25), None, None]);
        assert_eq!(anova.rows[2].group_id, 3);
        let one = compare_all(&samples[..2], &with_data[..2], &groups, None);
        assert!(matches!(one, Err(QuantError::TooFewCategories(_))));
    }
}
