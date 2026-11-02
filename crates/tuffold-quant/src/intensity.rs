//! Quantitation from MS1 precursor intensity.
//!
//! Input: per protein group, the MS1 peak area of each precursor (a unique
//! or razor peptide at one charge) in each biosample (`None` where it has
//! no feature there). "Peptide" below means such a precursor. Everything is compared on a log2 scale.
//!
//! - **Normalization** happens at the peptide level, before proteins are
//!   built: each pair of biosamples is compared on the peptides both
//!   measured (the half-sample mode of their log2 ratios), the pairwise
//!   ratios are fitted by least squares (the pairwise-ratio fit of MaxLFQ, Cox
//!   et al. 2014), and each biosample is shifted to the common level. That corrects for loading and
//!   spray differences under the assumption that the largest group of
//!   peptides doesn't change. A pair of biosamples sharing fewer than
//!   [`MIN_SHARED_FOR_NORMALIZATION`] peptides doesn't count, since a mode of
//!   a handful of peptides could be one regulated protein's; a biosample with
//!   no such pair is left unnormalized and listed as such. A mode rather than a median, because the
//!   peptides two samples share lean toward the ones up in either (those are
//!   identified more often), which drags a median along. Without any shared
//!   peptide, nothing is shifted, and the result says so.
//! - **Protein intensities** are Tukey median-polish profiles of the
//!   normalized peptides (`polish.rs`), so a sample that identified extra
//!   peptides of a protein doesn't look more abundant for it.
//! - **Missing values stay missing.** There is no imputation: a group absent
//!   from a biosample has no value there, its means use only the samples that
//!   have one, and Welch's t-test runs only with two or more values in each
//!   category. Biosamples with no intensities at all (no MS1 data, or nothing
//!   passing) are left out of comparisons.
//! - **Comparisons**: the log2 fold change is mean B − mean A of the log2
//!   protein intensities, and Welch's t-test is applied to the same values,
//!   with Benjamini-Hochberg q-values over the groups tested.

use crate::measures::top3;
use crate::pairwise::relative_levels;
use crate::polish::median_polish;

/// Peptides two biosamples must share for their ratio to normalize them.
pub const MIN_SHARED_FOR_NORMALIZATION: usize = 10;
use crate::stats::{benjamini_hochberg, welch};
use crate::{QuantError, Sample};

#[derive(Debug, Clone, PartialEq)]
pub struct ProteinPeptides {
    pub group_id: i64,
    pub accession: String,
    /// Per peptide, its summed peak area in each sample (parallel to the
    /// samples).
    pub peptides: Vec<Vec<Option<f64>>>,
    /// Parallel to `peptides`: whether a value includes a match-between-runs
    /// transfer. Such values are quantified but never normalized on.
    pub transferred: Vec<Vec<bool>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IntensityRow {
    pub group_id: i64,
    pub accession: String,
    /// Median-polish intensity per sample (normalized, linear scale).
    pub intensities: Vec<Option<f64>>,
    /// The same, log2.
    pub log2: Vec<Option<f64>>,
    /// Peptides measured in each sample.
    pub peptides: Vec<usize>,
    /// Top 3 per sample: the mean normalized area of its three most intense
    /// precursors ([`crate::measures::top3`]), linear scale.
    pub top3: Vec<Option<f64>>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct IntensityQuantities {
    pub samples: Vec<Sample>,
    /// log2 shift applied to each sample's peptides (0 when not normalized).
    pub shifts: Vec<f64>,
    /// Whether each sample was normalized (it shares enough peptides with
    /// another).
    pub normalized: Vec<bool>,
    /// Peptides measured in two or more samples, which normalization rests
    /// on; 0 means nothing was normalized.
    pub shared_peptides: usize,
    pub rows: Vec<IntensityRow>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IntensityComparisonRow {
    pub group_id: i64,
    pub accession: String,
    /// Mean log2 intensity, over the samples measuring it.
    pub mean_a: Option<f64>,
    pub mean_b: Option<f64>,
    pub measured_a: usize,
    pub measured_b: usize,
    /// mean B − mean A; `None` unless both categories measured it.
    pub log2_fold: Option<f64>,
    /// `None` unless both categories have two or more values (that vary).
    pub p: Option<f64>,
    pub q: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IntensityComparison {
    pub a: String,
    pub b: String,
    pub samples_a: Vec<String>,
    pub samples_b: Vec<String>,
    pub rows: Vec<IntensityComparisonRow>,
}

/// Normalize the peptides and build every group's protein intensities.
pub fn quantify_intensities(
    samples: Vec<Sample>,
    proteins: Vec<ProteinPeptides>,
) -> Result<IntensityQuantities, QuantError> {
    let n = samples.len();
    if let Some(bad) = proteins.iter().find(|p| {
        p.transferred.len() != p.peptides.len()
            || p.peptides.iter().any(|peptide| peptide.len() != n)
            || p.transferred.iter().any(|flags| flags.len() != n)
    }) {
        return Err(QuantError::Misaligned {
            accession: bad.accession.clone(),
            found: bad
                .peptides
                .iter()
                .map(Vec::len)
                .find(|&l| l != n)
                .unwrap_or(0),
            expected: n,
        });
    }
    // Non-finite areas (corrupt input) are treated as unmeasured.
    let logged = |value: Option<f64>| value.filter(|v| v.is_finite() && *v > 0.0).map(f64::log2);
    // Normalization rests on measured values only: a transfer is biased low
    // (the precursor was not identified there, so it is usually weaker than
    // where it was) and unevenly spread over the samples, which would tilt
    // the shifts.
    let for_normalization: Vec<Vec<Option<f64>>> = proteins
        .iter()
        .flat_map(|p| p.peptides.iter().zip(&p.transferred))
        .map(|(peptide, transferred)| {
            peptide
                .iter()
                .zip(transferred)
                .map(|(&v, &transfer)| logged(v).filter(|_| !transfer))
                .collect()
        })
        .collect();
    let shared_peptides = for_normalization
        .iter()
        .filter(|p| p.iter().flatten().count() >= 2)
        .count();
    let levels = relative_levels(&for_normalization, n, MIN_SHARED_FOR_NORMALIZATION);
    let normalized: Vec<bool> = levels.iter().map(Option::is_some).collect();
    let shifts: Vec<f64> = levels
        .into_iter()
        .map(|level| -level.unwrap_or(0.0))
        .collect();
    let rows = proteins
        .into_iter()
        .map(|protein| {
            let normalized: Vec<Vec<Option<f64>>> = protein
                .peptides
                .iter()
                .map(|peptide| {
                    peptide
                        .iter()
                        .zip(&shifts)
                        .map(|(&v, shift)| logged(v).map(|l| l + shift))
                        .collect()
                })
                .collect();
            let log2 = median_polish(&normalized, n);
            IntensityRow {
                group_id: protein.group_id,
                accession: protein.accession,
                intensities: log2.iter().map(|v| v.map(f64::exp2)).collect(),
                peptides: (0..n)
                    .map(|s| normalized.iter().filter(|p| p[s].is_some()).count())
                    .collect(),
                top3: (0..n)
                    .map(|s| top3(normalized.iter().filter_map(|p| p[s]).map(f64::exp2)))
                    .collect(),
                log2,
            }
        })
        .collect();
    Ok(IntensityQuantities {
        samples,
        shifts,
        normalized,
        shared_peptides,
        rows,
    })
}

impl IntensityQuantities {
    /// Per-sample log2 median-polish intensity of every group, and which
    /// samples have any (the values ANOVA and the heat map use).
    pub fn log2_values(&self) -> (Vec<bool>, Vec<crate::anova::GroupValues>) {
        let with_data: Vec<bool> = (0..self.samples.len())
            .map(|s| self.rows.iter().any(|row| row.log2[s].is_some()))
            .collect();
        let groups = self
            .rows
            .iter()
            .map(|row| crate::anova::GroupValues {
                group_id: row.group_id,
                accession: row.accession.clone(),
                log2: row.log2.clone(),
            })
            .collect();
        (with_data, groups)
    }

    /// Sample QC on the log2 median-polish intensities.
    pub fn sample_qc(&self) -> crate::qc::SampleQc {
        let (with_data, groups) = self.log2_values();
        crate::qc::sample_qc(&self.samples, &with_data, &groups)
    }

    /// Every category compared at once (one-way ANOVA on log2 median-polish
    /// intensities); samples without any intensity are left out.
    pub fn anova(&self) -> Result<crate::anova::AnovaComparison, QuantError> {
        let (with_data, groups) = self.log2_values();
        crate::anova::compare_all(&self.samples, &with_data, &groups, None)
    }

    /// Distinct categories, in first-seen order.
    pub fn categories(&self) -> Vec<String> {
        let mut seen: Vec<String> = Vec::new();
        for category in self.samples.iter().filter_map(|s| s.category.clone()) {
            if !seen.contains(&category) {
                seen.push(category);
            }
        }
        seen
    }

    fn members(&self, category: &str) -> Result<Vec<usize>, QuantError> {
        let members: Vec<usize> = (0..self.samples.len())
            .filter(|&s| self.samples[s].category.as_deref() == Some(category))
            .collect();
        if members.is_empty() {
            return Err(QuantError::EmptyCategory {
                category: category.to_string(),
                known: self.categories(),
            });
        }
        let measured: Vec<usize> = members
            .into_iter()
            .filter(|&s| self.rows.iter().any(|row| row.log2[s].is_some()))
            .collect();
        if measured.is_empty() {
            return Err(QuantError::NoIntensities(category.to_string()));
        }
        Ok(measured)
    }

    /// Compare category `b` against category `a`, most significant first
    /// (untested groups last, by fold change).
    pub fn compare(&self, a: &str, b: &str) -> Result<IntensityComparison, QuantError> {
        if a == b {
            return Err(QuantError::SameCategory(a.to_string()));
        }
        let (members_a, members_b) = (self.members(a)?, self.members(b)?);
        let mut rows: Vec<IntensityComparisonRow> = self
            .rows
            .iter()
            .map(|row| {
                let values = |members: &[usize]| -> Vec<f64> {
                    members.iter().filter_map(|&s| row.log2[s]).collect()
                };
                let (va, vb) = (values(&members_a), values(&members_b));
                let mean =
                    |v: &[f64]| (!v.is_empty()).then(|| v.iter().sum::<f64>() / v.len() as f64);
                let (mean_a, mean_b) = (mean(&va), mean(&vb));
                IntensityComparisonRow {
                    group_id: row.group_id,
                    accession: row.accession.clone(),
                    mean_a,
                    mean_b,
                    measured_a: va.len(),
                    measured_b: vb.len(),
                    log2_fold: mean_a.zip(mean_b).map(|(ma, mb)| mb - ma),
                    p: welch(&va, &vb).map(|(_, _, p)| p),
                    q: None,
                }
            })
            .collect();
        let q = benjamini_hochberg(&rows.iter().map(|r| r.p).collect::<Vec<_>>());
        for (row, q) in rows.iter_mut().zip(q) {
            row.q = q;
        }
        let fold = |r: &IntensityComparisonRow| r.log2_fold.map_or(-1.0, f64::abs);
        rows.sort_by(|x, y| {
            x.q.unwrap_or(2.0)
                .total_cmp(&y.q.unwrap_or(2.0))
                .then(fold(y).total_cmp(&fold(x)))
                .then_with(|| x.accession.cmp(&y.accession))
        });
        let names = |members: &[usize]| {
            members
                .iter()
                .map(|&s| self.samples[s].name.clone())
                .collect()
        };
        Ok(IntensityComparison {
            a: a.to_string(),
            b: b.to_string(),
            samples_a: names(&members_a),
            samples_b: names(&members_b),
            rows,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str, category: &str) -> Sample {
        Sample {
            name: name.into(),
            category: Some(category.into()),
        }
    }

    /// A one-peptide protein: its intensity is the peptide's.
    fn protein(id: i64, intensities: &[Option<f64>]) -> ProteinPeptides {
        ProteinPeptides {
            group_id: id,
            accession: format!("P{id}"),
            peptides: vec![intensities.to_vec()],
            transferred: vec![vec![false; intensities.len()]],
        }
    }

    fn four() -> Vec<Sample> {
        vec![
            sample("c1", "control"),
            sample("c2", "control"),
            sample("t1", "treated"),
            sample("t2", "treated"),
        ]
    }

    #[test]
    fn peptides_are_shifted_to_a_common_median() {
        // c2 was loaded twice as much: each of 12 shared peptides reads 2×
        // there (log2 ratio 1 for every pair), so c1 shifts +0.5 and c2 −0.5.
        let mut proteins: Vec<ProteinPeptides> = (0..12)
            .map(|i| {
                let v = f64::from(1 << (10 + i));
                protein(i64::from(i), &[Some(v), Some(2.0 * v)])
            })
            .collect();
        proteins.push(protein(99, &[None, Some(64.0)]));
        let q = quantify_intensities(
            vec![sample("c1", "control"), sample("c2", "control")],
            proteins,
        )
        .unwrap();
        assert_eq!(q.shifts, [0.5, -0.5]);
        assert_eq!(
            (q.normalized.clone(), q.shared_peptides),
            (vec![true, true], 12)
        );
        assert_eq!(q.rows[0].log2, [Some(10.5), Some(10.5)]);
        // Missing stays missing; P99 is shifted like the rest of c2.
        assert_eq!(q.rows[12].log2, [None, Some(5.5)]);
        assert_eq!(q.rows[12].peptides, [0, 1]);
    }

    #[test]
    fn transferred_values_are_quantified_but_never_normalized_on() {
        // 12 peptides measured in both samples at the same level, and 20
        // more measured in c1 and transferred into c2 at 4× (log2 +2). If
        // the transfers counted, 20 of the 32 shared peptides would say c2
        // is 4× and the shifts would be ±1; they don't: the 12 measured
        // ones say 0, so nothing is shifted, 12 peptides were shared, and
        // the transferred values keep their +2 (c1 = 2^10 = 1024, c2 =
        // 2^12 = 4096).
        let mut proteins: Vec<ProteinPeptides> = (0..12)
            .map(|i| protein(i, &[Some(1024.0), Some(1024.0)]))
            .collect();
        proteins.extend((12..32).map(|i| ProteinPeptides {
            transferred: vec![vec![false, true]],
            ..protein(i, &[Some(1024.0), Some(4096.0)])
        }));
        let q = quantify_intensities(
            vec![sample("c1", "control"), sample("c2", "control")],
            proteins,
        )
        .unwrap();
        assert_eq!(q.shifts, [0.0, 0.0]);
        assert_eq!(q.shared_peptides, 12);
        assert_eq!(q.rows[0].log2, [Some(10.0), Some(10.0)]);
        assert_eq!(q.rows[20].log2, [Some(10.0), Some(12.0)]);
        // A flag list that doesn't match the peptides is refused.
        let mut bad = protein(1, &[Some(1.0), Some(2.0)]);
        bad.transferred = Vec::new();
        assert!(matches!(
            quantify_intensities(vec![sample("a", "x"), sample("b", "x")], vec![bad]),
            Err(QuantError::Misaligned { .. })
        ));
    }

    #[test]
    fn too_few_shared_peptides_normalize_nothing() {
        // Two shared peptides, both doubled: that could be one regulated
        // protein, so nothing is shifted, and the samples say so.
        let q = quantify_intensities(
            vec![sample("c1", "control"), sample("t1", "treated")],
            vec![
                protein(1, &[Some(8.0), Some(16.0)]),
                protein(2, &[Some(32.0), Some(64.0)]),
                protein(3, &[None, Some(16.0)]),
            ],
        )
        .unwrap();
        assert_eq!(q.shifts, [0.0, 0.0]);
        assert_eq!(q.normalized, [false, false]);
        assert_eq!(q.rows[2].log2, [None, Some(4.0)]);
    }

    #[test]
    fn replicates_get_welch_and_a_fold_change() {
        // P2, P4, and P5 are unchanged, so the medians (of P1, P2, P4, P5)
        // don't move and nothing is shifted; P1 is 8× in treated.
        let unchanged = |id: i64, v: f64| protein(id, &[Some(v); 4]);
        let q = quantify_intensities(
            four(),
            vec![
                protein(1, &[Some(100.0), Some(120.0), Some(800.0), Some(960.0)]),
                unchanged(2, 1e6),
                protein(3, &[Some(50.0), None, Some(60.0), None]),
                unchanged(4, 2e6),
                unchanged(5, 3e6),
            ],
        )
        .unwrap();
        assert_eq!(q.shifts, [0.0; 4]);
        let comparison = q.compare("control", "treated").unwrap();
        let p1 = comparison.rows.iter().find(|r| r.group_id == 1).unwrap();
        assert!((p1.log2_fold.unwrap() - 3.0).abs() < 1e-9);
        assert!(p1.p.unwrap() < 0.05 && p1.q.is_some());
        assert_eq!(comparison.rows[0].group_id, 1);
        // One value per category: a fold change, but no test.
        let p3 = comparison.rows.iter().find(|r| r.group_id == 3).unwrap();
        assert_eq!((p3.measured_a, p3.measured_b, p3.p), (1, 1, None));
        assert!(p3.log2_fold.is_some() && p3.q.is_none());
    }

    #[test]
    fn categories_without_intensities_are_refused() {
        let q = quantify_intensities(
            vec![sample("c1", "control"), sample("t1", "treated")],
            vec![protein(1, &[Some(5.0), None])],
        )
        .unwrap();
        assert_eq!(
            q.compare("control", "treated").unwrap_err(),
            QuantError::NoIntensities("treated".into())
        );
        let misaligned = quantify_intensities(four(), vec![protein(1, &[Some(1.0)])]);
        assert!(matches!(misaligned, Err(QuantError::Misaligned { .. })));
    }
}
