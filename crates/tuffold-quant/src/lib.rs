//! Tuffold quantitation from spectral counts (this module), and from MS1
//! precursor intensity ([`ms1`] traces peaks, [`polish`] builds proteins from
//! peptides, [`intensity`] normalizes and compares them).
//!
//! Input: per protein group, the passing spectra of its unique and razor
//! peptides in each biosample (what the store's protein list reports).
//!
//! - **Normalized counts** scale each biosample to the mean total, so a
//!   sample that yielded more spectra overall does not look more abundant.
//!   Being shares of a total, they are compositional: when many proteins
//!   rise, unchanged ones look slightly down.
//! - **NSAF** (normalized spectral abundance factor, Zybailov et al. 2006):
//!   `(count / length) / Σ (count / length)` per biosample, for proteins with
//!   a known sequence length; `None` otherwise (never guessed).
//! - **Comparisons** between two categories of biosamples (e.g. control vs
//!   treated): log2 fold change of mean normalized counts (B over A, with a
//!   pseudocount of 0.5 so absent proteins stay finite); Fisher's exact test
//!   on summed counts (protein vs every other spectrum, A vs B — works with
//!   one sample per category); Welch's t-test on log2(normalized + 1) per
//!   sample when both categories have at least two samples, each category's
//!   variance floored at its Poisson counting noise
//!   ([`stats::counting_floor`]): counts of 1, 1 against 0, 0 have no
//!   spread, but aren't evidence of a change. Both tests get
//!   Benjamini-Hochberg q-values. With replicates in both categories the
//!   **primary** test is Welch's (it sees replicate variation); pooled
//!   counts treat every spectrum as independent and overstate significance
//!   there, so Fisher is primary only without replicates.
//! - Biosamples with no counted spectra at all (a failed run, or thresholds
//!   nothing passes) are left out of means and t-tests rather than counted
//!   as zeros.

pub mod align;
pub mod anova;
pub mod heatmap;
pub mod intensity;
pub mod measures;
pub mod ms1;
pub mod pairwise;
mod polish;
pub mod purity;
pub mod qc;
pub mod reporter;
pub mod stats;

use stats::{benjamini_hochberg, counting_floor, fisher_exact, mann_whitney, welch_floored};

/// Added to mean normalized counts before the fold change.
pub const PSEUDOCOUNT: f64 = 0.5;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sample {
    pub name: String,
    /// The group it belongs to for comparisons (e.g. "control").
    pub category: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ProteinCounts {
    pub group_id: i64,
    pub accession: String,
    /// Sequence length of the group's leader, when known.
    pub length: Option<usize>,
    /// Parallel to the samples.
    pub counts: Vec<usize>,
    /// Spectra of unique peptides only, per sample (empty when unknown).
    pub exclusive_counts: Vec<usize>,
    /// Distinct unique peptide sequences, per sample (empty when unknown).
    pub exclusive_peptides: Vec<usize>,
    /// Distinct unique and razor peptide sequences, per sample (empty when
    /// unknown), emPAI's observed peptides.
    pub peptides: Vec<usize>,
    /// The leader's observable tryptic peptides, when its sequence is known.
    pub observable: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct QuantRow {
    pub group_id: i64,
    pub accession: String,
    pub counts: Vec<usize>,
    pub normalized: Vec<f64>,
    /// `None` without a sequence length.
    pub nsaf: Option<Vec<f64>>,
    /// Spectra of unique peptides only, per sample.
    pub exclusive_counts: Vec<usize>,
    /// Distinct unique peptide sequences, per sample.
    pub exclusive_peptides: Vec<usize>,
    /// emPAI per sample ([`measures::empai`]); `None` without the
    /// leader's sequence or peptide counts.
    pub empai: Option<Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Quantities {
    pub samples: Vec<Sample>,
    /// Total counted spectra per sample (the normalization basis).
    pub totals: Vec<usize>,
    pub rows: Vec<QuantRow>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ComparisonRow {
    pub group_id: i64,
    pub accession: String,
    pub mean_a: f64,
    pub mean_b: f64,
    /// log2((mean B + 0.5) / (mean A + 0.5)) of normalized counts.
    pub log2_fold: f64,
    pub fisher_p: f64,
    pub fisher_q: f64,
    /// `None` unless both categories have two or more samples (and vary).
    pub welch_p: Option<f64>,
    pub welch_q: Option<f64>,
    /// Mann–Whitney U on the samples' normalized counts; `None` unless both
    /// categories have two or more samples.
    pub mann_whitney_p: Option<f64>,
    pub mann_whitney_q: Option<f64>,
}

/// Which test decides significance in a comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimaryTest {
    /// Both categories have replicates.
    Welch,
    /// A category has a single sample: pooled counts are all there is.
    Fisher,
    /// Rank-based, chosen by the user (never the default): no normality
    /// assumption, but with few replicates it cannot reach significance
    /// ([`stats::mann_whitney_floor`]).
    MannWhitney,
}

impl PrimaryTest {
    pub fn describe(self) -> &'static str {
        match self {
            Self::Welch => "Welch's t-test across replicates",
            Self::Fisher => "Fisher's exact test on pooled counts (no replicates)",
            Self::MannWhitney => "Mann–Whitney U test across replicates (rank-based)",
        }
    }
}

impl ComparisonRow {
    /// The q-value of the primary test (1 when it could not run).
    pub fn q(&self, primary: PrimaryTest) -> f64 {
        match primary {
            PrimaryTest::Welch => self.welch_q.unwrap_or(1.0),
            PrimaryTest::Fisher => self.fisher_q,
            PrimaryTest::MannWhitney => self.mann_whitney_q.unwrap_or(1.0),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Comparison {
    pub a: String,
    pub b: String,
    pub primary: PrimaryTest,
    pub samples_a: Vec<String>,
    pub samples_b: Vec<String>,
    pub rows: Vec<ComparisonRow>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QuantError {
    #[error("no biosample is in category '{category}' (categories: {})", known.join(", "))]
    EmptyCategory {
        category: String,
        known: Vec<String>,
    },
    #[error("no biosample in category '{0}' has any counted spectra")]
    NoSpectra(String),
    #[error(
        "no biosample in category '{0}' has any MS1 intensities; add its spectra files \
         (mzML with MS1 scans)"
    )]
    NoIntensities(String),
    #[error(
        "comparing every category needs two or more with data; found {}",
        if .0.is_empty() { "none".to_string() } else { .0.join(", ") }
    )]
    TooFewCategories(Vec<String>),
    #[error("compare two different categories, not '{0}' with itself")]
    SameCategory(String),
    #[error("protein {accession} has {found} counts for {expected} biosamples")]
    Misaligned {
        accession: String,
        found: usize,
        expected: usize,
    },
    #[error(
        "plex {0} has no reference channel; mark the channel holding the pooled reference \
         sample to combine plexes"
    )]
    NoReference(String),
    #[error("plex {0} holds only reference channels; nothing in it to compare")]
    OnlyReference(String),
}

/// Normalized counts and NSAF for every protein.
pub fn quantify(
    samples: Vec<Sample>,
    proteins: Vec<ProteinCounts>,
) -> Result<Quantities, QuantError> {
    if let Some(bad) = proteins.iter().find(|p| p.counts.len() != samples.len()) {
        return Err(QuantError::Misaligned {
            accession: bad.accession.clone(),
            found: bad.counts.len(),
            expected: samples.len(),
        });
    }
    let totals: Vec<usize> = (0..samples.len())
        .map(|s| proteins.iter().map(|p| p.counts[s]).sum())
        .collect();
    let nonzero: Vec<usize> = totals.iter().copied().filter(|&t| t > 0).collect();
    let mean_total = if nonzero.is_empty() {
        0.0
    } else {
        nonzero.iter().sum::<usize>() as f64 / nonzero.len() as f64
    };
    // Σ count / length per sample, over proteins with a length.
    let saf_totals: Vec<f64> = (0..samples.len())
        .map(|s| {
            proteins
                .iter()
                .filter_map(|p| {
                    p.length
                        .filter(|&l| l > 0)
                        .map(|l| p.counts[s] as f64 / l as f64)
                })
                .sum()
        })
        .collect();
    let rows = proteins
        .into_iter()
        .map(|protein| {
            let normalized = protein
                .counts
                .iter()
                .zip(&totals)
                .map(|(&count, &total)| {
                    if total == 0 {
                        0.0
                    } else {
                        count as f64 * mean_total / total as f64
                    }
                })
                .collect();
            let nsaf = protein.length.filter(|&l| l > 0).map(|length| {
                protein
                    .counts
                    .iter()
                    .zip(&saf_totals)
                    .map(|(&count, &saf)| {
                        if saf > 0.0 {
                            count as f64 / length as f64 / saf
                        } else {
                            0.0
                        }
                    })
                    .collect()
            });
            let empai = protein
                .observable
                .filter(|_| protein.peptides.len() == protein.counts.len())
                .and_then(|observable| {
                    protein
                        .peptides
                        .iter()
                        .map(|&observed| measures::empai(observed, observable))
                        .collect()
                });
            QuantRow {
                group_id: protein.group_id,
                accession: protein.accession,
                counts: protein.counts,
                normalized,
                nsaf,
                exclusive_counts: protein.exclusive_counts,
                exclusive_peptides: protein.exclusive_peptides,
                empai,
            }
        })
        .collect();
    Ok(Quantities {
        samples,
        totals,
        rows,
    })
}

impl Quantities {
    /// Per-sample log2(normalized + 1) of every group, and which samples
    /// have counted spectra (the values ANOVA and the heat map use).
    pub fn log2_values(&self) -> (Vec<bool>, Vec<anova::GroupValues>) {
        let with_data: Vec<bool> = self.totals.iter().map(|&t| t > 0).collect();
        let groups = self
            .rows
            .iter()
            .map(|row| anova::GroupValues {
                group_id: row.group_id,
                accession: row.accession.clone(),
                log2: row
                    .normalized
                    .iter()
                    .zip(&with_data)
                    .map(|(&v, &data)| data.then(|| (v + 1.0).log2()))
                    .collect(),
            })
            .collect();
        (with_data, groups)
    }

    /// Sample QC on log2(normalized + 1).
    pub fn sample_qc(&self) -> qc::SampleQc {
        let (with_data, groups) = self.log2_values();
        qc::sample_qc(&self.samples, &with_data, &groups)
    }

    /// Every category compared at once (one-way ANOVA on log2(normalized
    /// + 1)); samples without counted spectra are left out.
    pub fn anova(&self) -> Result<anova::AnovaComparison, QuantError> {
        let (with_data, groups) = self.log2_values();
        anova::compare_all(&self.samples, &with_data, &groups, Some(counting_floor))
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
        // A sample without any counted spectra carries no information.
        let counted: Vec<usize> = members
            .into_iter()
            .filter(|&s| self.totals[s] > 0)
            .collect();
        if counted.is_empty() {
            return Err(QuantError::NoSpectra(category.to_string()));
        }
        Ok(counted)
    }

    /// Compare category `b` against category `a`, most significant first.
    pub fn compare(&self, a: &str, b: &str) -> Result<Comparison, QuantError> {
        if a == b {
            return Err(QuantError::SameCategory(a.to_string()));
        }
        let (members_a, members_b) = (self.members(a)?, self.members(b)?);
        let total = |members: &[usize]| members.iter().map(|&s| self.totals[s] as u64).sum::<u64>();
        let (total_a, total_b) = (total(&members_a), total(&members_b));
        let mut rows: Vec<ComparisonRow> = self
            .rows
            .iter()
            .map(|row| {
                let sum =
                    |members: &[usize]| members.iter().map(|&s| row.counts[s] as u64).sum::<u64>();
                let mean = |members: &[usize]| {
                    members.iter().map(|&s| row.normalized[s]).sum::<f64>() / members.len() as f64
                };
                let logged = |members: &[usize]| -> Vec<f64> {
                    members
                        .iter()
                        .map(|&s| (row.normalized[s] + 1.0).log2())
                        .collect()
                };
                let (count_a, count_b) = (sum(&members_a), sum(&members_b));
                let (mean_a, mean_b) = (mean(&members_a), mean(&members_b));
                ComparisonRow {
                    group_id: row.group_id,
                    accession: row.accession.clone(),
                    mean_a,
                    mean_b,
                    log2_fold: ((mean_b + PSEUDOCOUNT) / (mean_a + PSEUDOCOUNT)).log2(),
                    fisher_p: fisher_exact(count_a, count_b, total_a - count_a, total_b - count_b),
                    fisher_q: 1.0,
                    welch_p: {
                        let (a, b) = (logged(&members_a), logged(&members_b));
                        welch_floored(&a, &b, counting_floor(&a), counting_floor(&b))
                            .map(|(_, _, p)| p)
                    },
                    welch_q: None,
                    mann_whitney_p: (members_a.len() >= 2 && members_b.len() >= 2)
                        .then(|| {
                            let values = |members: &[usize]| -> Vec<f64> {
                                members.iter().map(|&s| row.normalized[s]).collect()
                            };
                            mann_whitney(&values(&members_a), &values(&members_b))
                        })
                        .flatten(),
                    mann_whitney_q: None,
                }
            })
            .collect();
        let primary = if members_a.len() >= 2 && members_b.len() >= 2 {
            PrimaryTest::Welch
        } else {
            PrimaryTest::Fisher
        };
        let fisher_q =
            benjamini_hochberg(&rows.iter().map(|r| Some(r.fisher_p)).collect::<Vec<_>>());
        let welch_q = benjamini_hochberg(&rows.iter().map(|r| r.welch_p).collect::<Vec<_>>());
        let mann_whitney_q =
            benjamini_hochberg(&rows.iter().map(|r| r.mann_whitney_p).collect::<Vec<_>>());
        for (((row, fisher), welch), rank) in rows
            .iter_mut()
            .zip(fisher_q)
            .zip(welch_q)
            .zip(mann_whitney_q)
        {
            row.fisher_q = fisher.unwrap_or(1.0);
            row.welch_q = welch;
            row.mann_whitney_q = rank;
        }
        rows.sort_by(|x, y| {
            x.q(primary)
                .total_cmp(&y.q(primary))
                .then(x.fisher_q.total_cmp(&y.fisher_q))
                .then(y.log2_fold.abs().total_cmp(&x.log2_fold.abs()))
                .then_with(|| x.accession.cmp(&y.accession))
        });
        let names = |members: &[usize]| {
            members
                .iter()
                .map(|&s| self.samples[s].name.clone())
                .collect()
        };
        Ok(Comparison {
            a: a.to_string(),
            b: b.to_string(),
            primary,
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

    fn protein(id: i64, length: Option<usize>, counts: &[usize]) -> ProteinCounts {
        ProteinCounts {
            group_id: id,
            accession: format!("P{id}"),
            length,
            counts: counts.to_vec(),
            ..ProteinCounts::default()
        }
    }

    /// Two samples (c1 control, t1 treated). Totals: c1 = 10 + 30 = 40,
    /// t1 = 30 + 50 = 80; mean total 60.
    fn two_samples() -> Quantities {
        quantify(
            vec![sample("c1", "control"), sample("t1", "treated")],
            vec![
                protein(1, Some(100), &[10, 30]),
                protein(2, None, &[30, 50]),
            ],
        )
        .unwrap()
    }

    #[test]
    fn counts_scale_to_the_mean_total() {
        let q = two_samples();
        assert_eq!(q.totals, [40, 80]);
        // P1: 10 · 60/40 = 15; 30 · 60/80 = 22.5.
        assert_eq!(q.rows[0].normalized, [15.0, 22.5]);
        // P2: 30 · 1.5 = 45; 50 · 0.75 = 37.5.
        assert_eq!(q.rows[1].normalized, [45.0, 37.5]);
    }

    #[test]
    fn nsaf_needs_a_length_and_sums_over_known_ones() {
        let q = two_samples();
        // P1 is the only protein with a length, so its NSAF is 1 in both.
        assert_eq!(q.rows[0].nsaf, Some(vec![1.0, 1.0]));
        assert_eq!(q.rows[1].nsaf, None);
        let q = quantify(
            vec![sample("s", "x")],
            vec![protein(1, Some(100), &[10]), protein(2, Some(50), &[10])],
        )
        .unwrap();
        // SAF 0.1 and 0.2 of 0.3: NSAF 1/3 and 2/3.
        let nsaf: Vec<f64> = q.rows.iter().map(|r| r.nsaf.as_ref().unwrap()[0]).collect();
        assert!((nsaf[0] - 1.0 / 3.0).abs() < 1e-12 && (nsaf[1] - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn one_sample_per_category_gets_fisher_but_not_welch() {
        let comparison = two_samples().compare("control", "treated").unwrap();
        let p1 = comparison.rows.iter().find(|r| r.group_id == 1).unwrap();
        // log2((22.5 + 0.5) / (15 + 0.5)) = log2(23 / 15.5).
        assert!((p1.log2_fold - (23.0_f64 / 15.5).log2()).abs() < 1e-12);
        // Table [[10, 30], [30, 50]]: protein vs the rest, control vs treated.
        assert_eq!(p1.fisher_p, fisher_exact(10, 30, 30, 50));
        assert_eq!((p1.welch_p, p1.welch_q), (None, None));
        assert_eq!(comparison.samples_a, ["c1"]);
        assert_eq!(comparison.primary, PrimaryTest::Fisher);
    }

    #[test]
    fn replicates_get_a_welch_test() {
        let q = quantify(
            vec![
                sample("c1", "control"),
                sample("c2", "control"),
                sample("t1", "treated"),
                sample("t2", "treated"),
            ],
            vec![
                protein(1, None, &[5, 6, 20, 22]),
                protein(2, None, &[50, 50, 50, 50]),
            ],
        )
        .unwrap();
        let comparison = q.compare("control", "treated").unwrap();
        assert_eq!(comparison.primary, PrimaryTest::Welch);
        let changed = comparison.rows.iter().find(|r| r.group_id == 1).unwrap();
        assert!(changed.welch_p.is_some() && changed.welch_q.is_some());
        assert!(changed.log2_fold > 1.0);
    }

    #[test]
    fn samples_without_spectra_are_left_out_of_means() {
        // t2 counted nothing: treated is t1 alone (mean 20 normalized, not 10).
        let q = quantify(
            vec![
                sample("c1", "control"),
                sample("c2", "control"),
                sample("t1", "treated"),
                sample("t2", "treated"),
            ],
            vec![
                protein(1, None, &[10, 10, 20, 0]),
                protein(2, None, &[10, 10, 20, 0]),
            ],
        )
        .unwrap();
        let comparison = q.compare("control", "treated").unwrap();
        assert_eq!(comparison.samples_b, ["t1"]);
        assert_eq!(comparison.primary, PrimaryTest::Fisher);
        let all_empty = quantify(
            vec![sample("c1", "control"), sample("t1", "treated")],
            vec![protein(1, None, &[5, 0])],
        )
        .unwrap();
        assert_eq!(
            all_empty.compare("control", "treated").unwrap_err(),
            QuantError::NoSpectra("treated".into())
        );
    }

    #[test]
    fn categories_must_exist_and_differ() {
        let q = two_samples();
        assert_eq!(q.categories(), ["control", "treated"]);
        let error = q.compare("control", "knockout").unwrap_err();
        assert_eq!(
            error.to_string(),
            "no biosample is in category 'knockout' (categories: control, treated)"
        );
        assert!(matches!(
            q.compare("control", "control"),
            Err(QuantError::SameCategory(_))
        ));
        let misaligned = quantify(vec![sample("s", "x")], vec![protein(1, None, &[1, 2])]);
        assert!(matches!(misaligned, Err(QuantError::Misaligned { .. })));
    }
}
