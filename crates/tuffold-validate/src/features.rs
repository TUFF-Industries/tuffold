//! Per-engine feature matrices for rescoring.
//!
//! Features are the engine's own scores (only those reported on every PSM of
//! the engine, so nothing is imputed) followed by derived ones Percolator also
//! uses: precursor charge, absolute mass error, peptide length, missed
//! cleavages, and modification count. E-values and p-values are compared on
//! a `-log10` scale. Constant columns carry no information and are dropped.
//! Significance thresholds are features but not engine scores: an engine is
//! never ranked by one alone. The retention-time error
//! ([`crate::retention`]) is added by the caller, as it needs the confident
//! PSMs these features rank.

use std::collections::BTreeMap;

use crate::input::{PeptideRow, PsmRow};

/// Floor for `-log10` of e-values and p-values reported as 0.
const SMALLEST_PROBABILITY: f64 = 1e-300;
/// Charges at or above this share one indicator.
const CHARGE_CAP: i32 = 5;

/// Column-major feature values for one engine's PSMs (rows are local).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Features {
    pub names: Vec<String>,
    /// Columns `0..engine_scores` are engine scores, usable on their own to
    /// rank PSMs; the rest are derived.
    pub engine_scores: usize,
    pub columns: Vec<Vec<f64>>,
    pub rows: usize,
}

impl Features {
    pub fn build(rows: &[&PsmRow], peptides: &[&PeptideRow], score_names: &[String]) -> Self {
        let mut features = Self {
            names: Vec::new(),
            engine_scores: 0,
            columns: Vec::new(),
            rows: rows.len(),
        };
        let (thresholds, scores): (Vec<_>, Vec<_>) = engine_scores(rows, score_names)
            .into_iter()
            .partition(|(name, _)| is_threshold(name));
        for (name, column) in scores {
            features.push(name, column);
        }
        features.engine_scores = features.columns.len();
        for (name, column) in thresholds.into_iter().chain(derived(rows, peptides)) {
            features.push(name, column);
        }
        features
    }

    /// Features from ready columns: the first is the one score that may rank
    /// rows alone, the rest are derived.
    pub fn from_columns(columns: Vec<(String, Vec<f64>)>) -> Self {
        let mut columns = columns.into_iter();
        let mut features = Self {
            names: Vec::new(),
            engine_scores: 0,
            columns: Vec::new(),
            rows: 0,
        };
        if let Some((name, column)) = columns.next() {
            features.rows = column.len();
            features.push(name, column);
        }
        features.engine_scores = features.columns.len();
        for (name, column) in columns {
            features.push(name, column);
        }
        features
    }

    pub fn value(&self, column: usize, row: usize) -> f64 {
        self.columns[column][row]
    }

    /// Adds a derived column (kept only when finite and varying).
    pub fn add_derived(&mut self, name: &str, column: Vec<f64>) {
        self.push(name.to_string(), column);
    }

    /// Keeps a column only when it is finite and varies.
    fn push(&mut self, name: String, column: Vec<f64>) {
        let finite = column.iter().all(|value| value.is_finite());
        let varies = column.windows(2).any(|pair| pair[0] != pair[1]);
        if finite && varies {
            self.names.push(name);
            self.columns.push(column);
        }
    }
}

/// Scores reported on every row, by name, on a comparable scale.
fn engine_scores(rows: &[&PsmRow], score_names: &[String]) -> Vec<(String, Vec<f64>)> {
    let mut seen: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    for row in rows {
        for (name, _) in &row.scores {
            let entry = seen
                .entry(score_names[*name].as_str())
                .or_insert((*name, 0));
            entry.1 += 1;
        }
    }
    seen.into_iter()
        .filter(|(_, (_, count))| *count == rows.len())
        .map(|(name, (index, _))| {
            let probability_like = is_probability_like(name);
            let column = rows
                .iter()
                .map(|row| {
                    let value = row
                        .scores
                        .iter()
                        .find(|(score, _)| *score == index)
                        .map_or(f64::NAN, |(_, value)| *value);
                    if probability_like {
                        neg_log10_probability(value)
                    } else {
                        value
                    }
                })
                .collect();
            let label = if probability_like {
                format!("-log10({name})")
            } else {
                name.to_string()
            };
            (label, column)
        })
        .collect()
}

/// `-log10` of an e-value or p-value; 0 is floored at
/// [`SMALLEST_PROBABILITY`]. A negative (or NaN) value is not a probability
/// and gives NaN, which drops the whole column rather than scoring it as
/// perfect.
fn neg_log10_probability(value: f64) -> f64 {
    if value >= 0.0 {
        -value.max(SMALLEST_PROBABILITY).log10()
    } else {
        f64::NAN
    }
}

/// A significance threshold the engine computed for the spectrum (Mascot's
/// identity and homology thresholds): it says how hard the spectrum is, not
/// how good the match is, so it helps a model but never ranks on its own
/// (on a Scaffold file it out-ranked the ion score, by selection, not
/// merit).
fn is_threshold(name: &str) -> bool {
    name.to_lowercase().ends_with("threshold")
}

/// E-values and p-values span many orders of magnitude.
fn is_probability_like(name: &str) -> bool {
    let lowered = name.to_lowercase();
    ["expect", "evalue", "e-value", "pvalue", "p-value"]
        .iter()
        .any(|marker| lowered.contains(marker))
}

fn derived(rows: &[&PsmRow], peptides: &[&PeptideRow]) -> Vec<(String, Vec<f64>)> {
    let mut columns = Vec::new();
    let mut charges: Vec<i32> = rows.iter().map(|row| row.charge.min(CHARGE_CAP)).collect();
    charges.sort_unstable();
    charges.dedup();
    if charges.len() > 1 {
        for charge in charges {
            let name = if charge == CHARGE_CAP {
                format!("charge {CHARGE_CAP}+")
            } else {
                format!("charge {charge}")
            };
            let column = rows
                .iter()
                .map(|row| f64::from(u8::from(row.charge.min(CHARGE_CAP) == charge)))
                .collect();
            columns.push((name, column));
        }
    }
    let mass_errors: Option<Vec<f64>> = rows.iter().map(|row| mass_error_ppm(row)).collect();
    if let Some(column) = mass_errors {
        columns.push(("|mass error| (ppm)".to_string(), column));
    }
    columns.push((
        "peptide length".to_string(),
        peptides.iter().map(|p| p.sequence.len() as f64).collect(),
    ));
    columns.push((
        "missed cleavages".to_string(),
        peptides
            .iter()
            .map(|p| f64::from(missed_cleavages(&p.sequence)))
            .collect(),
    ));
    columns.push((
        "modifications".to_string(),
        peptides
            .iter()
            .map(|p| f64::from(p.modifications))
            .collect(),
    ));
    columns
}

/// Spacing of a precursor's isotope peaks in neutral mass (13C - 12C).
const ISOTOPE_SPACING: f64 = 1.003_354_835;
/// Isotope peaks engines may have picked instead of the monoisotopic one.
const ISOTOPE_OFFSETS: std::ops::RangeInclusive<i32> = -1..=3;

/// |mass error| in ppm against the nearest isotope peak: engines search
/// isotope errors (picking the 13C peak is common), and a +1 isotope match
/// is no worse an identification for being ~1 Da off.
fn mass_error_ppm(row: &PsmRow) -> Option<f64> {
    let calculated = row.calculated_mz.filter(|mz| *mz > 0.0)?;
    let experimental = row.experimental_mz?;
    let charge = f64::from(row.charge.max(1));
    ISOTOPE_OFFSETS
        .map(|k| {
            let shifted = calculated + f64::from(k) * ISOTOPE_SPACING / charge;
            ((experimental - shifted) / calculated * 1e6).abs()
        })
        .min_by(f64::total_cmp)
}

/// Internal K or R not followed by P (trypsin's rule).
pub(crate) fn missed_cleavages(sequence: &str) -> u32 {
    let residues = sequence.as_bytes();
    let count = residues
        .windows(2)
        .filter(|pair| matches!(pair[0], b'K' | b'R') && pair[1] != b'P')
        .count();
    u32::try_from(count).unwrap_or(u32::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_errors_measure_to_the_nearest_isotope_peak() {
        let row = |experimental: f64, charge: i32| PsmRow {
            experimental_mz: Some(experimental),
            ..psm(charge, Vec::new())
        };
        // 500.001 is 2 ppm off the monoisotopic peak.
        let direct = mass_error_ppm(&row(500.001, 2)).unwrap();
        assert!((direct - 2.0).abs() < 1e-6, "{direct}");
        // The +1 isotope at charge 2 sits 1.003354835 / 2 above: 2 ppm off it.
        let shifted = mass_error_ppm(&row(500.501_677_417_5 + 0.001, 2)).unwrap();
        assert!((shifted - 2.0).abs() < 1e-6, "{shifted}");
        // Halfway between peaks is far from both.
        assert!(mass_error_ppm(&row(500.25, 2)).unwrap() > 400.0);
        let unknown = PsmRow {
            experimental_mz: None,
            ..psm(2, Vec::new())
        };
        assert_eq!(mass_error_ppm(&unknown), None);
    }

    fn psm(charge: i32, scores: Vec<(usize, f64)>) -> PsmRow {
        PsmRow {
            id: 0,
            spectrum_id: 0,
            peptide_id: 0,
            engine: 0,
            charge,
            calculated_mz: Some(500.0),
            experimental_mz: Some(500.001),
            run_id: 1,
            retention_time: None,
            scores,
        }
    }

    fn peptide(sequence: &str, modifications: u32) -> PeptideRow {
        PeptideRow {
            id: 0,
            sequence: sequence.to_string(),
            modifications,
            is_decoy: false,
        }
    }

    #[test]
    fn thresholds_are_features_but_never_engine_scores() {
        let names = ["Mascot:score", "Mascot:homology threshold"].map(String::from);
        let rows = [
            psm(2, vec![(0, 40.0), (1, 30.0)]),
            psm(3, vec![(0, 20.0), (1, 25.0)]),
        ];
        let rows: Vec<&PsmRow> = rows.iter().collect();
        let peptides = [peptide("PEPTIDEK", 0), peptide("PEPTIDER", 1)];
        let peptides: Vec<&PeptideRow> = peptides.iter().collect();
        let features = Features::build(&rows, &peptides, &names);
        assert_eq!(features.engine_scores, 1);
        assert_eq!(features.names[0], "Mascot:score");
        assert_eq!(features.names[1], "Mascot:homology threshold");
    }

    #[test]
    fn missed_cleavages_follow_the_trypsin_rule() {
        assert_eq!(
            missed_cleavages("PEPTIDEK"),
            0,
            "C-terminal K is the cleavage"
        );
        assert_eq!(missed_cleavages("PEPKTIDER"), 1);
        assert_eq!(missed_cleavages("PEPKPTIDER"), 0, "KP is not cleaved");
        assert_eq!(missed_cleavages("KRK"), 2);
    }

    #[test]
    fn only_scores_on_every_psm_become_features_and_e_values_are_logged() {
        let names = vec![
            "xcorr".to_string(),
            "expect".to_string(),
            "partial".to_string(),
        ];
        let rows = [
            psm(2, vec![(0, 2.5), (1, 1e-3), (2, 7.0)]),
            psm(3, vec![(0, 1.5), (1, 10.0)]),
        ];
        let peptides = [peptide("PEPTIDEK", 0), peptide("PEPKTIDER", 1)];
        let features = Features::build(
            &rows.iter().collect::<Vec<_>>(),
            &peptides.iter().collect::<Vec<_>>(),
            &names,
        );
        assert_eq!(features.engine_scores, 2);
        assert_eq!(
            features.names,
            [
                "-log10(expect)",
                "xcorr",
                "charge 2",
                "charge 3",
                "peptide length",
                "missed cleavages",
                "modifications",
            ],
            "the mass error is constant (2 ppm both) and dropped"
        );
        // -log10(1e-3) = 3, -log10(10) = -1.
        assert!((features.value(0, 0) - 3.0).abs() < 1e-12);
        assert!((features.value(0, 1) + 1.0).abs() < 1e-12);
        assert_eq!(features.columns[2], [1.0, 0.0]);
        assert_eq!(features.columns[4], [8.0, 9.0]);
        assert_eq!(features.columns[5], [0.0, 1.0]);
    }

    #[test]
    fn invalid_e_values_drop_the_column_instead_of_scoring_perfectly() {
        assert_eq!(neg_log10_probability(0.0), 300.0);
        assert_eq!(neg_log10_probability(0.01), 2.0);
        assert!(neg_log10_probability(-1e-5).is_nan());
        assert!(neg_log10_probability(f64::NAN).is_nan());
        let names = vec!["expect".to_string()];
        let rows = [psm(2, vec![(0, -0.5)]), psm(2, vec![(0, 0.01)])];
        let peptides = [peptide("AK", 0), peptide("AK", 0)];
        let features = Features::build(
            &rows.iter().collect::<Vec<_>>(),
            &peptides.iter().collect::<Vec<_>>(),
            &names,
        );
        assert_eq!(features.engine_scores, 0);
    }

    #[test]
    fn non_finite_scores_are_not_features() {
        let names = vec!["score".to_string()];
        let rows = [psm(2, vec![(0, f64::INFINITY)]), psm(2, vec![(0, 1.0)])];
        let peptides = [peptide("AK", 0), peptide("AK", 0)];
        let features = Features::build(
            &rows.iter().collect::<Vec<_>>(),
            &peptides.iter().collect::<Vec<_>>(),
            &names,
        );
        assert_eq!(features.engine_scores, 0);
    }
}
