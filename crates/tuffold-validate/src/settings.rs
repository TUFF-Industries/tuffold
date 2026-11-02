//! What the user chooses: whether to rescore, and the thresholds — as FDR
//! (q-values) or, Scaffold-style, as probabilities (1 − PEP).

/// How thresholds are read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ThresholdMode {
    /// Each level passes at q ≤ the threshold (an FDR).
    #[default]
    Fdr,
    /// Each level passes at probability (1 − PEP) ≥ the threshold. The FDR of
    /// what passes is then estimated from the decoys that pass too.
    Probability,
}

impl ThresholdMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fdr => "fdr",
            Self::Probability => "probability",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "fdr" => Some(Self::Fdr),
            "probability" => Some(Self::Probability),
            _ => None,
        }
    }
}

/// A level thresholds apply to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Level {
    Psm,
    Peptide,
    Protein,
}

impl Level {
    pub const ALL: [Level; 3] = [Self::Psm, Self::Peptide, Self::Protein];

    pub fn name(self) -> &'static str {
        match self {
            Self::Psm => "PSM",
            Self::Peptide => "peptide",
            Self::Protein => "protein",
        }
    }
}

/// Floating point slack so a probability of exactly 0.99 passes a 99%
/// threshold (1 − 0.99 is not exactly 0.01).
const SLACK: f64 = 1e-12;

/// What each level must meet to be reported.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Thresholds {
    pub mode: ThresholdMode,
    /// The FDR (fdr mode) or the probability (probability mode) per level.
    pub psm: f64,
    pub peptide: f64,
    pub protein: f64,
    /// Distinct peptide sequences (unique or razor) a protein group needs.
    pub min_peptides: u32,
}

impl Default for Thresholds {
    /// 1% FDR at every level.
    fn default() -> Self {
        Self {
            mode: ThresholdMode::Fdr,
            psm: 0.01,
            peptide: 0.01,
            protein: 0.01,
            min_peptides: 1,
        }
    }
}

impl Thresholds {
    /// The usual probability thresholds: 95% spectra and peptides, 99%
    /// protein groups.
    pub fn probability_default() -> Self {
        Self {
            mode: ThresholdMode::Probability,
            psm: 0.95,
            peptide: 0.95,
            protein: 0.99,
            min_peptides: 1,
        }
    }

    pub fn value(&self, level: Level) -> f64 {
        match level {
            Level::Psm => self.psm,
            Level::Peptide => self.peptide,
            Level::Protein => self.protein,
        }
    }

    /// The same thresholds with `level` set to `value`.
    pub fn with(self, level: Level, value: f64) -> Self {
        match level {
            Level::Psm => Self { psm: value, ..self },
            Level::Peptide => Self {
                peptide: value,
                ..self
            },
            Level::Protein => Self {
                protein: value,
                ..self
            },
        }
    }

    /// The largest q-value (fdr mode) or PEP (probability mode) that passes
    /// at `level`; stores filter on that column with this bound.
    pub fn bound(&self, level: Level) -> f64 {
        match self.mode {
            ThresholdMode::Fdr => self.value(level),
            ThresholdMode::Probability => 1.0 - self.value(level) + SLACK,
        }
    }

    /// Whether an entry with this q-value and PEP passes at `level`. Unknown
    /// values (an entry that lost its competition) never pass.
    pub fn passes(&self, level: Level, q_value: Option<f64>, pep: Option<f64>) -> bool {
        let measure = match self.mode {
            ThresholdMode::Fdr => q_value,
            ThresholdMode::Probability => pep,
        };
        measure.is_some_and(|value| value <= self.bound(level))
    }

    /// `Err` names the first threshold outside its range.
    pub fn check(&self) -> Result<(), String> {
        let kind = match self.mode {
            ThresholdMode::Fdr => "FDR",
            ThresholdMode::Probability => "probability",
        };
        for level in Level::ALL {
            let value = self.value(level);
            if !(value > 0.0 && value < 1.0) {
                return Err(format!(
                    "{} {kind} must be above 0% and below 100%, not {}",
                    level.name(),
                    format_rate(value)
                ));
            }
        }
        if self.min_peptides == 0 {
            return Err("protein groups need at least 1 peptide".to_string());
        }
        Ok(())
    }

    /// "1% FDR" or "99% probability" for `level`.
    pub fn describe(&self, level: Level) -> String {
        match self.mode {
            ThresholdMode::Fdr => format!("{} FDR", format_rate(self.value(level))),
            ThresholdMode::Probability => {
                format!("{} probability", format_rate(self.value(level)))
            }
        }
    }
}

/// A change to thresholds as a person (or assistant) asks for it: values by
/// FDR or by probability per level (`[psm, peptide, protein]`), optionally a
/// mode, and min peptides. Every front end goes through
/// [`Thresholds::adjusted`], so they all read requests the same way.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThresholdRequest {
    pub by: Option<ThresholdMode>,
    pub fdr: [Option<f64>; 3],
    pub probability: [Option<f64>; 3],
    pub min_peptides: Option<u32>,
}

impl Thresholds {
    /// These thresholds changed as `request` asks. FDR values and
    /// probability values pick the mode and can't be mixed; switching mode
    /// starts from that mode's defaults (an FDR is no probability), keeping
    /// min peptides. The result is checked.
    pub fn adjusted(self, request: &ThresholdRequest) -> Result<Thresholds, String> {
        let any = |values: &[Option<f64>; 3]| values.iter().any(Option::is_some);
        let by_values = match (any(&request.fdr), any(&request.probability)) {
            (true, true) => {
                return Err("give FDR thresholds or probability thresholds, not both".into());
            }
            (true, false) => Some(ThresholdMode::Fdr),
            (false, true) => Some(ThresholdMode::Probability),
            (false, false) => None,
        };
        let mode = match (request.by, by_values) {
            (Some(asked), Some(implied)) if asked != implied => {
                return Err(format!(
                    "thresholding by {} contradicts the {} values given",
                    asked.as_str(),
                    implied.as_str()
                ));
            }
            (asked, implied) => asked.or(implied).unwrap_or(self.mode),
        };
        let base = if mode == self.mode {
            self
        } else {
            let defaults = match mode {
                ThresholdMode::Fdr => Self::default(),
                ThresholdMode::Probability => Self::probability_default(),
            };
            Self {
                min_peptides: self.min_peptides,
                ..defaults
            }
        };
        let values = match mode {
            ThresholdMode::Fdr => request.fdr,
            ThresholdMode::Probability => request.probability,
        };
        let adjusted = Level::ALL.into_iter().zip(values).fold(
            base,
            |thresholds, (level, value)| match value {
                Some(value) => thresholds.with(level, value),
                None => thresholds,
            },
        );
        let adjusted = Self {
            min_peptides: request.min_peptides.unwrap_or(base.min_peptides),
            ..adjusted
        };
        adjusted.check()?;
        Ok(adjusted)
    }
}

/// A rate as a percentage for people: 0.01 -> "1%", 0.025 -> "2.5%".
pub fn format_rate(rate: f64) -> String {
    let text = format!("{:.3}", rate * 100.0);
    let trimmed = text.trim_end_matches('0').trim_end_matches('.');
    format!("{trimmed}%")
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    /// Learn a combined score per engine (semi-supervised, cross-validated);
    /// off ranks by each engine's best single score.
    pub rescore: bool,
    pub thresholds: Thresholds,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            rescore: true,
            thresholds: Thresholds::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rates_format_without_trailing_zeros() {
        assert_eq!(format_rate(0.01), "1%");
        assert_eq!(format_rate(0.025), "2.5%");
        assert_eq!(format_rate(0.001), "0.1%");
        assert_eq!(format_rate(1.0 / 3.0), "33.333%");
    }

    #[test]
    fn thresholds_must_be_proper_rates() {
        assert!(Thresholds::default().check().is_ok());
        assert!(Thresholds::probability_default().check().is_ok());
        for bad in [0.0, 1.0, -0.1, f64::NAN] {
            let thresholds = Thresholds::default().with(Level::Peptide, bad);
            let message = thresholds.check().unwrap_err();
            assert!(message.starts_with("peptide FDR"), "{message}");
            let message = Thresholds::probability_default()
                .with(Level::Protein, bad)
                .check()
                .unwrap_err();
            assert!(message.starts_with("protein probability"), "{message}");
        }
        let thresholds = Thresholds {
            min_peptides: 0,
            ..Thresholds::default()
        };
        assert!(thresholds.check().is_err());
    }

    #[test]
    fn requests_choose_the_mode() {
        let saved = Thresholds::default().with(Level::Protein, 0.05);
        let request = |f: fn(&mut ThresholdRequest)| {
            let mut request = ThresholdRequest::default();
            f(&mut request);
            saved.adjusted(&request)
        };
        assert_eq!(request(|_| {}), Ok(saved));
        assert_eq!(
            request(|r| r.fdr[0] = Some(0.02)),
            Ok(saved.with(Level::Psm, 0.02))
        );
        // A probability switches mode, from that mode's defaults.
        let switched = request(|r| {
            r.probability[2] = Some(0.9);
            r.min_peptides = Some(2);
        })
        .unwrap();
        assert_eq!(switched.mode, ThresholdMode::Probability);
        assert_eq!(
            (switched.psm, switched.protein, switched.min_peptides),
            (0.95, 0.9, 2)
        );
        assert_eq!(
            request(|r| r.by = Some(ThresholdMode::Probability)),
            Ok(Thresholds::probability_default())
        );
        assert!(
            request(|r| {
                r.fdr[0] = Some(0.01);
                r.probability[2] = Some(0.99);
            })
            .is_err()
        );
        assert!(
            request(|r| {
                r.by = Some(ThresholdMode::Fdr);
                r.probability[2] = Some(0.99);
            })
            .is_err()
        );
        assert!(request(|r| r.fdr[1] = Some(1.5)).is_err(), "checked");
    }

    #[test]
    fn modes_read_q_values_or_peps() {
        let fdr = Thresholds::default();
        assert!(fdr.passes(Level::Psm, Some(0.01), Some(0.9)));
        assert!(!fdr.passes(Level::Psm, Some(0.011), Some(0.0)));
        let probability = Thresholds::probability_default();
        // 99% protein probability: PEP up to 1%, exactly 0.01 included.
        assert!(probability.passes(Level::Protein, Some(0.5), Some(0.01)));
        assert!(!probability.passes(Level::Protein, Some(0.0), Some(0.0101)));
        assert!(probability.passes(Level::Peptide, None, Some(0.05)));
        // Unknown never passes.
        assert!(!probability.passes(Level::Protein, Some(0.0), None));
        assert!(!fdr.passes(Level::Protein, None, Some(0.0)));
        assert_eq!(probability.describe(Level::Protein), "99% probability");
        assert_eq!(fdr.describe(Level::Psm), "1% FDR");
        assert_eq!(
            ThresholdMode::parse("probability"),
            Some(ThresholdMode::Probability)
        );
    }
}
