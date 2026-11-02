//! Peptide scores that pool a peptide's spectra (iProphet's idea): a correct
//! peptide tends to be matched again, by other spectra and at other charges,
//! while a random match rarely repeats.
//!
//! Pooling happens at the peptide level only. At the spectrum level it
//! breaks target-decoy FDR control: a wrong spectrum that lands on a peptide
//! other spectra identify correctly is boosted by them, and no decoy can be
//! (decoy peptides are never present), so decoys stop modelling the errors.
//! On the synthetic searches that raised the true FDP at 1% from 1.2% to
//! 3.3%. At the peptide level such a spectrum adds no wrong peptide, and a
//! peptide that is absent gets only random matches, as a decoy does.
//!
//! The features per peptide are its best and second-best spectrum scores
//! and how many other spectra and charges back it; they are ranked by the
//! same semi-supervised, cross-validated rescoring engines get
//! ([`crate::rescore`]), which falls back to the best spectrum score (the
//! peptide score without pooling) whenever pooling does not accept more
//! peptides.
//!
//! Modified forms are not pooled: a wrong form of a peptide that is present
//! (a misplaced or misread modification) would be boosted by the right one,
//! and no decoy models that either. Sibling peptides of the same protein
//! are, but decoy-safe ([`crate::siblings`]).

use std::collections::{BTreeMap, BTreeSet};

use crate::features::Features;
use crate::rescore::{Method, Rows, score_engine};
use crate::siblings::Siblings;

/// One spectrum's vote for a peptide.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Support {
    pub spectrum_id: i64,
    pub charge: i32,
    pub score: f64,
}

/// One peptide and every spectrum assigned to it.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PeptideSpectra {
    pub peptide_id: i64,
    pub is_decoy: bool,
    /// Best first.
    pub spectra: Vec<Support>,
}

/// A score per peptide (higher is better) and how it was made, given each
/// peptide's sibling peptides ([`crate::siblings`]). Without rescoring, or
/// when pooling cannot be trained or does not help, the score is the best
/// spectrum's.
pub(crate) fn peptide_scores(
    peptides: &[PeptideSpectra],
    siblings: Option<&Siblings>,
    rescore: bool,
) -> (Vec<f64>, Method) {
    let mut columns = columns(peptides);
    if let Some(siblings) = siblings {
        let name = if siblings.imputed.is_some() {
            SIBLINGS_IMPUTED
        } else {
            SIBLINGS
        };
        columns.push((
            name.to_string(),
            siblings
                .counts
                .iter()
                .map(|&n| (n as f64).ln_1p())
                .collect(),
        ));
    }
    let features = Features::from_columns(columns);
    let best: Vec<f64> = peptides.iter().map(|p| p.spectra[0].score).collect();
    if features.engine_scores == 0 {
        // Every peptide scores the same: there is nothing to rank by.
        let method = Method::BestScore {
            score: BEST.to_string(),
            higher_is_better: true,
            why: crate::rescore::BestScoreReason::OneFeature,
        };
        return (best, method);
    }
    let ids: Vec<i64> = peptides.iter().map(|p| p.peptide_id).collect();
    let decoys: Vec<bool> = peptides.iter().map(|p| p.is_decoy).collect();
    let rows = Rows {
        features: &features,
        spectra: &ids,
        folds: &ids,
        decoys: &decoys,
    };
    match score_engine(&rows, rescore) {
        (scores, method @ Method::Rescored { .. }) => (scores, method),
        // The best spectrum score, unchanged, so peptides rank exactly as
        // without pooling.
        (_, method) => (best, method),
    }
}

const BEST: &str = "best spectrum score";
/// The sibling-peptide feature's name, which reports look for.
pub const SIBLINGS: &str = "log(1 + sibling peptides)";
/// Its name when decoys named no target and their counts were imputed from
/// low-scoring targets' ([`crate::siblings`]).
pub const SIBLINGS_IMPUTED: &str = "log(1 + sibling peptides; decoys' imputed)";

fn columns(peptides: &[PeptideSpectra]) -> Vec<(String, Vec<f64>)> {
    // A peptide with one spectrum gets the lowest spectrum score anywhere as
    // its second: no second spectrum is the weakest evidence there is.
    let floor = peptides
        .iter()
        .flat_map(|p| p.spectra.iter().map(|s| s.score))
        .min_by(f64::total_cmp)
        .unwrap_or(0.0);
    let column = |value: &dyn Fn(&PeptideSpectra) -> f64| peptides.iter().map(value).collect();
    vec![
        (BEST.to_string(), column(&|p| p.spectra[0].score)),
        (
            "second-best spectrum score".to_string(),
            column(&|p| p.spectra.get(1).map_or(floor, |s| s.score)),
        ),
        (
            "log(1 + other spectra)".to_string(),
            column(&|p| (distinct(p, |s| s.spectrum_id).saturating_sub(1) as f64).ln_1p()),
        ),
        (
            "other charges".to_string(),
            column(&|p| distinct(p, |s| s.charge).saturating_sub(1) as f64),
        ),
    ]
}

fn distinct<T: Ord>(peptide: &PeptideSpectra, key: impl Fn(&Support) -> T) -> usize {
    peptide
        .spectra
        .iter()
        .map(key)
        .collect::<BTreeSet<T>>()
        .len()
}

/// Spectra grouped by peptide, peptides in id order, spectra best first
/// (ties by spectrum id, so the order never depends on the input's).
pub(crate) fn group_by_peptide(
    supports: impl IntoIterator<Item = (i64, bool, Support)>,
) -> Vec<PeptideSpectra> {
    let mut grouped: BTreeMap<i64, PeptideSpectra> = BTreeMap::new();
    for (peptide_id, is_decoy, support) in supports {
        grouped
            .entry(peptide_id)
            .or_insert_with(|| PeptideSpectra {
                peptide_id,
                is_decoy,
                spectra: Vec::new(),
            })
            .spectra
            .push(support);
    }
    let mut peptides: Vec<PeptideSpectra> = grouped.into_values().collect();
    for peptide in &mut peptides {
        peptide.spectra.sort_by(|a, b| {
            b.score
                .total_cmp(&a.score)
                .then(a.spectrum_id.cmp(&b.spectrum_id))
        });
    }
    peptides
}

#[cfg(test)]
mod tests {
    use super::*;

    fn support(spectrum_id: i64, charge: i32, score: f64) -> Support {
        Support {
            spectrum_id,
            charge,
            score,
        }
    }

    #[test]
    fn features_count_other_spectra_and_charges() {
        let peptides = group_by_peptide([
            (1, false, support(10, 2, 1.0)),
            (1, false, support(11, 3, 4.0)),
            (1, false, support(12, 2, 2.0)),
            (2, true, support(13, 2, -1.0)),
        ]);
        assert_eq!(
            peptides[0]
                .spectra
                .iter()
                .map(|s| s.spectrum_id)
                .collect::<Vec<_>>(),
            [11, 12, 10],
            "best first"
        );
        let columns = columns(&peptides);
        let names: Vec<&str> = columns.iter().map(|(name, _)| name.as_str()).collect();
        assert_eq!(
            names,
            [
                "best spectrum score",
                "second-best spectrum score",
                "log(1 + other spectra)",
                "other charges"
            ]
        );
        assert_eq!(columns[0].1, [4.0, -1.0]);
        // Peptide 2 has one spectrum: its second is the lowest score, -1.
        assert_eq!(columns[1].1, [2.0, -1.0]);
        assert_eq!(columns[2].1, [2_f64.ln_1p(), 0.0]);
        assert_eq!(columns[3].1, [1.0, 0.0]);
    }

    #[test]
    fn without_rescoring_peptides_keep_their_best_spectrum_score() {
        let peptides = group_by_peptide([
            (1, false, support(10, 2, 1.0)),
            (1, false, support(11, 3, 4.0)),
            (2, true, support(13, 2, -1.0)),
        ]);
        let (scores, method) = peptide_scores(&peptides, None, false);
        assert_eq!(scores, [4.0, -1.0]);
        assert!(matches!(method, Method::BestScore { .. }), "{method:?}");
    }
}
