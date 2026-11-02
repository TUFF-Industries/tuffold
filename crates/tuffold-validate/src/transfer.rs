//! Match-between-runs transfers, judged by target-decoy competition.
//!
//! A precursor identified in other runs but not in this one is looked for
//! where it should elute here (aligned retention time) at its m/z; whatever
//! peak is found there may be it or something else. Each such *target*
//! transfer gets a paired *decoy*: the same m/z and charge looked for at the
//! wrong time, borrowed from another identified precursor's place in the run
//! (FlashLFQ 2.0's construction, Solivais et al. 2024), so it lands on real,
//! co-eluting signal the way a wrong transfer does. Decoys stand in for wrong
//! transfers as decoy peptides do for wrong matches.
//!
//! Peaks are scored by a model trained only on *known* transfers: the run's
//! own MS/MS-identified precursors, looked for at the time the other runs
//! predict (positives) and at a borrowed time (negatives), as IonQuant trains
//! on ions with MS/MS evidence (Yu et al. 2021). The transfers judged are never
//! trained on. Each target competes with its decoy (the higher score stands,
//! a missing peak scores nothing), and q-values follow over the winners, the
//! same competition and +1 correction as every other level ([`crate::fdr`]).

use crate::fdr;
use crate::features::Features;
use crate::logistic::Logistic;

/// What a transferred peak is judged by, in this order: retention-time error
/// (in units of the alignment's spread), m/z error from the run's own
/// identified peaks (ppm), isotope error (|ln of observed over expected
/// M+1/M|), log10 apex intensity, scans in the peak, and intensity error
/// (|log10 apex intensity − the donors' apex intensity, put on this run's
/// scale|).
pub const TRANSFER_FEATURES: [&str; 6] = [
    "rt_error",
    "ppm_error",
    "isotope_error",
    "log_intensity",
    "points",
    "intensity_error",
];

/// Fewest known positives and negatives (each) the model is trained on;
/// below that, transfers are ranked by retention-time error alone.
pub const MIN_KNOWN: usize = 30;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransferKind {
    /// The run's own identified precursor, looked for where the other runs
    /// predict it.
    Known,
    /// A known precursor looked for at a borrowed time.
    KnownDecoy,
    /// A precursor not identified in the run, looked for where predicted.
    Target,
    /// A target looked for at a borrowed time.
    Decoy,
}

/// One looked-for peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransferRow {
    /// Pairs a target with its decoy (and a known with its known decoy).
    pub pair: usize,
    pub kind: TransferKind,
    /// [`TRANSFER_FEATURES`]; `None` when no peak was found.
    pub features: Option<[f64; 6]>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Judgement {
    /// Per row: the q-value of a target that won its pair; `None` for
    /// everything else (decoys, known rows, targets without a peak or beaten
    /// by their decoy).
    pub q: Vec<Option<f64>>,
    /// Whether a model was trained (else retention-time error alone ranks).
    pub trained: bool,
    /// Known positives and negatives with a peak.
    pub known: (usize, usize),
    /// Pairs won by targets and by decoys.
    pub won: (usize, usize),
}

/// Score every row with a peak and compete each target with its decoy.
pub fn judge(rows: &[TransferRow]) -> Judgement {
    let measured: Vec<usize> = (0..rows.len())
        .filter(|&i| rows[i].features.is_some())
        .collect();
    let columns = TRANSFER_FEATURES
        .iter()
        .enumerate()
        .map(|(c, name)| {
            let values = measured
                .iter()
                .map(|&i| rows[i].features.map_or(0.0, |f| f[c]))
                .collect();
            (name.to_string(), values)
        })
        .collect();
    let features = Features::from_columns(columns);
    let of_kind = |kind: TransferKind| -> Vec<usize> {
        (0..measured.len())
            .filter(|&j| rows[measured[j]].kind == kind)
            .collect()
    };
    let (positives, negatives) = (
        of_kind(TransferKind::Known),
        of_kind(TransferKind::KnownDecoy),
    );
    let model = (positives.len() >= MIN_KNOWN && negatives.len() >= MIN_KNOWN)
        .then(|| Logistic::fit(&features, &positives, &negatives))
        .flatten();
    let mut score: Vec<Option<f64>> = vec![None; rows.len()];
    for (j, &i) in measured.iter().enumerate() {
        score[i] = Some(match &model {
            Some(model) => model.score(&features, j),
            None => -features.value(0, j),
        });
    }

    // Each pair's winner among its target and decoy (a tie goes to the
    // decoy: never let a target win by default).
    let mut best: std::collections::BTreeMap<usize, (f64, bool, usize)> = Default::default();
    for (i, row) in rows.iter().enumerate() {
        let is_decoy = match row.kind {
            TransferKind::Target => false,
            TransferKind::Decoy => true,
            TransferKind::Known | TransferKind::KnownDecoy => continue,
        };
        let Some(s) = score[i] else { continue };
        let entry = best.entry(row.pair).or_insert((s, is_decoy, i));
        if s > entry.0 || (s == entry.0 && is_decoy) {
            *entry = (s, is_decoy, i);
        }
    }
    let winners: Vec<(f64, bool, usize)> = best.into_values().collect();
    let scores: Vec<f64> = winners.iter().map(|w| w.0).collect();
    let decoys: Vec<bool> = winners.iter().map(|w| w.1).collect();
    let q = fdr::q_values(&scores, &decoys);
    let mut out = vec![None; rows.len()];
    for (w, q) in winners.iter().zip(q) {
        if !w.1 {
            out[w.2] = Some(q);
        }
    }
    let won_by_decoys = decoys.iter().filter(|&&d| d).count();
    Judgement {
        q: out,
        trained: model.is_some(),
        known: (positives.len(), negatives.len()),
        won: (decoys.len() - won_by_decoys, won_by_decoys),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(pair: usize, kind: TransferKind, rt_error: Option<f64>) -> TransferRow {
        TransferRow {
            pair,
            kind,
            features: rt_error.map(|e| [e, 1.0, 0.1, 6.0, 8.0, 0.2]),
        }
    }

    #[test]
    fn without_enough_known_transfers_rt_error_ranks() {
        use TransferKind::{Decoy, Target};
        // Pair 0: target 0.1 beats decoy 2.0. Pair 1: decoy 0.5 beats target
        // 3.0. Pair 2: target 0.2, decoy found nothing. Pair 3: target found
        // nothing (no entry from it). Winners by score (−rt error): T −0.1,
        // T −0.2, D −0.5. FDR at each: (0 + 1)/1 = 1, (0 + 1)/2 = 0.5,
        // (1 + 1)/2 = 1; monotone from the worst up: 1, then 0.5, then
        // min(0.5, 1) = 0.5. So pairs 0 and 2's targets get 0.5.
        let rows = [
            row(0, Target, Some(0.1)),
            row(0, Decoy, Some(2.0)),
            row(1, Target, Some(3.0)),
            row(1, Decoy, Some(0.5)),
            row(2, Target, Some(0.2)),
            row(2, Decoy, None),
            row(3, Target, None),
            row(3, Decoy, None),
        ];
        let judged = judge(&rows);
        assert!(!judged.trained);
        assert_eq!(
            judged.q,
            [Some(0.5), None, None, None, Some(0.5), None, None, None]
        );
        assert_eq!(judged.won, (2, 1));
    }

    #[test]
    fn known_transfers_train_the_model_and_are_never_judged() {
        use TransferKind::{Decoy, Known, KnownDecoy, Target};
        // Known positives sit near the predicted time (rt error 0.0..0.4),
        // known decoys far off (3.0..4.0): the model learns that a small
        // error means a right peak. Then a near target beats its far decoy
        // and a far target loses to its near decoy.
        let mut rows = Vec::new();
        for k in 0..40 {
            let f = f64::from(k) / 100.0;
            rows.push(row(100 + k as usize, Known, Some(f)));
            rows.push(row(100 + k as usize, KnownDecoy, Some(3.0 + f / 0.4)));
        }
        rows.push(row(0, Target, Some(0.05)));
        rows.push(row(0, Decoy, Some(3.5)));
        rows.push(row(1, Target, Some(3.2)));
        rows.push(row(1, Decoy, Some(0.1)));
        let judged = judge(&rows);
        assert!(judged.trained);
        assert_eq!(judged.known, (40, 40));
        assert_eq!(judged.won, (1, 1));
        let n = rows.len();
        // Known rows never get a q-value; pair 0's target won, pair 1's lost.
        assert!(judged.q[..n - 4].iter().all(Option::is_none));
        assert!(judged.q[n - 4].is_some());
        assert_eq!(judged.q[n - 2], None);
    }
}
