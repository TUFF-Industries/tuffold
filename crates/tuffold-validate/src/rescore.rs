//! Semi-supervised rescoring of one engine's PSMs (Percolator / mokapot
//! style), with two linear models: L2-regularized logistic regression and
//! Fisher's discriminant. Both are trained in every fold, and the one whose
//! fold models accept more targets in the folds they never saw ranks the
//! engine, in every fold (see `logistic.rs` for why neither is always
//! better; accepted training targets favored logistic regression, which
//! fits that very boundary, even where it did worse on held-out data).
//!
//! 1. Rows are split into [`FOLDS`] folds by their peptide's sequence
//!    ([`Rows::folds`]), deterministically.
//! 2. For each fold, a model is trained on the other folds: start from the
//!    engine score that accepts the most targets at [`TRAIN_FDR`]; then
//!    [`ITERATIONS`] times, take those targets as positives and every decoy
//!    as negatives, fit the model, and re-rank. When fewer than
//!    [`MIN_POSITIVES`] targets pass, positives come from the first looser
//!    rung of [`POSITIVE_FDR_LADDER`] that has enough (weak scores may only
//!    reach 1% together, never alone). A model that accepts fewer training
//!    targets at [`TRAIN_FDR`] than the starting score is not used, and then
//!    no fold's is: folds ranked by different kinds of score (a trained
//!    model in one, an integer ion count in another) do not share one scale
//!    even after calibration, and merging them lost a fifth of the accepted
//!    spectra on real Sage results.
//! 3. Each fold is scored by the model it never saw, then calibrated so the
//!    folds share one scale: `(s - median) / MAD` of the fold's decoy scores
//!    (the standard deviation instead of the MAD, for every fold, when any
//!    fold's decoys mostly tie).
//!    Decoys model the null distribution and are plentiful, so this anchor is
//!    stable even in small folds. (mokapot anchors on the score accepted at
//!    1% FDR instead, which a fold of weak scores can fail to reach by
//!    chance, leaving nothing to calibrate against.)
//! 4. The merged, calibrated scores must accept at least as many targets at
//!    [`TRAIN_FDR`] as the best single score does over the same rows, or the
//!    best single score ranks the engine: rescoring may never cost
//!    identifications.
//!
//! When rescoring is off or cannot be trusted (one feature, too few accepted
//! targets, a fold that cannot be calibrated) the engine is ranked by its
//! best single score instead, and the reason is reported.

use crate::fdr::{q_values, winners};
use crate::features::Features;
use crate::lda::Lda;
use crate::logistic::Logistic;

pub(crate) const FOLDS: u64 = 3;
pub(crate) const TRAIN_FDR: f64 = 0.01;
const ITERATIONS: usize = 10;
/// FDRs tried in turn for positive training examples.
const POSITIVE_FDR_LADDER: [f64; 3] = [TRAIN_FDR, 0.05, 0.10];
/// Targets (whole engine) passing the loosest rung of the ladder under the
/// best single score, needed before a model is trained.
const MIN_CONFIDENT: usize = 30;
/// Decoys needed (whole engine) before a model is trained: about 50 per fold
/// for the discriminant's negatives and each fold's calibration.
const MIN_DECOYS: usize = 150;
/// Positives needed in one training set to fit an iteration.
const MIN_POSITIVES: usize = 10;
const FOLD_SEED: u64 = 0x5455_4646;

/// How an engine's PSMs were ranked.
#[derive(Debug, Clone, PartialEq)]
pub enum Method {
    Rescored {
        features: Vec<String>,
        model: RescoreModel,
    },
    BestScore {
        score: String,
        higher_is_better: bool,
        why: BestScoreReason,
    },
    Unscored {
        why: UnscoredReason,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BestScoreReason {
    RescoringOff,
    OneFeature,
    TooFewConfident,
    CalibrationFailed,
    ModelNotBetter,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnscoredReason {
    NoDecoys,
    NoSharedScore,
}

impl BestScoreReason {
    pub fn describe(self) -> String {
        match self {
            Self::RescoringOff => "rescoring is off".to_string(),
            Self::OneFeature => {
                "only one feature is available, so there is nothing to combine".to_string()
            }
            Self::TooFewConfident => format!(
                "fewer than {MIN_CONFIDENT} targets pass {}% FDR or fewer than {MIN_DECOYS} \
                 decoys exist, too few to train on",
                POSITIVE_FDR_LADDER[POSITIVE_FDR_LADDER.len() - 1] * 100.0
            ),
            Self::CalibrationFailed => {
                "a cross-validation fold's decoy scores are all identical, so the folds could \
                 not be put on one scale"
                    .to_string()
            }
            Self::ModelNotBetter => format!(
                "a trained model accepted fewer targets at {}% FDR than this score",
                TRAIN_FDR * 100.0
            ),
        }
    }
}

impl UnscoredReason {
    pub fn describe(self) -> &'static str {
        match self {
            Self::NoDecoys => {
                "it reported no decoy PSMs (the search had no decoys, or its output was already \
                 filtered), so its FDR cannot be estimated"
            }
            Self::NoSharedScore => "no score is reported on every one of its PSMs",
        }
    }
}

/// The kind of model that rescored an engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RescoreModel {
    Logistic,
    Discriminant,
}

impl RescoreModel {
    /// Tried in this order; a tie keeps the first.
    const ALL: [Self; 2] = [Self::Logistic, Self::Discriminant];

    /// "logistic regression", for reports.
    pub fn describe(self) -> &'static str {
        match self {
            Self::Logistic => "logistic regression",
            Self::Discriminant => "linear discriminant analysis",
        }
    }

    fn fit(self, features: &Features, positives: &[usize], negatives: &[usize]) -> Option<Model> {
        match self {
            Self::Logistic => Logistic::fit(features, positives, negatives).map(Model::Logistic),
            Self::Discriminant => Lda::fit(features, positives, negatives).map(Model::Discriminant),
        }
    }
}

enum Model {
    Logistic(Logistic),
    Discriminant(Lda),
}

impl Model {
    fn score(&self, features: &Features, row: usize) -> f64 {
        match self {
            Self::Logistic(model) => model.score(features, row),
            Self::Discriminant(model) => model.score(features, row),
        }
    }
}

/// Each fold's test rows and their raw model scores.
type FoldScores = Vec<(Vec<usize>, Vec<f64>)>;

/// One way to rank rows: a single engine score, oriented.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Direction {
    column: usize,
    sign: f64,
}

impl Direction {
    fn score(self, features: &Features, row: usize) -> f64 {
        self.sign * features.value(self.column, row)
    }
}

/// Rows of one engine with what competition and labels need.
pub(crate) struct Rows<'a> {
    pub features: &'a Features,
    /// The spectrum each row competes in.
    pub spectra: &'a [i64],
    /// What puts a row in its cross-validation fold ([`fold_of`]): for PSMs
    /// their peptide's sequence ([`sequence_fold`]), so a peptide's
    /// spectra are never split between training and testing. Split, a
    /// model trained on some spectra of a repeated wrong match (an absent
    /// analyte hitting one decoy again and again) learns to push down the
    /// held-out ones, and decoys stop counting such errors (Freestone,
    /// Noble & Keich 2023 showed Percolator losing FDR control this way).
    pub folds: &'a [i64],
    pub decoys: &'a [bool],
}

/// Every row's best single engine score (higher is better): the score and
/// direction that accept the most targets, as the engine is ranked without
/// rescoring. `None` without an engine score.
pub(crate) fn best_single_scores(rows: &Rows<'_>) -> Option<Vec<f64>> {
    let all: Vec<usize> = (0..rows.features.rows).collect();
    let (best, _) = best_direction(rows, &all)?;
    Some(
        all.iter()
            .map(|&row| best.score(rows.features, row))
            .collect(),
    )
}

/// Scores for every row (higher is better) and how they were made. The
/// features must include at least one engine score and the rows a decoy.
pub(crate) fn score_engine(rows: &Rows<'_>, rescore: bool) -> (Vec<f64>, Method) {
    let all: Vec<usize> = (0..rows.features.rows).collect();
    let (best, _) = best_direction(rows, &all).expect("an engine score");
    let fallback = |why: BestScoreReason| {
        let scores = all
            .iter()
            .map(|&row| best.score(rows.features, row))
            .collect();
        let method = Method::BestScore {
            score: rows.features.names[best.column].clone(),
            higher_is_better: best.sign > 0.0,
            why,
        };
        (scores, method)
    };
    if !rescore {
        return fallback(BestScoreReason::RescoringOff);
    }
    if rows.features.columns.len() < 2 {
        return fallback(BestScoreReason::OneFeature);
    }
    let decoys = rows.decoys.iter().filter(|&&decoy| decoy).count();
    let best_scores: Vec<f64> = all
        .iter()
        .map(|&row| best.score(rows.features, row))
        .collect();
    let confident = positives(rows, &all, &best_scores).len();
    if confident < MIN_CONFIDENT || decoys < MIN_DECOYS {
        return fallback(BestScoreReason::TooFewConfident);
    }
    let splits: Vec<(Vec<usize>, Vec<usize>)> = (0..FOLDS)
        .map(|fold| {
            all.iter()
                .partition(|&&row| fold_of(rows.folds[row]) != fold)
        })
        .collect();
    // Each kind of model in every fold; a kind with a fold that did not beat
    // its start is out, since folds ranked by different kinds of score do
    // not merge (module docs). Of the rest, the kind accepting the most
    // held-out targets wins.
    // The six trainings are independent: one thread each (the result does
    // not depend on scheduling).
    let trained: Vec<Vec<Option<Model>>> = std::thread::scope(|scope| {
        let handles: Vec<Vec<_>> = RescoreModel::ALL
            .iter()
            .map(|&kind| {
                splits
                    .iter()
                    .map(|(train, _)| scope.spawn(move || train_fold(rows, train, kind)))
                    .collect()
            })
            .collect();
        handles
            .into_iter()
            .map(|kind| {
                kind.into_iter()
                    .map(|handle| {
                        handle
                            .join()
                            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
                    })
                    .collect()
            })
            .collect()
    });
    let mut chosen: Option<(RescoreModel, FoldScores, usize)> = None;
    for (kind, models) in RescoreModel::ALL.into_iter().zip(trained) {
        let Some(models) = models.into_iter().collect::<Option<Vec<Model>>>() else {
            continue;
        };
        let mut folds = Vec::new();
        let mut held_out = 0;
        for ((_, test), model) in splits.iter().zip(&models) {
            let raw: Vec<f64> = test
                .iter()
                .map(|&row| model.score(rows.features, row))
                .collect();
            held_out += accepted(rows, test, &raw).len();
            folds.push((test.clone(), raw));
        }
        if chosen.as_ref().is_none_or(|(_, _, best)| held_out > *best) {
            chosen = Some((kind, folds, held_out));
        }
    }
    let Some((kind, folds, _)) = chosen else {
        return fallback(BestScoreReason::ModelNotBetter);
    };
    // Every fold must use the same spread measure, or their scales differ.
    let spread = if folds
        .iter()
        .all(|(test, raw)| decoy_scale(rows, test, raw, Spread::Mad).is_some())
    {
        Spread::Mad
    } else {
        Spread::StandardDeviation
    };
    let mut scores = vec![0.0; rows.features.rows];
    for (test, raw) in &folds {
        let Some((center, scale)) = decoy_scale(rows, test, raw, spread) else {
            return fallback(BestScoreReason::CalibrationFailed);
        };
        for (row, score) in test.iter().zip(raw) {
            scores[*row] = (score - center) / scale;
        }
    }
    if accepted(rows, &all, &scores).len() < accepted(rows, &all, &best_scores).len() {
        return fallback(BestScoreReason::ModelNotBetter);
    }
    let method = Method::Rescored {
        features: rows.features.names.clone(),
        model: kind,
    };
    (scores, method)
}

/// A model of `kind` trained on `train`, or `None` when it accepts fewer
/// training targets than the starting score.
fn train_fold(rows: &Rows<'_>, train: &[usize], kind: RescoreModel) -> Option<Model> {
    let Some((start, start_accepted)) = best_direction(rows, train) else {
        unreachable!("score_engine checked for an engine score")
    };
    let mut scores: Vec<f64> = train
        .iter()
        .map(|&row| start.score(rows.features, row))
        .collect();
    let negatives: Vec<usize> = train
        .iter()
        .copied()
        .filter(|&row| rows.decoys[row])
        .collect();
    let mut model = None;
    for _ in 0..ITERATIONS {
        let positives = positives(rows, train, &scores);
        if positives.len() < MIN_POSITIVES {
            break;
        }
        let Some(fitted) = kind.fit(rows.features, &positives, &negatives) else {
            break;
        };
        scores = train
            .iter()
            .map(|&row| fitted.score(rows.features, row))
            .collect();
        model = Some(fitted);
    }
    model.filter(|_| accepted(rows, train, &scores).len() >= start_accepted)
}

/// Targets confident by the engine's best single score: the positives
/// rescoring starts from (the first rung of [`POSITIVE_FDR_LADDER`] with
/// enough). Empty without an engine score.
pub(crate) fn confident_by_best_score(rows: &Rows<'_>) -> Vec<usize> {
    let all: Vec<usize> = (0..rows.features.rows).collect();
    let Some((best, _)) = best_direction(rows, &all) else {
        return Vec::new();
    };
    let scores: Vec<f64> = all
        .iter()
        .map(|&row| best.score(rows.features, row))
        .collect();
    positives(rows, &all, &scores)
}

/// The engine score (and orientation) accepting the most targets at
/// [`TRAIN_FDR`] over `subset`, ties broken by how well it separates targets
/// from decoys (AUC), then by column order. Also returns the accepted count.
fn best_direction(rows: &Rows<'_>, subset: &[usize]) -> Option<(Direction, usize)> {
    let mut best: Option<(Direction, usize, f64)> = None;
    for column in 0..rows.features.engine_scores {
        let separation = auc(rows, subset, column);
        for (sign, separation) in [(1.0, separation), (-1.0, 1.0 - separation)] {
            let direction = Direction { column, sign };
            let scores: Vec<f64> = subset
                .iter()
                .map(|&row| direction.score(rows.features, row))
                .collect();
            let count = accepted(rows, subset, &scores).len();
            let better = best.is_none_or(|(_, best_count, best_separation)| {
                count > best_count || (count == best_count && separation > best_separation)
            });
            if better {
                best = Some((direction, count, separation));
            }
        }
    }
    best.map(|(direction, count, _)| (direction, count))
}

/// Target rows in `subset` that win their spectrum and pass [`TRAIN_FDR`].
/// `scores` is parallel to `subset`.
fn accepted(rows: &Rows<'_>, subset: &[usize], scores: &[f64]) -> Vec<usize> {
    accepted_at(rows, subset, scores, TRAIN_FDR)
}

/// Positive training examples: the targets accepted at the first rung of
/// [`POSITIVE_FDR_LADDER`] with at least [`MIN_POSITIVES`] (or the loosest).
fn positives(rows: &Rows<'_>, subset: &[usize], scores: &[f64]) -> Vec<usize> {
    let competed = competed_q(rows, subset, scores);
    let mut found = Vec::new();
    for fdr in POSITIVE_FDR_LADDER {
        found = targets_within(rows, subset, &competed, fdr);
        if found.len() >= MIN_POSITIVES {
            break;
        }
    }
    found
}

fn accepted_at(rows: &Rows<'_>, subset: &[usize], scores: &[f64], fdr: f64) -> Vec<usize> {
    targets_within(rows, subset, &competed_q(rows, subset, scores), fdr)
}

fn targets_within(rows: &Rows<'_>, subset: &[usize], q: &[Option<f64>], fdr: f64) -> Vec<usize> {
    subset
        .iter()
        .zip(q)
        .filter(|&(&row, q)| !rows.decoys[row] && q.is_some_and(|q| q <= fdr))
        .map(|(&row, _)| row)
        .collect()
}

/// q-values after competition per spectrum, parallel to `subset`; `None`
/// for rows that lost their spectrum.
pub(crate) fn competed_q(rows: &Rows<'_>, subset: &[usize], scores: &[f64]) -> Vec<Option<f64>> {
    let spectra: Vec<i64> = subset.iter().map(|&row| rows.spectra[row]).collect();
    let decoys: Vec<bool> = subset.iter().map(|&row| rows.decoys[row]).collect();
    let kept = winners(&spectra, scores, &decoys);
    let kept_scores: Vec<f64> = kept.iter().map(|&i| scores[i]).collect();
    let kept_decoys: Vec<bool> = kept.iter().map(|&i| decoys[i]).collect();
    let mut q = vec![None; subset.len()];
    for (index, value) in kept.iter().zip(q_values(&kept_scores, &kept_decoys)) {
        q[*index] = Some(value);
    }
    q
}

/// Probability that a random target outscores a random decoy on `column`
/// (ties count half); 0.5 when either class is missing.
fn auc(rows: &Rows<'_>, subset: &[usize], column: usize) -> f64 {
    let mut ordered: Vec<usize> = subset.to_vec();
    ordered.sort_by(|&a, &b| {
        rows.features
            .value(column, a)
            .total_cmp(&rows.features.value(column, b))
    });
    let (mut rank_sum, mut targets) = (0.0, 0.0);
    let mut start = 0;
    while start < ordered.len() {
        let value = rows.features.value(column, ordered[start]);
        let end = start
            + ordered[start..]
                .iter()
                .take_while(|&&row| rows.features.value(column, row) == value)
                .count();
        // Average 1-based rank of the tie group.
        let rank = (start + 1 + end) as f64 / 2.0;
        for &row in &ordered[start..end] {
            if !rows.decoys[row] {
                rank_sum += rank;
                targets += 1.0;
            }
        }
        start = end;
    }
    let decoys = ordered.len() as f64 - targets;
    if targets == 0.0 || decoys == 0.0 {
        return 0.5;
    }
    (rank_sum - targets * (targets + 1.0) / 2.0) / (targets * decoys)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Spread {
    /// Median absolute deviation: robust, but 0 when most decoys tie.
    Mad,
    StandardDeviation,
}

/// Scores on the shared scale: `(s - median) / spread` of the decoy scores,
/// the spread being the MAD, or the standard deviation when more than half
/// the decoys tie (integer scores such as matched-ion counts make the MAD 0).
/// `None` without decoys or when their scores do not vary at all.
pub(crate) fn calibrate(rows: &Rows<'_>, test: &[usize], raw: &[f64]) -> Option<Vec<f64>> {
    let (center, scale) = decoy_scale(rows, test, raw, Spread::Mad)
        .or_else(|| decoy_scale(rows, test, raw, Spread::StandardDeviation))?;
    Some(raw.iter().map(|score| (score - center) / scale).collect())
}

/// `(median, spread)` of the decoy scores; `None` when the spread is 0.
fn decoy_scale(rows: &Rows<'_>, test: &[usize], raw: &[f64], spread: Spread) -> Option<(f64, f64)> {
    let mut decoy_scores: Vec<f64> = test
        .iter()
        .zip(raw)
        .filter(|(row, _)| rows.decoys[**row])
        .map(|(_, score)| *score)
        .collect();
    decoy_scores.sort_by(f64::total_cmp);
    let center = median(&decoy_scores)?;
    let scale = match spread {
        Spread::Mad => {
            let mut deviations: Vec<f64> =
                decoy_scores.iter().map(|s| (s - center).abs()).collect();
            deviations.sort_by(f64::total_cmp);
            median(&deviations)?
        }
        Spread::StandardDeviation => {
            let mean = decoy_scores.iter().sum::<f64>() / decoy_scores.len() as f64;
            let variance = decoy_scores.iter().map(|s| (s - mean).powi(2)).sum::<f64>()
                / decoy_scores.len() as f64;
            variance.sqrt()
        }
    };
    (scale > 0.0).then_some((center, scale))
}

fn median(sorted: &[f64]) -> Option<f64> {
    let middle = sorted.len() / 2;
    match sorted.len() {
        0 => None,
        n if n % 2 == 1 => Some(sorted[middle]),
        _ => Some((sorted[middle - 1] + sorted[middle]) / 2.0),
    }
}

/// A PSM's fold key: its peptide's sequence (I and L alike, modifications
/// ignored), hashed (FNV-1a).
pub(crate) fn sequence_fold(sequence: &str) -> i64 {
    let hash = sequence
        .bytes()
        .fold(0xcbf2_9ce4_8422_2325_u64, |hash, byte| {
            let byte = if byte == b'I' { b'L' } else { byte };
            (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
        });
    hash as i64
}

/// Deterministic fold for a fold key (SplitMix64 of it).
pub(crate) fn fold_of(key: i64) -> u64 {
    let mut z = (key as u64 ^ FOLD_SEED).wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)) % FOLDS
}

#[cfg(test)]
mod tests {
    use super::*;

    fn features(columns: Vec<Vec<f64>>, engine_scores: usize) -> Features {
        let rows = columns[0].len();
        Features {
            names: (0..columns.len()).map(|i| format!("f{i}")).collect(),
            engine_scores,
            columns,
            rows,
        }
    }

    #[test]
    fn auc_counts_pairs_and_half_ties() {
        // Targets 3, 1; decoys 2, 1. Pairs (t > d): (3,2) (3,1) = 2, tie (1,1)
        // = 0.5, (1,2) = 0: (2 + 0.5) / 4 = 0.625.
        let f = features(vec![vec![3.0, 1.0, 2.0, 1.0]], 1);
        let rows = Rows {
            features: &f,
            spectra: &[1, 2, 3, 4],
            folds: &[1, 2, 3, 4],
            decoys: &[false, false, true, true],
        };
        assert!((auc(&rows, &[0, 1, 2, 3], 0) - 0.625).abs() < 1e-12);
    }

    #[test]
    fn the_best_direction_may_be_lower_is_better() {
        // An e-value-like column: targets low, decoys high. Column 1 is noise.
        let mut good = Vec::new();
        let mut noise = Vec::new();
        let mut decoys = Vec::new();
        for i in 0..200 {
            let decoy = i % 4 == 0;
            good.push(if decoy {
                10.0 + f64::from(i)
            } else {
                f64::from(i) * 0.01
            });
            noise.push(f64::from(i % 7));
            decoys.push(decoy);
        }
        let f = features(vec![noise, good], 2);
        let spectra: Vec<i64> = (0..200).collect();
        let rows = Rows {
            features: &f,
            spectra: &spectra,
            folds: &spectra,
            decoys: &decoys,
        };
        let all: Vec<usize> = (0..200).collect();
        let (direction, accepted) = best_direction(&rows, &all).unwrap();
        assert_eq!(
            direction,
            Direction {
                column: 1,
                sign: -1.0
            }
        );
        assert_eq!(accepted, 150, "every target outranks every decoy");
    }

    #[test]
    fn folds_are_deterministic_and_balanced() {
        let mut counts = [0_usize; FOLDS as usize];
        for spectrum in 0..30_000 {
            counts[usize::try_from(fold_of(spectrum)).unwrap()] += 1;
        }
        assert_eq!(fold_of(42), fold_of(42));
        for count in counts {
            assert!((9_500..10_500).contains(&count), "{counts:?}");
        }
    }

    #[test]
    fn calibration_centers_on_the_decoy_median_in_mad_units() {
        // Decoys -2, 0, 2, 4: median 1; deviations 3, 1, 1, 3 -> MAD 2.
        // Targets 9 and 1 map to (9 - 1) / 2 = 4 and 0.
        let column = vec![-2.0, 0.0, 2.0, 4.0, 9.0, 1.0];
        let decoys = [true, true, true, true, false, false];
        let f = features(vec![column.clone()], 1);
        let spectra: Vec<i64> = (0..6).collect();
        let rows = Rows {
            features: &f,
            spectra: &spectra,
            folds: &spectra,
            decoys: &decoys,
        };
        let test: Vec<usize> = (0..6).collect();
        let calibrated = calibrate(&rows, &test, &column).unwrap();
        assert_eq!(calibrated, [-1.5, -0.5, 0.5, 1.5, 4.0, 0.0]);

        // Mostly tied decoys (MAD 0) fall back to the standard deviation:
        // decoys 0, 0, 0, 4: median 0, deviations 0 0 0 4 -> MAD 0; mean 1,
        // variance (1 + 1 + 1 + 9) / 4 = 3 -> spread sqrt(3).
        let tied = vec![0.0, 0.0, 0.0, 4.0, 6.0];
        let f = features(vec![tied.clone()], 1);
        let rows = Rows {
            features: &f,
            spectra: &[0, 1, 2, 3, 4],
            folds: &[0, 1, 2, 3, 4],
            decoys: &[true, true, true, true, false],
        };
        let calibrated = calibrate(&rows, &[0, 1, 2, 3, 4], &tied).unwrap();
        assert!((calibrated[4] - 6.0 / 3_f64.sqrt()).abs() < 1e-12);

        // Identical decoy scores give no scale.
        let flat = vec![1.0, 1.0, 5.0];
        let f = features(vec![flat.clone()], 1);
        let rows = Rows {
            features: &f,
            spectra: &[0, 1, 2],
            folds: &[0, 1, 2],
            decoys: &[true, true, false],
        };
        assert!(calibrate(&rows, &[0, 1, 2], &flat).is_none());
    }
}
