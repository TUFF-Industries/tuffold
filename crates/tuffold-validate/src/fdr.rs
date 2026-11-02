//! Target-decoy false discovery rates and posterior error probabilities.
//!
//! q-values follow target-decoy competition (Elias & Gygi 2007) with the +1
//! correction (Levitsky et al. 2017), matching mokapot's `tdc`:
//! `FDR(t) = (decoys scoring >= t + 1) / targets scoring >= t`, evaluated at
//! the end of each group of tied scores, made monotone from the worst score
//! up, and capped at 1. Callers compete first so each spectrum, peptide, or
//! protein pair contributes one entry ([`winners`]).
//!
//! PEPs come from isotonic regression of the decoy indicator on the score
//! (pool-adjacent-violators). Under competition, decoys stand in for
//! incorrect targets, so a fitted decoy fraction `p` gives
//! `PEP = p / (1 - p)`, capped at 1.
//!
//! Higher scores are better everywhere; scores must be finite, and parallel
//! slices must have equal lengths (the functions panic otherwise).

/// q-value of every entry. `scores` and `is_decoy` are parallel.
pub fn q_values(scores: &[f64], is_decoy: &[bool]) -> Vec<f64> {
    assert_eq!(
        scores.len(),
        is_decoy.len(),
        "scores and labels differ in length"
    );
    let order = descending(scores);
    let mut q = vec![1.0; scores.len()];
    let mut ties = Vec::new();
    let (mut targets, mut decoys) = (0_usize, 0_usize);
    for group in tie_groups(&order, scores) {
        for &index in group {
            if is_decoy[index] {
                decoys += 1;
            } else {
                targets += 1;
            }
        }
        let fdr = if targets == 0 {
            1.0
        } else {
            (decoys + 1) as f64 / targets as f64
        };
        ties.push((group, fdr));
    }
    let mut lowest = 1.0_f64;
    for (group, fdr) in ties.iter().rev() {
        lowest = lowest.min(*fdr);
        for &index in *group {
            q[index] = lowest;
        }
    }
    q
}

/// Posterior error probability of every entry. `scores` and `is_decoy` are
/// parallel.
pub fn peps(scores: &[f64], is_decoy: &[bool]) -> Vec<f64> {
    assert_eq!(
        scores.len(),
        is_decoy.len(),
        "scores and labels differ in length"
    );
    let order = descending(scores);
    // Blocks of (decoys, entries, tie groups covered), merged until the decoy
    // fraction never decreases from the best score to the worst.
    let mut blocks: Vec<(f64, f64, Vec<&[usize]>)> = Vec::new();
    for group in tie_groups(&order, scores) {
        let decoys = group.iter().filter(|&&index| is_decoy[index]).count() as f64;
        blocks.push((decoys, group.len() as f64, vec![group]));
        while blocks.len() > 1 {
            let last = &blocks[blocks.len() - 1];
            let before = &blocks[blocks.len() - 2];
            if before.0 / before.1 <= last.0 / last.1 {
                break;
            }
            let (decoys, entries, groups) = blocks.pop().expect("two blocks");
            let merged = blocks.last_mut().expect("two blocks");
            merged.0 += decoys;
            merged.1 += entries;
            merged.2.extend(groups);
        }
    }
    let mut pep = vec![1.0; scores.len()];
    for (decoys, entries, groups) in blocks {
        let fraction = decoys / entries;
        let value = if fraction >= 1.0 {
            1.0
        } else {
            (fraction / (1.0 - fraction)).min(1.0)
        };
        for index in groups.into_iter().flatten() {
            pep[*index] = value;
        }
    }
    pep
}

/// PEPs for `query` scores read off the isotonic fit of a separate training
/// set: each query takes the PEP of the best training entry scoring at or
/// below it (the top entry's when it outscores them all). Fitting on other
/// data keeps an entry's own label out of its PEP. With no training data
/// every PEP is 1 (no evidence).
pub fn peps_from(train_scores: &[f64], train_decoys: &[bool], query: &[f64]) -> Vec<f64> {
    if train_scores.is_empty() {
        return vec![1.0; query.len()];
    }
    let fitted = peps(train_scores, train_decoys);
    let order = descending(train_scores);
    let sorted: Vec<f64> = order.iter().map(|&i| train_scores[i]).collect();
    query
        .iter()
        .map(|score| {
            let above = sorted.partition_point(|train| train > score);
            fitted[order[above.min(order.len() - 1)]]
        })
        .collect()
}

/// Indices of the entries that survive competition: the best score per key.
/// A tie goes to the decoy (the conservative choice), then to the earlier
/// entry. Returned in input order.
pub fn winners(keys: &[i64], scores: &[f64], is_decoy: &[bool]) -> Vec<usize> {
    use std::collections::HashMap;
    assert!(
        keys.len() == scores.len() && scores.len() == is_decoy.len(),
        "keys, scores, and labels differ in length"
    );
    let mut best: HashMap<i64, usize> = HashMap::new();
    for (index, key) in keys.iter().enumerate() {
        best.entry(*key)
            .and_modify(|current| {
                let better = scores[index] > scores[*current]
                    || (scores[index] == scores[*current]
                        && is_decoy[index]
                        && !is_decoy[*current]);
                if better {
                    *current = index;
                }
            })
            .or_insert(index);
    }
    let mut kept: Vec<usize> = best.into_values().collect();
    kept.sort_unstable();
    kept
}

/// Indices sorted by score, best first (stable for equal scores).
fn descending(scores: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    order
}

/// Consecutive runs of equal scores in `order`.
fn tie_groups<'a>(order: &'a [usize], scores: &'a [f64]) -> impl Iterator<Item = &'a [usize]> {
    order.chunk_by(move |&a, &b| scores[a] == scores[b])
}

/// `(q-value, PEP)` of every entry of a competed level, for reporting.
///
/// The PEP is floored at the q-value. An FDR is the average PEP of what
/// passes, so with PEPs rising down the list no entry's PEP can be below its
/// q-value; the isotonic fit alone can say 0 (no decoys yet), which would
/// read as certainty. The floor carries the q-value's +1 correction — the
/// small-sample caution — into probabilities.
pub fn confidence(scores: &[f64], is_decoy: &[bool]) -> Vec<(f64, f64)> {
    q_values(scores, is_decoy)
        .into_iter()
        .zip(peps(scores, is_decoy))
        .map(|(q, pep)| (q, pep.max(q)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() < 1e-12, "{actual:?} != {expected:?}");
        }
    }

    // The q-value expectations below were derived by hand; mokapot 0.10.0's
    // `qvalues.tdc` returns the same values for every one of these vectors.

    #[test]
    fn reported_peps_never_fall_below_the_q_value() {
        // Five targets and no decoy: the isotonic PEP is 0 for all, but the
        // +1-corrected q-value is 1/5 throughout, so no 100% probabilities.
        let found = confidence(&[5.0, 4.0, 3.0, 2.0, 1.0], &[false; 5]);
        assert!(
            found.iter().all(|&(q, pep)| q == 0.2 && pep == 0.2),
            "{found:?}"
        );
        // Above the floor the isotonic PEP stands: T T D D gives PEPs
        // 0, 0, 1, 1 and q-values 1/2, 1/2, 1, 1.
        let found = confidence(&[4.0, 3.0, 2.0, 1.0], &[false, false, true, true]);
        assert_eq!(found, [(0.5, 0.5), (0.5, 0.5), (1.0, 1.0), (1.0, 1.0)]);
    }

    #[test]
    fn q_values_follow_target_decoy_competition_with_the_plus_one_correction() {
        // Best first: T T D T D T.
        // FDR after each: 1/1, 1/2, 2/2, 2/3, 3/3, 3/4.
        // Monotone from the bottom: 3/4, 3/4, 2/3, 2/3, 1/2, 1/2.
        let scores = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0];
        let decoy = [false, false, true, false, true, false];
        close(
            &q_values(&scores, &decoy),
            &[0.5, 0.5, 2.0 / 3.0, 2.0 / 3.0, 0.75, 0.75],
        );
    }

    #[test]
    fn q_values_do_not_depend_on_input_order() {
        let scores = [6.0, 10.0, 8.0, 5.0, 9.0, 7.0];
        let decoy = [true, false, true, false, false, false];
        close(
            &q_values(&scores, &decoy),
            &[0.75, 0.5, 2.0 / 3.0, 0.75, 0.5, 2.0 / 3.0],
        );
    }

    #[test]
    fn tied_scores_share_the_fdr_at_the_end_of_their_group() {
        // Groups: {3, 3} = T D -> (1+1)/1 = 2; {2} T -> 2/2 = 1; {1} T -> 2/3.
        // A target tied with a decoy gets no credit for sorting first.
        let scores = [3.0, 3.0, 2.0, 1.0];
        let decoy = [false, true, false, false];
        close(&q_values(&scores, &decoy), &[2.0 / 3.0; 4]);
        // Group {3, 3} = T T -> 1/2; {2} D -> 2/2; {1} D -> 3/2, capped at 1.
        let scores = [3.0, 3.0, 2.0, 1.0];
        let decoy = [false, false, true, true];
        close(&q_values(&scores, &decoy), &[0.5, 0.5, 1.0, 1.0]);
    }

    #[test]
    fn q_values_are_one_until_a_target_is_seen() {
        close(&q_values(&[2.0, 1.0], &[true, false]), &[1.0, 1.0]);
        assert!(q_values(&[], &[]).is_empty());
    }

    #[test]
    fn peps_pool_adjacent_violators() {
        // Decoy indicator best first: 0 0 0 1 0 0 1 1. The 1 0 0 run violates
        // monotonicity and pools to 1/3: PEP = (1/3) / (2/3) = 0.5.
        let scores = [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let decoy = [false, false, false, true, false, false, true, true];
        close(
            &peps(&scores, &decoy),
            &[0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0],
        );
    }

    #[test]
    fn peps_treat_tied_scores_as_one_block() {
        // {2, 2} = T D (decoy fraction 1/2), then {1} = T (0): decreasing, so
        // they pool to 1/3 and all three get PEP (1/3) / (2/3) = 0.5.
        let scores = [2.0, 2.0, 1.0];
        let decoy = [false, true, false];
        close(&peps(&scores, &decoy), &[0.5, 0.5, 0.5]);
    }

    #[test]
    fn held_out_peps_read_the_training_curve() {
        // Training curve (see peps_pool_adjacent_violators): scores 8..=1 have
        // PEPs 0 0 0 0.5 0.5 0.5 1 1.
        let train = [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let decoy = [false, false, false, true, false, false, true, true];
        // 9 is above all -> score 8's 0; 5 -> itself, 0.5; 4.5 -> score 4,
        // 0.5; 2.5 -> score 2, 1; 0 is below all -> score 1's 1.
        let query = [9.0, 5.0, 4.5, 2.5, 0.0];
        close(
            &peps_from(&train, &decoy, &query),
            &[0.0, 0.5, 0.5, 1.0, 1.0],
        );
        assert_eq!(peps_from(&[], &[], &[3.0]), [1.0]);
    }

    #[test]
    fn competition_keeps_the_best_entry_per_key_and_gives_ties_to_decoys() {
        let keys = [7, 7, 8, 8, 9];
        let scores = [1.0, 2.0, 5.0, 5.0, 0.0];
        let decoy = [false, false, false, true, false];
        assert_eq!(winners(&keys, &scores, &decoy), vec![1, 3, 4]);
    }
}
