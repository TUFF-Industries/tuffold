//! Replicate spectra (iProphet's NRS), made safe for target-decoy FDR.
//!
//! A spectrum matched to a peptide that other spectra confidently identify is
//! more likely right. Counted plainly that lifts a *wrong* spectrum landing
//! on a present peptide too, and decoy peptides are never present, so decoys
//! stop modelling those errors (synthetic searches: true FDP 1.18% → 3.29% at
//! 1%). Entrapment can't see it either.
//!
//! So, as for sibling peptides, a decoy PSM also gets the count of its
//! *paired target peptide*, the target its sequence was made from. A wrong
//! match lands on that target as often as on the decoy (same mass, drawn the
//! same way), so wrong targets and decoys are lifted equally often and the
//! decoys count those errors again.
//!
//! Wrong matches also repeat: a wrong peptide resembling an abundant one can
//! catch dozens of its spectra, and once a few pass, the rest are lifted
//! (PXD006447 with entrapment: one twin peptide gained 41 spectra, spectrum
//! entrapment 1.02% → 1.75%). A decoy can repeat the same way, so replicates
//! are counted label-blind: a peptide's other spectra scoring at or above
//! the first pass's [`REPLICATE_FDR`] cut, decoy or not. A decoy gets the
//! larger of its own count and its paired target's.
//!
//! Pairing is by sequence, for the common decoy kinds: pseudo-reversed
//! (reversed but the C-terminal residue, as Sage makes them), reversed whole,
//! and peptides of a reversed protein (the reversed sequence shifted by one:
//! `FEDK` in `…RFEDKZ…` pairs with `DEFR`). Other decoys (shuffled, or
//! names only) don't pair, so the pairing is checked on the data as
//! siblings' is ([`MIN_PAIRED_RATIO`]). Failing it, decoys' counts are
//! imputed as siblings' are (`siblings.rs`): each decoy peptide, in the
//! order of a hash of its id, takes the next quantile of the counts of
//! targets scoring below the decoys' median (nearly all wrong; a wrong
//! match's count depends on the peptide it landed on, not its score), and
//! each of its PSMs the larger of that and its own count. One draw per
//! decoy peptide, as a wrong target's spectra share the peptide they landed
//! on.
//!
//! The cut ([`REPLICATE_FDR`]) comes from the rescoring before this
//! feature; a PSM's own spectrum is never counted.

use std::collections::{HashMap, HashSet};

use crate::input::PeptideRow;

pub const REPLICATES: &str = "log(1 + other confident spectra of the peptide)";
/// Its name when decoys' counts were imputed (module docs).
pub const REPLICATES_IMPUTED: &str =
    "log(1 + other confident spectra of the peptide; decoys' imputed)";

/// Replicates are spectra passing this FDR in the first pass. At 1% (the
/// rescoring's training cut) wrong peptides lifted each other: on Sage with
/// entrapment twins, a twin near the cut had replicates 20% of the time
/// against a decoy's 17%, and the twins' estimate went 1.11% → 1.23%; at
/// 0.5% the two match (8.4% vs 8.1%) and nothing is lost.
const REPLICATE_FDR: f64 = 0.005;

/// Decoys must have replicates at least this often relative to low-scoring
/// targets. Lower than siblings' 0.8: the right spectra among those targets
/// always have replicates (their peptide is present), which inflates the
/// target share (synthetic searches: 8.0% vs 10.4% with pairing intact),
/// while a failed pairing leaves decoys near 0%.
const MIN_PAIRED_RATIO: f64 = 0.5;

/// Why replicate counts were not used.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Unpaired {
    /// Share of decoy PSMs whose paired target has confident spectra.
    pub decoys: f64,
    /// Share of targets scoring below the decoys' median with replicates.
    pub low_targets: f64,
}

/// The replicate feature, `log(1 + replicates)` per row.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Replicates {
    pub column: Vec<f64>,
    /// Set when decoys were not paired and their counts were imputed.
    pub imputed: Option<Unpaired>,
}

/// `log(1 + replicates)` per row, decoys' from their paired targets or,
/// when they don't pair, imputed.
pub(crate) fn replicate_feature(
    peptides: &[&PeptideRow],
    spectra: &[i64],
    decoys: &[bool],
    scores: &[f64],
    q: &[Option<f64>],
) -> Replicates {
    // The score that passes REPLICATE_FDR, applied to decoys and targets alike.
    let cut = (0..peptides.len())
        .filter(|&row| !decoys[row] && q[row].is_some_and(|q| q <= REPLICATE_FDR))
        .map(|row| scores[row])
        .min_by(f64::total_cmp)
        .unwrap_or(f64::INFINITY);
    let mut confident: HashMap<i64, HashSet<i64>> = HashMap::new();
    for row in 0..peptides.len() {
        if q[row].is_some() && scores[row] >= cut {
            confident
                .entry(peptides[row].id)
                .or_default()
                .insert(spectra[row]);
        }
    }
    let targets: HashMap<&str, i64> = peptides
        .iter()
        .filter(|p| !p.is_decoy)
        .map(|p| (p.sequence.as_str(), p.id))
        .collect();
    // Targets by all but their last residue, for reversed-protein decoys.
    let mut by_stem: HashMap<&str, Vec<i64>> = HashMap::new();
    for peptide in peptides.iter().filter(|p| !p.is_decoy) {
        if let Some(stem) = peptide
            .sequence
            .get(..peptide.sequence.len().saturating_sub(1))
        {
            let ids = by_stem.entry(stem).or_default();
            if !ids.contains(&peptide.id) {
                ids.push(peptide.id);
            }
        }
    }
    let replicates = |peptide: i64, own: i64| {
        confident
            .get(&peptide)
            .map_or(0, |seen| seen.len() - usize::from(seen.contains(&own)))
    };
    let mut counts: Vec<usize> = (0..peptides.len())
        .map(|row| {
            let own = replicates(peptides[row].id, spectra[row]);
            if !decoys[row] {
                return own;
            }
            paired_targets(&peptides[row].sequence, &targets, &by_stem)
                .into_iter()
                .map(|target| replicates(target, spectra[row]))
                .fold(own, usize::max)
        })
        .collect();
    let imputed = check_pairing(scores, decoys, &counts).err();
    if imputed.is_some() {
        impute(peptides, decoys, scores, &mut counts);
    }
    Replicates {
        column: counts.iter().map(|&n| (n as f64).ln_1p()).collect(),
        imputed,
    }
}

/// Each decoy peptide takes a quantile of the counts of targets scoring
/// below the decoys' median (module docs); its PSMs keep their own count
/// when larger.
fn impute(peptides: &[&PeptideRow], decoys: &[bool], scores: &[f64], counts: &mut [usize]) {
    let Some(median) = decoy_median(scores, decoys) else {
        return;
    };
    let mut low: Vec<usize> = (0..counts.len())
        .filter(|&i| !decoys[i] && scores[i] < median)
        .map(|i| counts[i])
        .collect();
    low.sort_unstable();
    let mut ids: Vec<i64> = (0..counts.len())
        .filter(|&i| decoys[i])
        .map(|i| peptides[i].id)
        .collect();
    ids.sort_unstable();
    ids.dedup();
    ids.sort_by_key(|&id| (mix(id), id));
    let draw: HashMap<i64, usize> = ids
        .iter()
        .enumerate()
        .map(|(k, &id)| {
            let value = if low.is_empty() {
                0
            } else {
                low[k * low.len() / ids.len()]
            };
            (id, value)
        })
        .collect();
    for i in (0..counts.len()).filter(|&i| decoys[i]) {
        counts[i] = counts[i].max(draw[&peptides[i].id]);
    }
}

/// SplitMix64's finalizer: a fixed, well-spread order for ids.
fn mix(id: i64) -> u64 {
    let mut z = (id as u64).wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// The decoys' median score; `None` without decoys.
fn decoy_median(scores: &[f64], decoys: &[bool]) -> Option<f64> {
    let mut decoy_scores: Vec<f64> = (0..scores.len())
        .filter(|&i| decoys[i])
        .map(|i| scores[i])
        .collect();
    decoy_scores.sort_by(f64::total_cmp);
    decoy_scores.get(decoy_scores.len() / 2).copied()
}

/// The target peptides a decoy sequence may have been made from.
fn paired_targets(
    decoy: &str,
    targets: &HashMap<&str, i64>,
    by_stem: &HashMap<&str, Vec<i64>>,
) -> Vec<i64> {
    let reversed: String = decoy.chars().rev().collect();
    let mut found = Vec::new();
    if let Some((body, last)) = decoy.char_indices().last().map(|(i, c)| (&decoy[..i], c)) {
        let pseudo: String = body.chars().rev().chain(std::iter::once(last)).collect();
        found.extend(targets.get(pseudo.as_str()));
    }
    found.extend(targets.get(reversed.as_str()));
    // A reversed protein's peptide starts after its K/R: drop that, and the
    // target ends with the cleavage residue that preceded the decoy.
    if let Some(stem) = reversed.strip_prefix(['K', 'R']) {
        found.extend(by_stem.get(stem).into_iter().flatten());
    }
    found.sort_unstable();
    found.dedup();
    found
}

/// Decoys must have replicates about as often as wrong targets: targets
/// scoring below the decoys' median are nearly all wrong.
fn check_pairing(scores: &[f64], decoys: &[bool], counts: &[usize]) -> Result<(), Unpaired> {
    let Some(median) = decoy_median(scores, decoys) else {
        return Ok(());
    };
    let share = |rows: Vec<usize>| {
        let lifted = rows.iter().filter(|&&i| counts[i] > 0).count();
        lifted as f64 / rows.len().max(1) as f64
    };
    let decoy_share = share((0..scores.len()).filter(|&i| decoys[i]).collect());
    let low_target_share = share(
        (0..scores.len())
            .filter(|&i| !decoys[i] && scores[i] < median)
            .collect(),
    );
    if decoy_share >= MIN_PAIRED_RATIO * low_target_share {
        Ok(())
    } else {
        Err(Unpaired {
            decoys: decoy_share,
            low_targets: low_target_share,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peptide(id: i64, sequence: &str, is_decoy: bool) -> PeptideRow {
        PeptideRow {
            id,
            sequence: sequence.into(),
            modifications: 0,
            is_decoy,
        }
    }

    #[test]
    fn decoys_pair_with_the_target_they_were_made_from() {
        let rows = [
            peptide(1, "PEPTIDEK", false),
            peptide(2, "DEFR", false),
            peptide(3, "ACDK", false),
        ];
        let targets: HashMap<&str, i64> =
            rows.iter().map(|p| (p.sequence.as_str(), p.id)).collect();
        let mut by_stem: HashMap<&str, Vec<i64>> = HashMap::new();
        for p in &rows {
            by_stem
                .entry(&p.sequence[..p.sequence.len() - 1])
                .or_default()
                .push(p.id);
        }
        // Pseudo-reversed, wholly reversed, and from a reversed protein.
        assert_eq!(paired_targets("EDITPEPK", &targets, &by_stem), [1]);
        assert_eq!(paired_targets("KEDITPEP", &targets, &by_stem), [1]);
        assert_eq!(paired_targets("FEDK", &targets, &by_stem), [2]);
        assert!(paired_targets("WWWWK", &targets, &by_stem).is_empty());
    }

    #[test]
    fn a_decoy_counts_its_paired_targets_confident_spectra_never_its_own() {
        let target = peptide(1, "PEPTIDEK", false);
        let decoy = peptide(2, "EDITPEPK", true);
        let peptides = [&target, &target, &target, &decoy];
        let spectra = [10, 11, 12, 13];
        let decoys = [false, false, false, true];
        let scores = [5.0, 4.0, 1.0, 0.5];
        let q = [Some(0.0), Some(0.0), Some(0.5), None];
        let replicates = replicate_feature(&peptides, &spectra, &decoys, &scores, &q);
        assert_eq!(replicates.imputed, None);
        let feature = replicates.column;
        // Spectra 10 and 11 are confident: each has one other; 12 has two;
        // the decoy gets its target's two.
        assert_eq!(
            feature,
            [1f64.ln_1p(), 1f64.ln_1p(), 2f64.ln_1p(), 2f64.ln_1p()]
        );
    }

    #[test]
    fn decoys_repeating_on_their_own_count_like_targets() {
        let decoy = peptide(2, "EDITPEPK", true);
        let other = peptide(3, "ACDEFGHK", false);
        // Two decoy spectra score above the cut set by the confident target.
        let peptides = [&other, &decoy, &decoy];
        let spectra = [1, 2, 3];
        let decoys = [false, true, true];
        let scores = [3.0, 4.0, 5.0];
        let q = [Some(0.0), Some(0.5), Some(0.5)];
        let replicates = replicate_feature(&peptides, &spectra, &decoys, &scores, &q);
        assert_eq!(replicates.imputed, None);
        let feature = replicates.column;
        assert_eq!(feature, [0.0, 1f64.ln_1p(), 1f64.ln_1p()]);
    }

    #[test]
    fn unpaired_decoys_take_low_scoring_targets_counts() {
        let present = peptide(1, "PEPTIDEK", false);
        let shuffled = peptide(2, "TPEDIPEK", true);
        // Low-scoring targets land on the present peptide; the decoys pair
        // with nothing.
        let peptides = [&present, &present, &present, &shuffled, &shuffled];
        let spectra = [1, 2, 3, 4, 5];
        let decoys = [false, false, false, true, true];
        let scores = [9.0, 0.1, 0.2, 1.0, 1.0];
        let q = [Some(0.0), Some(0.5), Some(0.5), None, None];
        let replicates = replicate_feature(&peptides, &spectra, &decoys, &scores, &q);
        let unpaired = replicates.imputed.expect("imputed");
        assert_eq!(unpaired.decoys, 0.0);
        assert_eq!(unpaired.low_targets, 1.0);
        // Both low targets (0.1 and 0.2, below the decoys' median 1.0) have
        // spectrum 1 as their replicate: the decoy peptide takes 1, on both
        // its PSMs.
        let one = 1f64.ln_1p();
        assert_eq!(replicates.column, [0.0, one, one, one, one]);
    }
}
