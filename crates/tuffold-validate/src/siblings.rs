//! Sibling peptides (iProphet's NSP), made safe for target-decoy FDR.
//!
//! A peptide of a protein whose other peptides are identified is more likely
//! right. Counted plainly, that also lifts a *wrong* peptide that happens to
//! fall in a protein that is present, and no decoy is ever lifted (decoy
//! proteins are never present), so decoys would stop modelling those errors;
//! entrapment cannot see them either (its proteins are absent). iProphet's
//! NSP model is all of its lead over Tuffold on PXD006447.
//!
//! So a decoy peptide is given the siblings of its *paired target* protein
//! (`rev_X` counts the identified peptides of `X`, the pairing picked
//! protein FDR uses). A random match lands on `X`'s peptides as often as on
//! `rev_X`'s, so wrong targets and decoys are lifted equally often, and the
//! decoys count those errors again.
//!
//! Identified means passing [`TRAIN_FDR`] by the peptide's best spectrum,
//! before pooling; siblings are distinct sequences, never the peptide's own
//! (other modified forms of it are not siblings).
//!
//! Pairing goes by name, and names fail: TPP's `reverse-5287` names no
//! target, and a Scaffold file named its targets `A0QSM5` but their decoys
//! `tr|…|…-DECOY`. Decoys then get no siblings, the plain leak again (on the
//! synthetic searches with only some of a present protein's peptides in the
//! sample: 8.8% true FDP at 1%, against 0.95% paired). So the pairing is
//! checked on the data, by what it must achieve: decoys must have siblings
//! about as often as wrong targets. Targets scoring below the decoys' median
//! are nearly all wrong (the few right ones only raise their share, making
//! the check stricter); unless decoys reach [`MIN_PAIRED_RATIO`] of their
//! share, siblings are not used. Measured: 6.4% vs 5.6% (Sage, PXD006447),
//! 42% vs 46% (Scaffold Mascot file) pass; 0% vs 11.6% (TPP numbered
//! decoys) and 0.3% vs 29.5% (that Scaffold Sequest file) do not.
//!
//! When the pairing fails, decoys' counts are imputed instead: wrong targets
//! are what decoys must stand for, and a wrong match's count depends on the
//! protein it landed on, not on its score, so the targets scoring below the
//! decoys' median give the counts wrong matches get. Each decoy, in an order
//! fixed by a hash of its id, takes the next quantile of those counts, so
//! decoys carry their distribution exactly (right matches among them only
//! lift decoys more, which is conservative). Checked with entrapment twins
//! embedded in their real proteins (four twin seeds, PXD006447 Comet):
//! twins whose real counterpart isn't identified pass at 0.74–1.00× the
//! decoys' rate (paired: 0.58–0.76×); twins of identified peptides, which
//! share their mass and take their spectra, pass more often, a class wrong
//! real peptides rarely fall in. Embedding each twin in another real
//! protein instead, so the protein whose spectra it takes isn't the one
//! whose siblings it gets, brings all twins to 0.94–1.17× (paired:
//! 0.91–1.10×). On
//! ProteomeTools pools its extra peptides were pool peptides or synthesis
//! by-products but for 3 of 16. Name pairing stays first: it is
//! conservative (reversed decoys catch their own targets' spectra and
//! inherit a present protein's siblings), which also covers wrong
//! near-copies of present peptides that imputation would not.

use std::collections::{BTreeSet, HashMap};

use tuffold_core::decoy::target_accession;

use crate::fdr::q_values;
use crate::input::Input;
use crate::pooling::PeptideSpectra;
use crate::rescore::TRAIN_FDR;

/// Decoys must have siblings at least this often relative to low-scoring
/// targets.
const MIN_PAIRED_RATIO: f64 = 0.8;

/// Why sibling counts were not used: decoys had siblings this much less
/// often than low-scoring targets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Unpaired {
    /// Share of decoy peptides with a sibling.
    pub decoys: f64,
    /// Share of targets scoring below the decoys' median with a sibling.
    pub low_targets: f64,
}

/// Sibling peptides of each peptide (parallel to `peptides`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Siblings {
    pub counts: Vec<usize>,
    /// Set when decoys were not paired with their targets by name and their
    /// counts were imputed from low-scoring targets' (module docs).
    pub imputed: Option<Unpaired>,
}

/// Sibling peptides of each peptide, decoys' from their paired targets or,
/// when names don't pair them, imputed.
pub(crate) fn sibling_counts(peptides: &[PeptideSpectra], input: &Input) -> Siblings {
    let sequence: HashMap<i64, &str> = input
        .peptides
        .iter()
        .map(|p| (p.id, p.sequence.as_str()))
        .collect();
    let mut proteins_of: HashMap<i64, Vec<i64>> = HashMap::new();
    for &(peptide, protein) in &input.peptide_proteins {
        proteins_of.entry(peptide).or_default().push(protein);
    }
    let is_decoy_protein: HashMap<i64, bool> =
        input.proteins.iter().map(|p| (p.id, p.is_decoy)).collect();
    let target_named: HashMap<&str, i64> = input
        .proteins
        .iter()
        .filter(|p| !p.is_decoy)
        .map(|p| (p.accession.as_str(), p.id))
        .collect();
    // A decoy protein's target, when its label names one the project holds.
    let pair: HashMap<i64, i64> = input
        .proteins
        .iter()
        .filter(|p| p.is_decoy)
        .filter_map(|p| {
            let target = target_accession(&p.accession)?;
            Some((p.id, *target_named.get(target)?))
        })
        .collect();

    let best: Vec<f64> = peptides.iter().map(|p| p.spectra[0].score).collect();
    let decoys: Vec<bool> = peptides.iter().map(|p| p.is_decoy).collect();
    let mut identified: HashMap<i64, BTreeSet<&str>> = HashMap::new();
    for (peptide, q) in peptides.iter().zip(q_values(&best, &decoys)) {
        if peptide.is_decoy || q > TRAIN_FDR {
            continue;
        }
        let seq = sequence
            .get(&peptide.peptide_id)
            .copied()
            .unwrap_or_default();
        for protein in proteins_of.get(&peptide.peptide_id).into_iter().flatten() {
            if is_decoy_protein.get(protein) == Some(&false) {
                identified.entry(*protein).or_default().insert(seq);
            }
        }
    }
    let siblings_in = |protein: i64, own: &str| {
        identified
            .get(&protein)
            .map_or(0, |seqs| seqs.len() - usize::from(seqs.contains(own)))
    };

    let mut counts: Vec<usize> = peptides
        .iter()
        .map(|peptide| {
            let own = sequence
                .get(&peptide.peptide_id)
                .copied()
                .unwrap_or_default();
            let proteins = proteins_of.get(&peptide.peptide_id).into_iter().flatten();
            if !peptide.is_decoy {
                return proteins
                    .filter(|protein| is_decoy_protein.get(protein) == Some(&false))
                    .map(|&protein| siblings_in(protein, own))
                    .max()
                    .unwrap_or(0);
            }
            // A target with no evidence in the project has no siblings.
            proteins
                .filter_map(|protein| pair.get(protein))
                .map(|&target| siblings_in(target, own))
                .max()
                .unwrap_or(0)
        })
        .collect();
    match check_pairing(&best, &decoys, &counts) {
        Ok(()) => Siblings {
            counts,
            imputed: None,
        },
        Err(unpaired) => {
            let ids: Vec<i64> = peptides.iter().map(|p| p.peptide_id).collect();
            impute(&best, &decoys, &ids, &mut counts);
            Siblings {
                counts,
                imputed: Some(unpaired),
            }
        }
    }
}

/// The decoys' median best score: targets at or below it are nearly all
/// wrong. `None` without decoys.
fn decoy_median(best: &[f64], decoys: &[bool]) -> Option<f64> {
    let mut scores: Vec<f64> = best
        .iter()
        .zip(decoys)
        .filter(|(_, decoy)| **decoy)
        .map(|(score, _)| *score)
        .collect();
    scores.sort_by(f64::total_cmp);
    scores.get(scores.len() / 2).copied()
}

/// Decoys take the counts of targets at or below the decoys' median: each
/// decoy, in the order of a hash of its id, the next quantile of them.
fn impute(best: &[f64], decoys: &[bool], ids: &[i64], counts: &mut [usize]) {
    let Some(median) = decoy_median(best, decoys) else {
        return;
    };
    let mut low: Vec<usize> = (0..counts.len())
        .filter(|&i| !decoys[i] && best[i] <= median)
        .map(|i| counts[i])
        .collect();
    low.sort_unstable();
    let mut order: Vec<usize> = (0..counts.len()).filter(|&i| decoys[i]).collect();
    order.sort_by_key(|&i| (mix(ids[i]), ids[i]));
    let n = order.len();
    for (k, i) in order.into_iter().enumerate() {
        counts[i] = if low.is_empty() {
            0
        } else {
            low[k * low.len() / n]
        };
    }
}

/// SplitMix64's finalizer: a fixed, well-spread order for ids.
fn mix(id: i64) -> u64 {
    let mut z = (id as u64).wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Decoys must have siblings about as often as wrong targets (module docs).
fn check_pairing(best: &[f64], decoys: &[bool], counts: &[usize]) -> Result<(), Unpaired> {
    let Some(median) = decoy_median(best, decoys) else {
        return Ok(());
    };
    let share = |keep: &dyn Fn(usize) -> bool| {
        let (with, all) = (0..counts.len())
            .filter(|&i| keep(i))
            .fold((0, 0), |(with, all), i| {
                (with + usize::from(counts[i] > 0), all + 1)
            });
        with as f64 / f64::from(u32::try_from(all).unwrap_or(u32::MAX).max(1))
    };
    let unpaired = Unpaired {
        decoys: share(&|i| decoys[i]),
        low_targets: share(&|i| !decoys[i] && best[i] <= median),
    };
    if unpaired.decoys >= MIN_PAIRED_RATIO * unpaired.low_targets {
        Ok(())
    } else {
        Err(unpaired)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::{PeptideRow, ProteinRow};
    use crate::pooling::{Support, group_by_peptide};

    fn peptide(id: i64, sequence: &str, is_decoy: bool) -> PeptideRow {
        PeptideRow {
            id,
            sequence: sequence.to_string(),
            modifications: 0,
            is_decoy,
        }
    }

    fn protein(id: i64, accession: &str, is_decoy: bool) -> ProteinRow {
        ProteinRow {
            id,
            accession: accession.to_string(),
            is_decoy,
        }
    }

    /// Protein 1 (`sp|P1|A`) has four strong peptides (10–13) and a weak one
    /// (14); its decoy `rev_sp|P1|A` (101) a weak peptide (110). Protein 2 has
    /// one weak peptide (20) and no other evidence. Protein 3 has 120 strong
    /// peptides, so that, with the +1 correction, the strong ones can pass
    /// 1% FDR (q = 1/124 at peptide 13). Decoys 120–125 are weak, from
    /// unrelated decoy proteins.
    fn project(decoy_accession: &str) -> (Vec<PeptideSpectra>, Input) {
        let mut input = Input {
            peptides: vec![
                peptide(10, "AAAK", false),
                peptide(11, "CCCK", false),
                peptide(12, "DDDK", false),
                peptide(13, "EEEK", false),
                peptide(14, "FFFK", false),
                peptide(20, "GGGK", false),
                peptide(110, "HHHK", true),
            ],
            proteins: vec![
                protein(1, "sp|P1|A", false),
                protein(2, "sp|P2|B", false),
                protein(101, decoy_accession, true),
            ],
            peptide_proteins: vec![
                (10, 1),
                (11, 1),
                (12, 1),
                (13, 1),
                (14, 1),
                (20, 2),
                (110, 101),
            ],
            ..Input::default()
        };
        let mut supports = vec![];
        for (id, score) in [(10, 9.0), (11, 8.0), (12, 7.0), (13, 6.0)] {
            supports.push((id, false, score));
        }
        // 14 and 20 fall below decoy 110: q = 2/125 and 2/126, not identified.
        supports.extend([(14, false, -0.5), (20, false, 0.4), (110, true, 0.45)]);
        input.proteins.push(protein(3, "sp|P3|C", false));
        for k in 0..120 {
            let id = 1000 + k;
            input
                .peptides
                .push(peptide(id, &format!("FILLER{k}K"), false));
            input.peptide_proteins.push((id, 3));
            supports.push((id, false, 20.0 + f64::from(i32::try_from(k).unwrap())));
        }
        for k in 0..6 {
            let id = 120 + k;
            input.peptides.push(peptide(id, "KKKK", true));
            input
                .proteins
                .push(protein(200 + k, &format!("rev_sp|Q{k}|Z"), true));
            input.peptide_proteins.push((id, 200 + k));
            supports.push((id, true, -1.0 - f64::from(i32::try_from(k).unwrap())));
        }
        let grouped = group_by_peptide(supports.into_iter().enumerate().map(
            |(spectrum, (id, decoy, score))| {
                let spectrum_id = i64::try_from(spectrum).unwrap();
                let support = Support {
                    spectrum_id,
                    charge: 2,
                    score,
                };
                (id, decoy, support)
            },
        ));
        (grouped, input)
    }

    #[test]
    fn a_decoy_counts_the_identified_peptides_of_its_target() {
        let (grouped, input) = project("rev_sp|P1|A");
        let siblings = sibling_counts(&grouped, &input);
        assert_eq!(siblings.imputed, None, "paired by name");
        let counts = siblings.counts;
        let by_id: HashMap<i64, usize> = grouped.iter().map(|p| p.peptide_id).zip(counts).collect();
        // 10–13 pass 1%: each has the three others as siblings.
        assert_eq!(by_id[&10], 3);
        // The weak 14 of the same protein has all four.
        assert_eq!(by_id[&14], 4);
        // Protein 2 has nothing identified.
        assert_eq!(by_id[&20], 0);
        // The decoy of protein 1 gets protein 1's four, as a wrong target
        // landing there would.
        assert_eq!(by_id[&110], 4);
        // Decoys whose targets have no evidence get none.
        assert_eq!(by_id[&120], 0);
        // Protein 3's peptides have its 119 others.
        assert_eq!(by_id[&1000], 119);
    }

    #[test]
    fn a_decoy_naming_no_target_in_the_project_gets_no_siblings() {
        // TPP's numbered `reverse-1` names no protein the project holds.
        let (grouped, input) = project("reverse-1");
        let siblings = sibling_counts(&grouped, &input);
        assert_eq!(siblings.imputed, None, "too few to judge");
        let counts = siblings.counts;
        let decoy = grouped.iter().position(|p| p.peptide_id == 110).unwrap();
        assert_eq!(counts[decoy], 0);
    }

    #[test]
    fn decoys_must_have_siblings_about_as_often_as_low_scoring_targets() {
        // Decoys score 0..=3 (median 2). Targets at 0.5 and 1.5 are below it,
        // one with siblings: 50%. The target at 5 is above and not counted.
        let best = [0.0, 1.0, 2.0, 3.0, 0.5, 1.5, 5.0];
        let decoys = [true, true, true, true, false, false, false];
        let check = |decoy_counts: [usize; 4]| {
            let mut counts = decoy_counts.to_vec();
            counts.extend([1, 0, 3]);
            check_pairing(&best, &decoys, &counts)
        };
        assert_eq!(
            check([0, 0, 0, 0]),
            Err(Unpaired {
                decoys: 0.0,
                low_targets: 0.5
            })
        );
        // 25% is under 0.8 × 50%; 50% is not.
        assert!(check([1, 0, 0, 0]).is_err());
        assert_eq!(check([1, 1, 0, 0]), Ok(()));
    }

    #[test]
    fn unpaired_decoys_take_the_low_scoring_targets_counts() {
        // Decoys score 0..=3 (median 2); targets at 0.5, 1.0, 1.5 and 2.0
        // are at or below it, with 4, 0, 2 and 0 siblings; the target at 5
        // (7 siblings) is above and doesn't count.
        let best = [0.0, 1.0, 2.0, 3.0, 0.5, 1.0, 1.5, 2.0, 5.0];
        let decoys = [true, true, true, true, false, false, false, false, false];
        let ids = [11, 12, 13, 14, 21, 22, 23, 24, 25];
        let mut counts = vec![0, 0, 0, 0, 4, 0, 2, 0, 7];
        impute(&best, &decoys, &ids, &mut counts);
        let mut taken = counts[..4].to_vec();
        taken.sort_unstable();
        assert_eq!(taken, [0, 0, 2, 4], "each low target's count once");
        assert_eq!(counts[4..], [4, 0, 2, 0, 7], "targets keep theirs");
        // The same ids always get the same counts.
        let mut again = vec![0, 0, 0, 0, 4, 0, 2, 0, 7];
        impute(&best, &decoys, &ids, &mut again);
        assert_eq!(again, counts);
    }

    #[test]
    fn decoys_naming_no_target_are_imputed_when_low_targets_have_siblings() {
        // Weak targets of protein 1 sit among the decoys, so low-scoring
        // targets have siblings and decoys naming nothing would have none.
        let (mut grouped, input) = project("reverse-1");
        for peptide in &mut grouped {
            if peptide.peptide_id == 14 {
                peptide.spectra[0].score = -5.5;
            }
        }
        let siblings = sibling_counts(&grouped, &input);
        let unpaired = siblings.imputed.expect("imputed");
        assert_eq!(unpaired.decoys, 0.0);
        assert!(unpaired.low_targets > 0.0);
        // Peptide 14 (4 siblings) is the only low target: every decoy takes 4.
        for (peptide, count) in grouped.iter().zip(&siblings.counts) {
            if peptide.is_decoy {
                assert_eq!(*count, 4, "decoy {}", peptide.peptide_id);
            }
        }
    }
}
