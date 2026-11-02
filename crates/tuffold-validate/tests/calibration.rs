//! FDR control on average: across many seeded synthetic searches, the true
//! false discovery proportion among accepted targets must not exceed the
//! nominal 1% (beyond sampling noise). One run can land above 1% by chance;
//! a systematic bias shows up only in the mean.

mod support;

use std::collections::HashSet;

use support::{
    generate_partial, generate_with, only_engine, with_retention_times, with_unnamed_decoys,
    with_unpaired_decoy_sequences,
};
use tuffold_validate::{Settings, validate};

const SEEDS: u64 = 40;
const SPECTRA: u64 = 3000;
const FDR: f64 = 0.01;

/// Accepted at 1% FDR (q-value) or at 99% probability (PEP ≤ 1%).
#[derive(Clone, Copy)]
enum Cut {
    Fdr,
    Probability,
}

impl Cut {
    fn accepts(self, q_value: f64, pep: f64) -> bool {
        match self {
            Self::Fdr => q_value <= FDR,
            Self::Probability => pep <= FDR,
        }
    }
}

/// Mean and standard error of the true FDP at spectrum and peptide level.
fn measure(engines: usize, near_miss: f64) -> [(f64, f64); 2] {
    measure_cut(engines, near_miss, Cut::Fdr)
}

fn measure_cut(engines: usize, near_miss: f64, cut: Cut) -> [(f64, f64); 2] {
    measure_partial(engines, near_miss, cut, support::PEPTIDES_PER_PROTEIN)
}

fn measure_partial(
    engines: usize,
    near_miss: f64,
    cut: Cut,
    present_peptides: u64,
) -> [(f64, f64); 2] {
    measure_searches(engines, cut, |seed| {
        generate_partial(seed, SPECTRA, engines, near_miss, present_peptides)
    })
}

fn measure_searches(
    engines: usize,
    cut: Cut,
    search: impl Fn(u64) -> support::Synthetic,
) -> [(f64, f64); 2] {
    let mut spectra = Vec::new();
    let mut peptides = Vec::new();
    for seed in 1..=SEEDS {
        let mut data = search(seed);
        if engines == 1 {
            data = only_engine(data, 0);
        }
        let validation = validate(&data.input, &Settings::default()).unwrap();
        let accepted: Vec<_> = validation
            .spectra
            .iter()
            .filter(|m| !m.is_decoy && cut.accepts(m.q_value, m.pep))
            .collect();
        let wrong = accepted
            .iter()
            .filter(|m| {
                data.true_peptides[usize::try_from(m.spectrum_id - 1).unwrap()] != m.peptide_id
            })
            .count();
        spectra.push(wrong as f64 / accepted.len().max(1) as f64);
        let truth: HashSet<i64> = data.true_peptides.iter().copied().collect();
        let accepted: Vec<_> = validation
            .peptides
            .iter()
            .filter(|p| !p.is_decoy && cut.accepts(p.q_value, p.pep))
            .collect();
        let wrong = accepted
            .iter()
            .filter(|p| !truth.contains(&p.peptide_id))
            .count();
        peptides.push(wrong as f64 / accepted.len().max(1) as f64);
    }
    [spectra, peptides].map(|values| {
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
        (mean, (variance / n).sqrt())
    })
}

fn assert_controlled(label: &str, engines: usize, near_miss: f64) {
    assert_cut_controlled(label, engines, near_miss, Cut::Fdr);
}

fn assert_cut_controlled(label: &str, engines: usize, near_miss: f64, cut: Cut) {
    for (level, (mean, error)) in ["spectra", "peptides"]
        .iter()
        .zip(measure_cut(engines, near_miss, cut))
    {
        eprintln!(
            "{label} {level}: mean FDP {:.3}% ± {:.3}%",
            mean * 100.0,
            error * 100.0
        );
        assert!(
            mean <= FDR + 2.0 * error,
            "{label} {level}: mean FDP {mean:.4} above {FDR} (± {error:.4})"
        );
    }
}

#[test]
fn one_engine_controls_the_fdr_on_average() {
    assert_controlled("one engine", 1, 0.0);
}

#[test]
fn two_independent_engines_control_the_fdr_on_average() {
    assert_controlled("two engines", 2, 0.0);
}

/// A 99% probability threshold (PEP ≤ 1%) must not let more than 1% wrong
/// identifications through on average either: what passes it is no worse
/// than what passes 1% FDR (reported PEPs never fall below the q-value).
#[test]
fn a_99_percent_probability_threshold_is_at_least_as_strict() {
    assert_cut_controlled("one engine, 99% probability", 1, 0.0, Cut::Probability);
    assert_cut_controlled("two engines, 99% probability", 2, 0.0, Cut::Probability);
}

/// Engines agreeing on the same wrong target (a near-miss) is an error decoys
/// cannot model: with 5% shared near-misses one engine alone already exceeds
/// 1% — a limit of target-decoy validation itself. Taking the strongest vote
/// lets two engines' evidence for one spectrum count once, which raises that
/// excess only slightly (200 seeds: 2.11% -> 2.32%; summing votes gave 2.7%).
/// The bound catches a regression towards summing.
#[test]
fn combining_engines_barely_amplifies_shared_near_misses() {
    let one = measure(1, 0.05);
    let two = measure(2, 0.05);
    for (level, ((one_mean, one_error), (two_mean, two_error))) in
        ["spectra", "peptides"].iter().zip(one.into_iter().zip(two))
    {
        eprintln!(
            "near-miss {level}: one engine {:.3}%, two engines {:.3}%",
            one_mean * 100.0,
            two_mean * 100.0
        );
        let noise = 2.0 * (one_error.powi(2) + two_error.powi(2)).sqrt();
        assert!(
            two_mean <= 1.15 * one_mean + noise,
            "{level}: two engines {two_mean:.4} vs one {one_mean:.4}"
        );
    }
}

/// Sibling peptides lift every peptide of a protein with other identified
/// peptides, wrong ones too. With only 2 of each present protein's 5
/// peptides in the sample, wrong matches often land on the other 3, are
/// lifted, and only decoys given their paired target's siblings can count
/// them (`siblings.rs`). Counted plainly, decoys were never lifted.
#[test]
fn sibling_peptides_do_not_lift_wrong_peptides_past_the_decoys() {
    for engines in [1, 2] {
        for (level, (mean, error)) in
            ["spectra", "peptides"]
                .iter()
                .zip(measure_partial(engines, 0.0, Cut::Fdr, 2))
        {
            eprintln!(
                "{engines} engine(s), 2 of 5 peptides present, {level}: mean FDP {:.3}% ± {:.3}%",
                mean * 100.0,
                error * 100.0
            );
            assert!(
                mean <= FDR + 2.0 * error,
                "{level}: mean FDP {mean:.4} above {FDR} (± {error:.4})"
            );
        }
    }
}

/// The same searches with decoy proteins that name no target: name pairing
/// fails, and decoys' sibling counts are imputed from low-scoring targets'
/// (`siblings.rs`). Before, siblings were left out there.
#[test]
fn imputed_sibling_counts_keep_the_fdr_controlled() {
    for engines in [1, 2] {
        let measured = measure_searches(engines, Cut::Fdr, |seed| {
            with_unnamed_decoys(generate_partial(seed, SPECTRA, engines, 0.0, 2))
        });
        for (level, (mean, error)) in ["spectra", "peptides"].iter().zip(measured) {
            eprintln!(
                "{engines} engine(s), 2 of 5 peptides present, decoys unnamed, {level}: \
                 mean FDP {:.3}% ± {:.3}%",
                mean * 100.0,
                error * 100.0
            );
            assert!(
                mean <= FDR + 2.0 * error,
                "{level}: mean FDP {mean:.4} above {FDR} (± {error:.4})"
            );
        }
    }
}

/// Decoys whose sequences name no target (rotated, as shuffled decoys) get
/// replicate counts drawn from low-scoring targets' (`replicates.rs`).
/// Before, replicates were left out there.
#[test]
fn imputed_replicate_counts_keep_the_fdr_controlled() {
    for engines in [1, 2] {
        let measured = measure_searches(engines, Cut::Fdr, |seed| {
            with_unpaired_decoy_sequences(generate_with(seed, SPECTRA, engines, 0.0))
        });
        for (level, (mean, error)) in ["spectra", "peptides"].iter().zip(measured) {
            eprintln!(
                "{engines} engine(s), decoy sequences unpaired, {level}: mean FDP {:.3}% ± {:.3}%",
                mean * 100.0,
                error * 100.0
            );
            assert!(
                mean <= FDR + 2.0 * error,
                "{level}: mean FDP {mean:.4} above {FDR} (± {error:.4})"
            );
        }
    }
}

/// A retention-time error lifts wrong targets and decoys alike (neither
/// elutes where its sequence says), so the FDR must hold with it too.
#[test]
fn retention_times_keep_the_fdr_controlled() {
    for engines in [1, 2] {
        let measured = measure_searches(engines, Cut::Fdr, |seed| {
            with_retention_times(generate_with(seed, SPECTRA, engines, 0.0), seed)
        });
        for (level, (mean, error)) in ["spectra", "peptides"].iter().zip(measured) {
            eprintln!(
                "{engines} engine(s) with retention times, {level}: mean FDP {:.3}% ± {:.3}%",
                mean * 100.0,
                error * 100.0
            );
            assert!(
                mean <= FDR + 2.0 * error,
                "{level}: mean FDP {mean:.4} above {FDR} (± {error:.4})"
            );
        }
    }
}
