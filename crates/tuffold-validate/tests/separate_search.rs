//! Separate target and decoy searches (Proteome Discoverer's tables): each
//! spectrum has a best target and a best decoy match, both rank 1, and
//! validation competes them per spectrum after rescoring (what PD imports
//! do). That must keep the true FDP at 1% on average, as must competing them
//! first by the engine's main score (one rank 1, as a concatenated search
//! ranks them). Measured at 40 seeds, spectra then peptides: 0.50% ±
//! 0.09% and 0.31% ± 0.05% after rescoring (271 spectra accepted per
//! search), 0.57% ± 0.11% and 0.26% ± 0.06% by main score first (249).

mod support;

use std::collections::HashSet;

use support::{SplitMix, generate};
use tuffold_validate::{PsmRow, Settings, validate};

const SEEDS: u64 = 40;
const SPECTRA: u64 = 3000;
const FDR: f64 = 0.01;
const DECOY_OFFSET: i64 = 1_000_000;

/// One separate search: a correct spectrum's target is its true peptide
/// (scores shifted up as in `support`), an incorrect one's a random target;
/// every decoy match is a random decoy. Returns the rows and each
/// spectrum's true peptide.
fn separate(seed: u64) -> (tuffold_validate::Input, Vec<i64>) {
    let data = generate(seed, SPECTRA, 1);
    let mut input = data.input;
    let mut rng = SplitMix::new(seed ^ 0x5EA7_A7E5);
    let targets = (support::PROTEINS * support::PEPTIDES_PER_PROTEIN) as i64;
    let mut psms = Vec::new();
    for (index, truth) in data.true_peptides.iter().enumerate() {
        let spectrum_id = index as i64 + 1;
        let correct = input.psms[index].peptide_id == *truth;
        for decoy in [false, true] {
            let shift = if correct && !decoy { 1.0 } else { 0.0 };
            let peptide = if correct && !decoy {
                *truth
            } else {
                let id = rng.index(targets as u64) as i64 + 1;
                if decoy { DECOY_OFFSET + id } else { id }
            };
            psms.push(PsmRow {
                id: psms.len() as i64 + 1,
                spectrum_id,
                peptide_id: peptide,
                engine: 0,
                charge: 2,
                calculated_mz: None,
                experimental_mz: None,
                run_id: 1,
                retention_time: None,
                scores: vec![
                    (0, rng.normal() + 1.6 * shift),
                    (1, rng.normal() + 1.2 * shift),
                    (2, 10_f64.powf(-(rng.normal() + shift + 2.0))),
                ],
            });
        }
    }
    input.psms = psms;
    (input, data.true_peptides)
}

/// Keeps each spectrum's better match by the main score (f1).
fn concatenated(mut input: tuffold_validate::Input) -> tuffold_validate::Input {
    let rows = std::mem::take(&mut input.psms);
    input.psms = rows
        .chunks(2)
        .map(|pair| {
            if pair[1].scores[0].1 > pair[0].scores[0].1 {
                pair[1].clone()
            } else {
                pair[0].clone()
            }
        })
        .collect();
    input
}

/// Mean and standard error of the true FDP of spectra and of peptides, and
/// the spectra accepted per search.
fn measure(joint: bool) -> ([(f64, f64); 2], f64) {
    let mut spectra = Vec::new();
    let mut peptides = Vec::new();
    let mut accepted_total = 0;
    for seed in 1..=SEEDS {
        let (input, truth) = separate(seed);
        let input = if joint { concatenated(input) } else { input };
        let validation = validate(&input, &Settings::default()).unwrap();
        let accepted: Vec<_> = validation
            .spectra
            .iter()
            .filter(|m| !m.is_decoy && m.q_value <= FDR)
            .collect();
        let wrong = accepted
            .iter()
            .filter(|m| truth[usize::try_from(m.spectrum_id - 1).unwrap()] != m.peptide_id)
            .count();
        accepted_total += accepted.len();
        spectra.push(wrong as f64 / accepted.len().max(1) as f64);
        let present: HashSet<i64> = truth.iter().copied().collect();
        let accepted: Vec<_> = validation
            .peptides
            .iter()
            .filter(|p| !p.is_decoy && p.q_value <= FDR)
            .collect();
        let wrong = accepted
            .iter()
            .filter(|p| !present.contains(&p.peptide_id))
            .count();
        peptides.push(wrong as f64 / accepted.len().max(1) as f64);
    }
    let stats = [spectra, peptides].map(|values| {
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
        (mean, (variance / n).sqrt())
    });
    (stats, accepted_total as f64 / SEEDS as f64)
}

#[test]
fn both_competitions_control_the_fdr() {
    for (label, joint) in [("after rescoring", false), ("by main score first", true)] {
        let (stats, accepted) = measure(joint);
        eprintln!("{label}: {accepted:.0} spectra per search");
        for (level, (mean, error)) in ["spectra", "peptides"].iter().zip(stats) {
            eprintln!(
                "{label} {level}: mean FDP {:.3}% ± {:.3}%",
                mean * 100.0,
                error * 100.0
            );
            assert!(
                mean <= FDR + 2.0 * error,
                "{label} {level}: {mean:.4} ± {error:.4}"
            );
        }
    }
}
