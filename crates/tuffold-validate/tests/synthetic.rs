//! End-to-end validation of synthetic searches with known ground truth.
//!
//! The claims under test: q-values control the true false discovery
//! proportion, rescoring beats the best single score, engines agreeing
//! reinforce one another, and protein groups stay mostly true. mokapot
//! numbers come from `tools/crosscheck/synthetic.py` on the same data.

mod support;

use std::collections::{BTreeMap, HashMap, HashSet};

use support::{PRESENT, elution_times, generate, only_engine, with_retention_times};
use tuffold_validate::fdr::q_values;
use tuffold_validate::{
    BestScoreReason, Method, NoteCode, Settings, Thresholds, Validation, validate,
};

const SEED: u64 = 20_261_008;
const SPECTRA: u64 = 6000;

/// Accepted target spectra at `fdr` and how many of them are wrong.
fn accepted_spectra(validation: &Validation, true_peptides: &[i64], fdr: f64) -> (usize, usize) {
    let accepted: Vec<_> = validation
        .spectra
        .iter()
        .filter(|m| !m.is_decoy && m.q_value <= fdr)
        .collect();
    let wrong = accepted
        .iter()
        .filter(|m| {
            let index = usize::try_from(m.spectrum_id - 1).unwrap();
            true_peptides[index] != m.peptide_id
        })
        .count();
    (accepted.len(), wrong)
}

fn rescoring(rescore: bool) -> Settings {
    Settings {
        rescore,
        thresholds: Thresholds::default(),
    }
}

#[test]
fn one_engine_rescored_controls_the_true_fdr_and_beats_its_best_score() {
    let data = only_engine(generate(SEED, SPECTRA, 1), 0);
    let rescored = validate(&data.input, &rescoring(true)).unwrap();
    let single = validate(&data.input, &rescoring(false)).unwrap();

    match &rescored.engines[0].method {
        Method::Rescored { features, .. } => {
            assert_eq!(
                features,
                &["-log10(expect)", "f1", "f2", tuffold_validate::REPLICATES]
            );
        }
        other => panic!("expected rescoring, got {other:?}"),
    }
    assert!(rescored.notes.is_empty(), "{:?}", rescored.notes);
    match &single.engines[0].method {
        Method::BestScore {
            score,
            why,
            higher_is_better,
        } => {
            // f1 has the largest shift (1.6 standard deviations).
            assert_eq!(
                (score.as_str(), *why, *higher_is_better),
                ("f1", BestScoreReason::RescoringOff, true)
            );
        }
        other => panic!("expected the best score, got {other:?}"),
    }

    let (rescored_count, rescored_wrong) = accepted_spectra(&rescored, &data.true_peptides, 0.01);
    let (single_count, single_wrong) = accepted_spectra(&single, &data.true_peptides, 0.01);
    eprintln!(
        "rescored {rescored_count} ({rescored_wrong} wrong); f1 {single_count} ({single_wrong} wrong)"
    );
    assert!(
        rescored_count > single_count + single_count / 5,
        "rescoring should add >20%"
    );
    for (count, wrong) in [
        (rescored_count, rescored_wrong),
        (single_count, single_wrong),
    ] {
        let proportion = wrong as f64 / count as f64;
        assert!(proportion <= 0.02, "true FDP {proportion} at 1% FDR");
    }
}

#[test]
fn one_engine_holds_up_against_mokapot_and_the_oracle() {
    // `python tools/crosscheck/synthetic.py` on this exact data (its generator
    // is a draw-for-draw port; the checksum line confirms it):
    // - mokapot 0.10.0 (linear SVM, 3 folds) over seeds 1..10 accepts
    //   1018..1195 target PSMs at 1% FDR, median 1123;
    // - the oracle score (the generator's true weights) accepts 1328. That is
    //   the most one spectrum's own scores can give; replicate spectra add
    //   evidence from the peptide's other spectra, so Tuffold passes it, and
    //   only its true FDP (here and in tests/calibration.rs) says whether
    //   that is honest.
    const MOKAPOT_MEDIAN: usize = 1123;
    const ORACLE: usize = 1328;
    let data = only_engine(generate(SEED, SPECTRA, 1), 0);
    let validation = validate(&data.input, &rescoring(true)).unwrap();
    let (accepted, wrong) = accepted_spectra(&validation, &data.true_peptides, 0.01);
    assert!(
        accepted * 100 >= MOKAPOT_MEDIAN * 95,
        "Tuffold {accepted} vs mokapot median {MOKAPOT_MEDIAN}"
    );
    assert!(accepted > ORACLE, "Tuffold {accepted} vs oracle {ORACLE}");
    assert!(
        wrong as f64 / accepted as f64 <= 0.015,
        "true FDP {wrong}/{accepted} at 1% FDR"
    );
}

#[test]
fn agreeing_engines_reinforce_one_another() {
    let both = generate(SEED, SPECTRA, 2);
    let true_peptides = both.true_peptides.clone();
    let combined = validate(&both.input, &rescoring(true)).unwrap();
    let alone = validate(
        &only_engine(generate(SEED, SPECTRA, 2), 0).input,
        &rescoring(true),
    )
    .unwrap();

    let (combined_count, combined_wrong) = accepted_spectra(&combined, &true_peptides, 0.01);
    let (alone_count, _) = accepted_spectra(&alone, &true_peptides, 0.01);
    eprintln!("combined {combined_count} ({combined_wrong} wrong); engine 0 alone {alone_count}");
    assert!(
        combined_count > alone_count,
        "two engines should accept more"
    );
    assert!(combined_wrong as f64 / combined_count as f64 <= 0.02);
    // Correct spectra are found by both engines; most accepted matches
    // carry two votes.
    let two_votes = combined
        .spectra
        .iter()
        .filter(|m| !m.is_decoy && m.q_value <= 0.01 && m.engines == 2)
        .count();
    assert!(
        two_votes * 10 > combined_count * 9,
        "{two_votes} of {combined_count}"
    );
    assert_eq!(combined.engines.len(), 2);
    assert_eq!(combined.psms.len(), both.input.psms.len());
}

#[test]
fn peptide_and_protein_levels_stay_mostly_true() {
    let data = generate(SEED, SPECTRA, 2);
    let validation = validate(&data.input, &rescoring(true)).unwrap();
    let true_peptides: HashSet<i64> = data.true_peptides.iter().copied().collect();
    let accepted: Vec<_> = validation
        .peptides
        .iter()
        .filter(|p| !p.is_decoy && p.q_value <= 0.01)
        .collect();
    let wrong = accepted
        .iter()
        .filter(|p| !true_peptides.contains(&p.peptide_id))
        .count();
    eprintln!("peptides: {} accepted, {wrong} wrong", accepted.len());
    assert!(wrong as f64 / accepted.len() as f64 <= 0.03);

    let present = i64::try_from(PRESENT).unwrap();
    let groups: Vec<_> = validation
        .proteins
        .groups
        .iter()
        .filter(|g| !g.is_decoy && g.q_value.is_some_and(|q| q <= 0.01))
        .collect();
    let wrong = groups.iter().filter(|g| g.leader > present).count();
    eprintln!("protein groups: {} accepted, {wrong} wrong", groups.len());
    // 200 proteins are present; most are found (179 at the time of writing)
    // and few absent ones pass.
    assert!(groups.len() >= 170, "{}", groups.len());
    assert!(wrong as f64 / groups.len() as f64 <= 0.05);
    let leaders: HashMap<i64, usize> = groups
        .iter()
        .map(|g| (g.leader, g.peptides.len()))
        .collect();
    assert!(leaders.keys().filter(|&&leader| leader <= present).count() >= 170);
}

/// Correct peptides are matched by several spectra, random ones rarely:
/// pooling a peptide's spectra must find more peptides than its best
/// spectrum alone (q-values over best-spectrum scores, recomputed here), and
/// stay true (`tests/calibration.rs` checks the mean over many seeds).
#[test]
fn pooling_a_peptides_spectra_finds_more_peptides() {
    let data = generate(SEED, SPECTRA, 2);
    let validation = validate(&data.input, &rescoring(true)).unwrap();
    match &validation.peptide_method {
        // Every spectrum is 2+ here, so "other charges" is constant and
        // dropped.
        Some(Method::Rescored { features, .. }) => assert_eq!(
            features,
            &[
                "best spectrum score",
                "second-best spectrum score",
                "log(1 + other spectra)",
                "log(1 + sibling peptides)"
            ]
        ),
        other => panic!("expected pooled peptides, got {other:?}"),
    }
    let mut best: BTreeMap<i64, (f64, bool)> = BTreeMap::new();
    for found in &validation.spectra {
        let entry = best
            .entry(found.peptide_id)
            .or_insert((found.score, found.is_decoy));
        entry.0 = entry.0.max(found.score);
    }
    let scores: Vec<f64> = best.values().map(|(score, _)| *score).collect();
    let decoys: Vec<bool> = best.values().map(|(_, decoy)| *decoy).collect();
    let unpooled = q_values(&scores, &decoys)
        .iter()
        .zip(&decoys)
        .filter(|(q, decoy)| **q <= 0.01 && !**decoy)
        .count();
    let truth: HashSet<i64> = data.true_peptides.iter().copied().collect();
    let pooled: Vec<i64> = validation
        .peptides
        .iter()
        .filter(|p| !p.is_decoy && p.q_value <= 0.01)
        .map(|p| p.peptide_id)
        .collect();
    let wrong = pooled.iter().filter(|id| !truth.contains(id)).count();
    eprintln!(
        "pooled {} peptides ({wrong} wrong); best spectrum alone {unpooled}",
        pooled.len()
    );
    assert!(
        pooled.len() * 100 >= unpooled * 110,
        "{} vs {unpooled}",
        pooled.len()
    );
    assert!(wrong as f64 / pooled.len() as f64 <= 0.03);
}

/// Correct spectra elute where their peptide's sequence says, wrong ones
/// anywhere: the retention-time error must become a feature and add
/// identifications.
#[test]
fn retention_times_add_identifications() {
    let plain = only_engine(generate(SEED, SPECTRA, 1), 0);
    let timed = with_retention_times(only_engine(generate(SEED, SPECTRA, 1), 0), SEED);
    let without = validate(&plain.input, &rescoring(true)).unwrap();
    let with = validate(&timed.input, &rescoring(true)).unwrap();
    match &with.engines[0].method {
        Method::Rescored { features, .. } => assert!(
            features.iter().any(|f| f.starts_with("|RT error|")),
            "{features:?}"
        ),
        other => panic!("expected rescoring, got {other:?}"),
    }
    let (before, _) = accepted_spectra(&without, &plain.true_peptides, 0.01);
    let (after, wrong) = accepted_spectra(&with, &timed.true_peptides, 0.01);
    eprintln!("with retention times {after} ({wrong} wrong), without {before}");
    assert!(after > before, "{after} vs {before}");
    assert!(wrong as f64 / after as f64 <= 0.02);
}

/// Spectra of peptides other spectra confidently identify are lifted; wrong
/// spectra on present peptides are lifted as often for decoys, whose paired
/// targets' spectra count, so the error rate holds.
#[test]
fn replicate_spectra_add_identifications() {
    let data = only_engine(generate(SEED, SPECTRA, 1), 0);
    let with = validate(&data.input, &rescoring(true)).unwrap();
    match &with.engines[0].method {
        Method::Rescored { features, .. } => {
            assert!(
                features.iter().any(|f| f == tuffold_validate::REPLICATES),
                "{features:?}"
            );
        }
        other => panic!("expected rescoring, got {other:?}"),
    }
    let (after, wrong) = accepted_spectra(&with, &data.true_peptides, 0.01);
    eprintln!("with replicates {after} ({wrong} wrong)");
    assert!(wrong as f64 / after as f64 <= 0.02);
}

/// Decoys that can't be paired with a target (shuffled, not reversed) get
/// replicate counts drawn from low-scoring targets', and say so.
#[test]
fn unpaired_decoys_get_imputed_replicates() {
    let data = support::with_unpaired_decoy_sequences(only_engine(generate(SEED, SPECTRA, 1), 0));
    let validation = validate(&data.input, &rescoring(true)).unwrap();
    match &validation.engines[0].method {
        Method::Rescored { features, .. } => {
            assert!(
                !features.iter().any(|f| f == tuffold_validate::REPLICATES),
                "{features:?}"
            );
            assert!(
                features
                    .iter()
                    .any(|f| f == tuffold_validate::REPLICATES_IMPUTED),
                "{features:?}"
            );
        }
        other => panic!("expected rescoring, got {other:?}"),
    }
    assert!(
        validation
            .notes
            .iter()
            .any(|note| note.code == tuffold_validate::NoteCode::ReplicatesImputed),
        "{:?}",
        validation.notes
    );
    let (accepted, wrong) = accepted_spectra(&validation, &data.true_peptides, 0.01);
    assert!(
        wrong as f64 / accepted as f64 <= 0.02,
        "{wrong} of {accepted}"
    );
}

/// A retention-time predictor that knows each peptide's own elution time,
/// as a fine-tuned deep model nearly does.
struct Oracle(std::collections::HashMap<i64, f64>);

struct Trained<'a>(&'a std::collections::HashMap<i64, f64>);

impl tuffold_validate::TrainedRetention for Trained<'_> {
    fn predict(&self, peptides: &[i64]) -> Vec<Option<f64>> {
        peptides.iter().map(|id| self.0.get(id).copied()).collect()
    }
}

impl tuffold_validate::RetentionPredictor for Oracle {
    fn train(&self, _: &[(i64, f64)]) -> Option<Box<dyn tuffold_validate::TrainedRetention + '_>> {
        Some(Box::new(Trained(&self.0)))
    }
}

/// A better retention-time model replaces the classical one (its feature
/// named as such), adds identifications, and keeps the error rate.
#[test]
fn a_deep_retention_model_replaces_the_classical_one() {
    let timed = with_retention_times(only_engine(generate(SEED, SPECTRA, 1), 0), SEED);
    let oracle = Oracle(elution_times(&timed, SEED));
    let classical = validate(&timed.input, &rescoring(true)).unwrap();
    let deep =
        tuffold_validate::validate_with(&timed.input, &rescoring(true), Some(&oracle)).unwrap();
    match &deep.engines[0].method {
        Method::Rescored { features, .. } => {
            assert!(
                features
                    .iter()
                    .any(|f| f == tuffold_validate::DEEP_RT_ERROR),
                "{features:?}"
            );
            assert!(
                !features.iter().any(|f| f == tuffold_validate::RT_ERROR),
                "{features:?}"
            );
        }
        other => panic!("expected rescoring, got {other:?}"),
    }
    let (before, _) = accepted_spectra(&classical, &timed.true_peptides, 0.01);
    let (after, wrong) = accepted_spectra(&deep, &timed.true_peptides, 0.01);
    eprintln!("deep model {after} ({wrong} wrong), classical {before}");
    assert!(after >= before, "{after} vs {before}");
    assert!(wrong as f64 / after as f64 <= 0.02);
}

/// Features predicted outside validation join rescoring; PSMs without them
/// get the lowest value, and the error rate holds.
#[test]
fn predicted_features_join_rescoring() {
    let mut data = only_engine(generate(SEED, SPECTRA, 1), 0);
    let plain = validate(&data.input, &rescoring(true)).unwrap();
    // A noisy signal of correctness (decoys are never correct), missing for
    // one PSM in twenty.
    let mut rng = support::SplitMix::new(SEED ^ 0xFEED);
    let mut values = std::collections::HashMap::new();
    for psm in &data.input.psms {
        let correct =
            data.true_peptides[usize::try_from(psm.spectrum_id - 1).unwrap()] == psm.peptide_id;
        let value = if correct { 1.0 } else { 0.0 } + 0.6 * rng.normal();
        if psm.id % 20 != 0 {
            values.insert(psm.id, vec![value]);
        }
    }
    data.input.predicted = Some(tuffold_validate::PredictedFeatures {
        names: vec!["predicted signal".into()],
        values,
    });
    let with = validate(&data.input, &rescoring(true)).unwrap();
    match &with.engines[0].method {
        Method::Rescored { features, .. } => {
            assert!(
                features.iter().any(|f| f == "predicted signal"),
                "{features:?}"
            );
        }
        other => panic!("expected rescoring, got {other:?}"),
    }
    let (before, _) = accepted_spectra(&plain, &data.true_peptides, 0.01);
    let (after, wrong) = accepted_spectra(&with, &data.true_peptides, 0.01);
    eprintln!("with predicted features {after} ({wrong} wrong), without {before}");
    assert!(after > before, "{after} vs {before}");
    assert!(wrong as f64 / after as f64 <= 0.02);
}

#[test]
fn without_decoys_nothing_is_estimated() {
    let mut data = generate(SEED, 500, 1);
    let decoys: HashSet<i64> = data
        .input
        .peptides
        .iter()
        .filter(|p| p.is_decoy)
        .map(|p| p.id)
        .collect();
    data.input
        .psms
        .retain(|psm| !decoys.contains(&psm.peptide_id));
    let validation = validate(&data.input, &Settings::default()).unwrap();
    assert!(validation.spectra.is_empty() && validation.peptides.is_empty());
    let codes: Vec<NoteCode> = validation.notes.iter().map(|n| n.code).collect();
    assert_eq!(codes, [NoteCode::EngineUnscored, NoteCode::NoDecoys]);
    assert!(matches!(
        validation.engines[0].method,
        Method::Unscored { .. }
    ));
}

#[test]
fn small_searches_fall_back_to_the_best_score_and_say_so() {
    let data = only_engine(generate(SEED, 40, 1), 0);
    let validation = validate(&data.input, &Settings::default()).unwrap();
    assert!(matches!(
        validation.engines[0].method,
        Method::BestScore {
            why: BestScoreReason::TooFewConfident,
            ..
        }
    ));
    assert_eq!(validation.notes[0].code, NoteCode::RescoringFallback);
}
