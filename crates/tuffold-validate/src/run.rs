//! One validation pass over a project's PSMs.
//!
//! 1. **Engines.** Each search engine's rank-1 PSMs are ranked on their own
//!    ([`crate::rescore`]): competed per spectrum, with q-values and PEPs.
//! 2. **Spectra.** Engines vote for a peptide per spectrum with their
//!    calibrated score: `(s - median) / MAD` of the engine's decoy scores
//!    ([`crate::rescore`] calibrates rescored folds that way, and engines
//!    ranked by one score get the same treatment), so every engine speaks on
//!    one decoy-anchored scale. A peptide's support is its *strongest* vote;
//!    the best-supported peptide wins (ties to a decoy), and spectrum-level
//!    q-values and PEPs follow target-decoy competition over the winners.
//!
//!    Votes are not summed. Engines that agree on a wrong target (a
//!    near-miss: a homologous peptide, a misplaced modification) would
//!    reinforce each other, and decoys can never model that, so summing
//!    breaks FDR control. Measured on synthetic data with 5% shared
//!    near-misses (200 seeds): one engine alone reaches a 2.1% true FDP at 1%
//!    FDR — a limit of target-decoy validation itself — the strongest vote
//!    2.3%, summing 2.7% (see `tests/calibration.rs`). Engines still fill
//!    each other's gaps: a spectrum only one engine identifies keeps its vote. (Votes are scores rather than PEPs because isotonic PEPs
//!    are flat at the top, which would tie every confident spectrum.)
//! 3. **Peptides.** Each peptide is scored by its best spectrum; q-values
//!    and PEPs follow over peptides.
//! 4. **Proteins.** [`crate::inference`] at the peptide threshold.

use std::collections::{BTreeMap, HashMap};

use crate::fdr::{confidence, peps_from};
use crate::features::Features;
use crate::inference::{PeptideEvidence, ProteinInference, infer_proteins};
use crate::input::{Input, PeptideRow, PredictedFeatures, PsmRow};
use crate::notes::{Note, NoteCode};
use crate::pooling::{Support, group_by_peptide, peptide_scores};
use crate::replicates::{REPLICATES, REPLICATES_IMPUTED, replicate_feature};
use crate::rescore::{
    BestScoreReason, FOLDS, Method, Rows, TRAIN_FDR, UnscoredReason, best_single_scores, calibrate,
    fold_of, score_engine, sequence_fold,
};
use crate::retention::{
    DEEP_RT_ERROR, LOG_DEEP_RT_ERROR, LOG_RT_ERROR, RT_ERROR, RetentionPredictor, rt_errors,
};
use crate::settings::Settings;
use crate::siblings::sibling_counts;

#[derive(Debug, Clone, PartialEq)]
pub struct EngineSummary {
    pub name: String,
    /// Rank-1 PSMs, decoys included.
    pub psms: usize,
    pub decoys: usize,
    pub method: Method,
    /// Targets passing 1% FDR on this engine alone.
    pub accepted: usize,
}

/// A rank-1 PSM ranked within its engine.
#[derive(Debug, Clone, PartialEq)]
pub struct PsmResult {
    pub psm_id: i64,
    pub score: f64,
    /// `None` when another PSM of the same engine won the spectrum.
    pub q_value: Option<f64>,
    pub pep: Option<f64>,
}

/// The peptide a spectrum is assigned after the engines vote.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumMatch {
    pub spectrum_id: i64,
    pub peptide_id: i64,
    pub is_decoy: bool,
    pub score: f64,
    pub q_value: f64,
    pub pep: f64,
    /// Engines whose top PSM is this peptide.
    pub engines: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PeptideResult {
    pub peptide_id: i64,
    pub is_decoy: bool,
    pub score: f64,
    /// The score its protein group ranks by: `score`, except when decoys'
    /// sibling counts were imputed, then the peptide's score without siblings.
    /// An imputed count stands for a wrong match landing in a present
    /// protein, which joins that protein's group; on a decoy it would make a
    /// group of its own that no wrong target could make.
    pub group_score: f64,
    pub q_value: f64,
    pub pep: f64,
    pub best_spectrum_id: i64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Validation {
    pub engines: Vec<EngineSummary>,
    pub psms: Vec<PsmResult>,
    pub spectra: Vec<SpectrumMatch>,
    pub peptides: Vec<PeptideResult>,
    /// How peptides were scored ([`crate::pooling`]); `None` when no engine
    /// was ranked.
    pub peptide_method: Option<Method>,
    pub proteins: ProteinInference,
    pub notes: Vec<Note>,
    /// PTM site localization outcomes for validated target spectra with
    /// ambiguous variable modifications (`(spectrum_id, outcome)`).
    pub localizations: Vec<(i64, crate::localize::SpectrumLocalization)>,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ValidationError {
    #[error("{0}")]
    Thresholds(String),
    #[error("PSM {psm} refers to {what} that the project does not contain")]
    Inconsistent { psm: i64, what: String },
}

/// An engine's winners after ranking, as indices into `Input::psms`.
struct Ranked {
    rows: Vec<usize>,
    scores: Vec<f64>,
    q: Vec<Option<f64>>,
    pep: Vec<Option<f64>>,
    /// Scores and q-values before replicate spectra ([`crate::replicates`]),
    /// which peptides are scored from: a wrong spectrum on a present peptide
    /// is no wrong peptide, so decoys lifted for that would only make
    /// peptides harder to accept.
    peptide_scores: Vec<f64>,
    peptide_q: Vec<Option<f64>>,
}

/// Which of an engine's scores a vote uses.
#[derive(Clone, Copy)]
enum Pass {
    Spectra,
    Peptides,
}

pub fn validate(input: &Input, settings: &Settings) -> Result<Validation, ValidationError> {
    validate_with(input, settings, None)
}

/// [`validate`], with a deep-learned retention-time model to use in place of
/// the classical one.
pub fn validate_with(
    input: &Input,
    settings: &Settings,
    retention: Option<&dyn RetentionPredictor>,
) -> Result<Validation, ValidationError> {
    settings
        .thresholds
        .check()
        .map_err(ValidationError::Thresholds)?;
    let peptides = index_peptides(input)?;
    let mut validation = Validation::default();
    let mut ranked = Vec::new();
    for (engine, name) in input.engines.iter().enumerate() {
        let rows: Vec<usize> = (0..input.psms.len())
            .filter(|&row| input.psms[row].engine == engine)
            .collect();
        if rows.is_empty() {
            continue;
        }
        let (summary, result, note) =
            rank_engine(input, &peptides, name, rows, settings.rescore, retention);
        validation.notes.extend(note);
        validation
            .notes
            .extend(engine_note(&summary, settings.rescore));
        validation.notes.extend(few_decoys_note(&summary));
        validation.engines.push(summary);
        ranked.extend(result);
    }
    let no_decoys_anywhere = validation.engines.iter().all(|engine| {
        matches!(
            engine.method,
            Method::Unscored {
                why: UnscoredReason::NoDecoys
            }
        )
    });
    if ranked.is_empty() && no_decoys_anywhere {
        validation.notes.push(Note::new(
            NoteCode::NoDecoys,
            "No engine has decoy PSMs with a usable score, so no FDR can be estimated. Search \
             against a target-decoy database and import the unfiltered results.",
        ));
    }
    if ranked.is_empty() {
        return Ok(validation);
    }
    validation.psms = psm_results(input, &ranked);
    validation.spectra = spectrum_matches(input, &peptides, &ranked, Pass::Spectra);
    let peptide_votes = spectrum_matches(input, &peptides, &ranked, Pass::Peptides);
    let (peptide_results, peptide_method, note) =
        peptide_results(input, &peptide_votes, settings.rescore);
    validation.notes.extend(note);
    validation.peptides = peptide_results;
    validation.peptide_method = Some(peptide_method);
    let evidence: Vec<PeptideEvidence<'_>> = validation
        .peptides
        .iter()
        .map(|result| PeptideEvidence {
            id: result.peptide_id,
            sequence: &peptides[&result.peptide_id].sequence,
            is_decoy: result.is_decoy,
            score: result.group_score,
            q_value: result.q_value,
            pep: result.pep,
        })
        .collect();
    validation.proteins = infer_proteins(
        &evidence,
        &input.peptide_proteins,
        &input.proteins,
        &settings.thresholds,
    );
    Ok(validation)
}

/// A rank-1 PSM's best single engine score: higher is better within its
/// engine (an index into [`Input::engines`]); engines aren't on one scale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EngineScore {
    pub engine: usize,
    pub score: f64,
}

/// Each rank-1 PSM's [`EngineScore`], by PSM id. The decoys choose only
/// which score and direction, as when an engine is ranked without
/// rescoring; engines without decoys or a shared score have none.
pub fn engine_scores(input: &Input) -> Result<HashMap<i64, EngineScore>, ValidationError> {
    let peptides = index_peptides(input)?;
    let mut scores = HashMap::new();
    for engine in 0..input.engines.len() {
        let psms: Vec<&PsmRow> = input
            .psms
            .iter()
            .filter(|psm| psm.engine == engine)
            .collect();
        let peptide_rows: Vec<&PeptideRow> =
            psms.iter().map(|psm| peptides[&psm.peptide_id]).collect();
        let decoys: Vec<bool> = peptide_rows.iter().map(|p| p.is_decoy).collect();
        if !decoys.contains(&true) {
            continue;
        }
        let spectra: Vec<i64> = psms.iter().map(|psm| psm.spectrum_id).collect();
        let folds: Vec<i64> = peptide_rows
            .iter()
            .map(|p| sequence_fold(&p.sequence))
            .collect();
        let features = Features::build(&psms, &peptide_rows, &input.score_names);
        if features.engine_scores == 0 {
            continue;
        }
        let view = Rows {
            features: &features,
            spectra: &spectra,
            folds: &folds,
            decoys: &decoys,
        };
        let Some(best) = best_single_scores(&view) else {
            continue;
        };
        scores.extend(
            psms.iter()
                .zip(best)
                .map(|(psm, score)| (psm.id, EngineScore { engine, score })),
        );
    }
    Ok(scores)
}

fn index_peptides(input: &Input) -> Result<HashMap<i64, &PeptideRow>, ValidationError> {
    let peptides: HashMap<i64, &PeptideRow> = input.peptides.iter().map(|p| (p.id, p)).collect();
    for psm in &input.psms {
        let missing = if !peptides.contains_key(&psm.peptide_id) {
            Some(format!("peptide {}", psm.peptide_id))
        } else if psm.engine >= input.engines.len() {
            Some(format!("engine {}", psm.engine))
        } else {
            psm.scores
                .iter()
                .find(|(name, _)| *name >= input.score_names.len())
                .map(|(name, _)| format!("score {name}"))
        };
        if let Some(what) = missing {
            return Err(ValidationError::Inconsistent { psm: psm.id, what });
        }
    }
    Ok(peptides)
}

/// Predicted features as derived columns; a PSM without them gets each
/// column's lowest value (they are all higher-is-better).
fn add_predicted(features: &mut Features, psms: &[&PsmRow], predicted: &PredictedFeatures) {
    for (column, name) in predicted.names.iter().enumerate() {
        let values: Vec<Option<f64>> = psms
            .iter()
            .map(|psm| {
                predicted
                    .values
                    .get(&psm.id)
                    .and_then(|values| values.get(column))
                    .copied()
                    .filter(|value| value.is_finite())
            })
            .collect();
        let Some(lowest) = values.iter().flatten().copied().min_by(f64::total_cmp) else {
            continue;
        };
        features.add_derived(
            name,
            values.into_iter().map(|v| v.unwrap_or(lowest)).collect(),
        );
    }
}

fn rank_engine(
    input: &Input,
    peptides: &HashMap<i64, &PeptideRow>,
    name: &str,
    rows: Vec<usize>,
    rescore: bool,
    retention: Option<&dyn RetentionPredictor>,
) -> (EngineSummary, Option<Ranked>, Option<Note>) {
    let psms: Vec<_> = rows.iter().map(|&row| &input.psms[row]).collect();
    let peptide_rows: Vec<&PeptideRow> = psms.iter().map(|psm| peptides[&psm.peptide_id]).collect();
    let decoys: Vec<bool> = peptide_rows.iter().map(|p| p.is_decoy).collect();
    let spectra: Vec<i64> = psms.iter().map(|psm| psm.spectrum_id).collect();
    let folds: Vec<i64> = peptide_rows
        .iter()
        .map(|p| sequence_fold(&p.sequence))
        .collect();
    let decoy_count = decoys.iter().filter(|&&decoy| decoy).count();
    let psm_count = rows.len();
    let summary = |method: Method, accepted: usize| EngineSummary {
        name: name.to_string(),
        psms: psm_count,
        decoys: decoy_count,
        method,
        accepted,
    };
    let mut features = Features::build(&psms, &peptide_rows, &input.score_names);
    if rescore && let Some(predicted) = &input.predicted {
        add_predicted(&mut features, &psms, predicted);
    }
    let unscored = if decoy_count == 0 {
        Some(UnscoredReason::NoDecoys)
    } else if features.engine_scores == 0 {
        Some(UnscoredReason::NoSharedScore)
    } else {
        None
    };
    if let Some(why) = unscored {
        return (summary(Method::Unscored { why }, 0), None, None);
    }
    if rescore {
        let view = Rows {
            features: &features,
            spectra: &spectra,
            folds: &folds,
            decoys: &decoys,
        };
        if let Some((errors, deep)) = rt_errors(&view, &psms, &peptide_rows, retention) {
            let logged = errors.iter().map(|error| error.ln_1p()).collect();
            let (name, log_name) = if deep {
                (DEEP_RT_ERROR, LOG_DEEP_RT_ERROR)
            } else {
                (RT_ERROR, LOG_RT_ERROR)
            };
            features.add_derived(name, errors);
            features.add_derived(log_name, logged);
        }
    }
    let view = Rows {
        features: &features,
        spectra: &spectra,
        folds: &folds,
        decoys: &decoys,
    };
    let (scores, method) = score_engine(&view, rescore);
    let all: Vec<usize> = (0..rows.len()).collect();
    let first_scores = match method {
        Method::BestScore { .. } => {
            calibrate(&view, &all, &scores).unwrap_or_else(|| scores.clone())
        }
        Method::Rescored { .. } | Method::Unscored { .. } => scores.clone(),
    };
    let first_q = crate::rescore::competed_q(&view, &all, &first_scores);
    let (features, scores, method, note) = if rescore {
        with_replicates(
            features,
            &peptide_rows,
            &spectra,
            &folds,
            &decoys,
            (scores, method),
            name,
        )
    } else {
        (features, scores, method, None)
    };
    let view = Rows {
        features: &features,
        spectra: &spectra,
        folds: &folds,
        decoys: &decoys,
    };
    // Rescored scores are calibrated per fold already; put single-score
    // rankings on the same decoy-anchored scale so engines can be compared.
    let scores = match method {
        Method::BestScore { .. } => calibrate(&view, &all, &scores).unwrap_or(scores),
        Method::Rescored { .. } | Method::Unscored { .. } => scores,
    };
    let q = crate::rescore::competed_q(&view, &all, &scores);
    let pep = cross_fit_peps(&decoys, &scores, &q, &folds);
    let accepted = (0..rows.len())
        .filter(|&i| !decoys[i] && q[i].is_some_and(|q| q <= TRAIN_FDR))
        .count();
    let ranked = Ranked {
        rows,
        scores,
        q,
        pep,
        peptide_scores: first_scores,
        peptide_q: first_q,
    };
    (summary(method, accepted), Some(ranked), note)
}

/// Rescore again with each PSM's replicate spectra (decoys: their paired
/// target's too); the first pass stands only when rescoring falls back.
fn with_replicates(
    features: Features,
    peptides: &[&PeptideRow],
    spectra: &[i64],
    folds: &[i64],
    decoys: &[bool],
    (scores, method): (Vec<f64>, Method),
    engine: &str,
) -> (Features, Vec<f64>, Method, Option<Note>) {
    let all: Vec<usize> = (0..spectra.len()).collect();
    let view = Rows {
        features: &features,
        spectra,
        folds,
        decoys,
    };
    let q = crate::rescore::competed_q(&view, &all, &scores);
    let replicates = replicate_feature(peptides, spectra, decoys, &scores, &q);
    let note = replicates.imputed.map(|unpaired| {
        Note::new(
            NoteCode::ReplicatesImputed,
            format!(
                "{engine}: decoys' sequences don't name their target peptides: only {:.1}% of \
                 decoys had a paired target with confident spectra, against {:.1}% of \
                 low-scoring targets. So decoys' replicate counts were drawn from those of \
                 low-scoring targets, which are nearly all wrong, for the decoys to count the \
                 errors that evidence adds.",
                100.0 * unpaired.decoys,
                100.0 * unpaired.low_targets
            ),
        )
    });
    let name = if note.is_some() {
        REPLICATES_IMPUTED
    } else {
        REPLICATES
    };
    let mut pooled = features.clone();
    pooled.add_derived(name, replicates.column);
    let view = Rows {
        features: &pooled,
        spectra,
        folds,
        decoys,
    };
    // Not chosen by which pass accepts more: that picks favourable decoy
    // noise (synthetic one-engine true FDP 1.16% → 1.26%).
    let (pooled_scores, pooled_method) = score_engine(&view, true);
    if matches!(pooled_method, Method::Rescored { .. }) {
        return (pooled, pooled_scores, pooled_method, note);
    }
    (features, scores, method, None)
}

/// PEPs of the competition winners (those with a q-value), each fitted on
/// the other folds' winners. An in-sample fit lets an entry's own label set
/// its PEP — targets pull theirs down, decoys up — which would bias the
/// engine vote towards targets.
fn cross_fit_peps(
    decoys: &[bool],
    scores: &[f64],
    q: &[Option<f64>],
    folds: &[i64],
) -> Vec<Option<f64>> {
    let winners: Vec<usize> = (0..scores.len()).filter(|&i| q[i].is_some()).collect();
    let mut pep = vec![None; scores.len()];
    for fold in 0..FOLDS {
        let (held_out, train): (Vec<usize>, Vec<usize>) =
            winners.iter().partition(|&&i| fold_of(folds[i]) == fold);
        let train_scores: Vec<f64> = train.iter().map(|&i| scores[i]).collect();
        let train_decoys: Vec<bool> = train.iter().map(|&i| decoys[i]).collect();
        let query: Vec<f64> = held_out.iter().map(|&i| scores[i]).collect();
        for (index, value) in held_out
            .iter()
            .zip(peps_from(&train_scores, &train_decoys, &query))
        {
            pep[*index] = Some(value);
        }
    }
    pep
}

/// Decoys below this share of an engine's rank-1 PSMs suggest results
/// already filtered. Unfiltered searches have a decoy in roughly a third of
/// rank-1 PSMs (Sage on PXD006447 41%, MSFragger 27%, Mascot via Scaffold
/// 25%); a list cut at 1% FDR has about 1% (Scaffold Q+S MSFragger file,
/// PXD068782: 35 of 21,300).
const FEW_DECOYS_SHARE: f64 = 0.05;

fn few_decoys_note(summary: &EngineSummary) -> Option<Note> {
    let share = summary.decoys as f64 / summary.psms.max(1) as f64;
    (summary.decoys > 0 && share < FEW_DECOYS_SHARE).then(|| {
        Note::new(
            NoteCode::FewDecoys,
            format!(
                "Only {} of {}'s {} rank-1 PSMs are decoys ({:.1}%; an unfiltered search has \
                 about a third): its results look already filtered by the program that wrote \
                 them, so Tuffold's FDR is estimated on what that filter kept and checks \
                 little. Import the unfiltered results where you can.",
                summary.decoys,
                summary.name,
                summary.psms,
                share * 100.0
            ),
        )
    })
}

fn engine_note(summary: &EngineSummary, rescore: bool) -> Option<Note> {
    match &summary.method {
        Method::Unscored { why } => Some(Note::new(
            NoteCode::EngineUnscored,
            format!(
                "{} PSMs were not validated: {}.",
                summary.name,
                why.describe()
            ),
        )),
        Method::BestScore { score, why, .. }
            if rescore && *why != BestScoreReason::RescoringOff =>
        {
            Some(Note::new(
                NoteCode::RescoringFallback,
                format!(
                    "{} was ranked by {score} alone: {}.",
                    summary.name,
                    why.describe()
                ),
            ))
        }
        Method::BestScore { .. } | Method::Rescored { .. } => None,
    }
}

fn psm_results(input: &Input, ranked: &[Ranked]) -> Vec<PsmResult> {
    let mut results: Vec<PsmResult> = ranked
        .iter()
        .flat_map(|engine| {
            engine.rows.iter().enumerate().map(|(i, &row)| PsmResult {
                psm_id: input.psms[row].id,
                score: engine.scores[i],
                q_value: engine.q[i],
                pep: engine.pep[i],
            })
        })
        .collect();
    results.sort_by_key(|result| result.psm_id);
    results
}

fn spectrum_matches(
    input: &Input,
    peptides: &HashMap<i64, &PeptideRow>,
    ranked: &[Ranked],
    pass: Pass,
) -> Vec<SpectrumMatch> {
    // (spectrum, peptide) -> (strongest vote, engines voting for it)
    let mut votes: BTreeMap<(i64, i64), (f64, u32)> = BTreeMap::new();
    for engine in ranked {
        let (scores, q) = match pass {
            Pass::Spectra => (&engine.scores, &engine.q),
            Pass::Peptides => (&engine.peptide_scores, &engine.peptide_q),
        };
        for (i, &row) in engine.rows.iter().enumerate() {
            // Rows that lost their spectrum within the engine do not vote.
            if q[i].is_none() {
                continue;
            }
            let psm = &input.psms[row];
            let support = scores[i];
            let vote = votes
                .entry((psm.spectrum_id, psm.peptide_id))
                .or_insert((0.0, 0));
            vote.0 = if vote.1 == 0 {
                support
            } else {
                vote.0.max(support)
            };
            vote.1 += 1;
        }
    }
    let mut chosen: BTreeMap<i64, SpectrumMatch> = BTreeMap::new();
    for ((spectrum_id, peptide_id), (score, engines)) in votes {
        let is_decoy = peptides[&peptide_id].is_decoy;
        let candidate = SpectrumMatch {
            spectrum_id,
            peptide_id,
            is_decoy,
            score,
            q_value: 1.0,
            pep: 1.0,
            engines,
        };
        chosen
            .entry(spectrum_id)
            .and_modify(|current| {
                if score > current.score
                    || (score == current.score && is_decoy && !current.is_decoy)
                {
                    *current = candidate.clone();
                }
            })
            .or_insert(candidate);
    }
    let mut matches: Vec<SpectrumMatch> = chosen.into_values().collect();
    let scores: Vec<f64> = matches.iter().map(|m| m.score).collect();
    let decoys: Vec<bool> = matches.iter().map(|m| m.is_decoy).collect();
    for (found, (q, pep)) in matches.iter_mut().zip(confidence(&scores, &decoys)) {
        found.q_value = q;
        found.pep = pep;
    }
    matches
}

/// Peptides scored by their spectra, pooled ([`crate::pooling`]) when
/// rescoring is on.
fn peptide_results(
    input: &Input,
    spectra: &[SpectrumMatch],
    rescore: bool,
) -> (Vec<PeptideResult>, Method, Option<Note>) {
    // A spectrum's charge as its engines reported it for the winning peptide.
    let mut charges: HashMap<(i64, i64), i32> = HashMap::new();
    for psm in &input.psms {
        charges
            .entry((psm.spectrum_id, psm.peptide_id))
            .or_insert(psm.charge);
    }
    let grouped = group_by_peptide(spectra.iter().map(|found| {
        let support = Support {
            spectrum_id: found.spectrum_id,
            charge: charges
                .get(&(found.spectrum_id, found.peptide_id))
                .copied()
                .unwrap_or(0),
            score: found.score,
        };
        (found.peptide_id, found.is_decoy, support)
    }));
    let siblings = sibling_counts(&grouped, input);
    let (scores, method) = peptide_scores(&grouped, Some(&siblings), rescore);
    let note = match siblings.imputed {
        Some(unpaired) if rescore => Some(Note::new(
            NoteCode::SiblingsImputed,
            format!(
                "Decoy proteins don't name their targets the way the targets are named (for \
                 example rev_X for X): only {:.1}% of decoy peptides had a sibling through \
                 their names, against {:.1}% of low-scoring targets. So decoys' sibling counts \
                 were drawn from those of low-scoring targets, which are nearly all wrong, for \
                 the decoys to count the errors sibling evidence adds.",
                unpaired.decoys * 100.0,
                unpaired.low_targets * 100.0
            ),
        )),
        _ => None,
    };
    let decoys: Vec<bool> = grouped.iter().map(|p| p.is_decoy).collect();
    let group_scores = if siblings.imputed.is_some() {
        peptide_scores(&grouped, None, rescore).0
    } else {
        scores.clone()
    };
    let results = grouped
        .iter()
        .zip(&scores)
        .zip(&group_scores)
        .zip(confidence(&scores, &decoys))
        .map(
            |(((peptide, &score), &group_score), (q_value, pep))| PeptideResult {
                peptide_id: peptide.peptide_id,
                is_decoy: peptide.is_decoy,
                score,
                group_score,
                q_value,
                pep,
                best_spectrum_id: peptide.spectra[0].spectrum_id,
            },
        )
        .collect();
    (results, method, note)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::PsmRow;

    fn psm(spectrum_id: i64, peptide_id: i64, engine: usize) -> PsmRow {
        PsmRow {
            id: spectrum_id * 10 + i64::try_from(engine).unwrap(),
            spectrum_id,
            peptide_id,
            engine,
            charge: 2,
            calculated_mz: None,
            experimental_mz: None,
            run_id: 1,
            retention_time: None,
            scores: Vec::new(),
        }
    }

    fn peptide(id: i64, is_decoy: bool) -> PeptideRow {
        PeptideRow {
            id,
            sequence: "PEPTIDEK".into(),
            modifications: 0,
            is_decoy,
        }
    }

    #[test]
    fn engine_scores_follow_the_best_single_score_and_skip_engines_without_decoys() {
        // Engine A: a lower-is-better score, 200 targets at 1..=200 and 20
        // decoys at 201..=220. Negated, all 200 targets pass 1% ((0 + 1) /
        // 200 = 0.5%); as given, none do (the decoys come first), so the
        // best single score is the negation. Engine B has no decoys and
        // gets no scores.
        let scored = |spectrum: i64, peptide_id: i64, engine: usize, value: f64| PsmRow {
            scores: vec![(0, value)],
            ..psm(spectrum, peptide_id, engine)
        };
        let mut psms: Vec<PsmRow> = (1..=220)
            .map(|i| scored(i, if i <= 200 { 10 } else { 20 }, 0, i as f64))
            .collect();
        psms.push(scored(1, 10, 1, 5.0));
        let input = Input {
            engines: vec!["A".into(), "B".into()],
            score_names: vec!["sprank".into()],
            psms,
            peptides: vec![peptide(10, false), peptide(20, true)],
            ..Input::default()
        };
        let scores = engine_scores(&input).unwrap();
        assert_eq!(scores.len(), 220);
        assert!(!scores.contains_key(&11));
        let best_target = scores[&10];
        assert_eq!(best_target.engine, 0);
        assert!(best_target.score > scores[&2000].score);
        assert!(scores[&2000].score > scores[&2200].score);
    }

    #[test]
    fn an_engine_with_few_decoys_is_flagged_as_filtered() {
        let summary = |psms, decoys| EngineSummary {
            name: "MSFragger".into(),
            psms,
            decoys,
            method: Method::BestScore {
                score: "hyperscore".into(),
                higher_is_better: true,
                why: BestScoreReason::RescoringOff,
            },
            accepted: 0,
        };
        let note = few_decoys_note(&summary(21_300, 35)).expect("flagged");
        assert_eq!(note.code, NoteCode::FewDecoys);
        assert!(
            note.message.starts_with("Only 35 of MSFragger's 21300"),
            "{}",
            note.message
        );
        // A third decoys is an ordinary search; none is a different note.
        assert_eq!(few_decoys_note(&summary(19_577, 5_201)), None);
        assert_eq!(few_decoys_note(&summary(1_000, 0)), None);
    }

    #[test]
    fn engines_vote_with_their_strongest_score_and_ties_go_to_decoys() {
        // Spectrum 1: engine 0 says 10 (score 2), engine 1 says decoy 20
        //   (score 1): 10 wins with 2.
        // Spectrum 2: both say 10 (scores 3 and 1): the strongest vote, 3,
        //   from two engines.
        // Spectrum 3: 30 (1) against decoy 20 (1): the tie goes to the decoy.
        let input = Input {
            engines: vec!["A".into(), "B".into()],
            psms: vec![
                psm(1, 10, 0),
                psm(1, 20, 1),
                psm(2, 10, 0),
                psm(2, 10, 1),
                psm(3, 30, 0),
                psm(3, 20, 1),
            ],
            peptides: vec![peptide(10, false), peptide(20, true), peptide(30, false)],
            ..Input::default()
        };
        let peptides: HashMap<i64, &PeptideRow> =
            input.peptides.iter().map(|p| (p.id, p)).collect();
        let engine = |rows: Vec<usize>, scores: Vec<f64>| Ranked {
            q: vec![Some(0.0); rows.len()],
            pep: vec![Some(0.5); rows.len()],
            peptide_scores: scores.clone(),
            peptide_q: vec![Some(0.0); rows.len()],
            scores,
            rows,
        };
        let ranked = [
            engine(vec![0, 2, 4], vec![2.0, 3.0, 1.0]),
            engine(vec![1, 3, 5], vec![1.0, 1.0, 1.0]),
        ];
        let matches = spectrum_matches(&input, &peptides, &ranked, Pass::Spectra);
        let summary: Vec<(i64, i64, u32)> = matches
            .iter()
            .map(|m| (m.spectrum_id, m.peptide_id, m.engines))
            .collect();
        assert_eq!(summary, [(1, 10, 1), (2, 10, 2), (3, 20, 1)]);
        let scores: Vec<f64> = matches.iter().map(|m| m.score).collect();
        for (score, expected) in scores.iter().zip([2.0, 3.0, 1.0]) {
            assert!((score - expected).abs() < 1e-12, "{scores:?}");
        }
        // Best first: s2 T (FDR 1/1), s1 T (1/2), s3 D (2/2); monotone from
        // the bottom: s3 1, s1 1/2, s2 1/2.
        let q: Vec<f64> = matches.iter().map(|m| m.q_value).collect();
        assert_eq!(q, [0.5, 0.5, 1.0]);

        // Peptide 10 is best supported by spectrum 2.
        let (peptides, _, _) = peptide_results(&input, &matches, false);
        assert_eq!(peptides[0].peptide_id, 10);
        assert_eq!(peptides[0].best_spectrum_id, 2);
        assert_eq!(peptides.len(), 2);
    }

    #[test]
    fn inconsistent_input_is_refused() {
        let input = Input {
            engines: vec!["A".into()],
            psms: vec![psm(1, 99, 0)],
            ..Input::default()
        };
        let error = validate(&input, &Settings::default()).unwrap_err();
        assert_eq!(
            error.to_string(),
            "PSM 10 refers to peptide 99 that the project does not contain"
        );
    }
}
