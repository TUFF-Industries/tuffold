//! Deep-learned rescoring inputs: how well each rank-1 PSM's spectrum
//! matches the fragment intensities AlphaPeptDeep predicts for its peptide,
//! and a retention-time model that validation fine-tunes per run and fold.
//! Used whenever the models are installed and rescoring is on; when they
//! aren't, validation says so in a note.

use std::collections::{HashMap, HashSet};

use tuffold_predict::{
    DeepRetention, FEATURE_NAMES, MIN_FIT, MODELS_URL, Models, Predictions, PsmSpectrum, encode,
};
use tuffold_store::{Peaks, Project, PsmPeaks, StoreError};
use tuffold_validate::{
    EngineScore, Input, Note, NoteCode, PredictedFeatures, PsmRow, engine_scores,
};

use crate::ValidateError;

/// PSMs read from the project at a time.
const CHUNK: usize = 20_000;
/// PSMs the fragment model's settings are calibrated on.
const CALIBRATION_SAMPLE: usize = 2_000;

/// A rank-1 PSM the models can take, with peaks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Candidate {
    psm: i64,
    spectrum: i64,
    peptide: i64,
    charge: u8,
}

/// The PSMs the fragment model is calibrated on, best-scoring first.
#[derive(Debug, Clone, PartialEq, Eq)]
struct CalibrationSample {
    psms: Vec<Candidate>,
    /// How many of `psms`, from the start, are the best-scoring ones; the
    /// rest are spread over PSMs without an engine score.
    best: usize,
}

/// A random sample's fit is set by the wrong matches when few are right (a
/// pull-down's elutions: under 0.6 though the best PSMs match every b/y
/// ion), so the sample is the best-scoring PSMs instead: each engine's PSMs
/// in order of its best single score, engines interleaved by rank fraction
/// (their scores share no scale), one PSM per peptide and charge so an
/// abundant peptide's many spectra don't fill it. Targets and decoys are
/// taken alike. PSMs without an engine score (no decoys to say which score
/// is which way) fill what is left, spread evenly, so a small scored engine
/// beside a large unscored one isn't calibrated on a few dozen spectra.
fn calibration_sample(
    candidates: &[Candidate],
    scores: &HashMap<i64, EngineScore>,
) -> CalibrationSample {
    let mut by_engine: HashMap<usize, Vec<(f64, Candidate)>> = HashMap::new();
    let mut unscored = Vec::new();
    for candidate in candidates {
        match scores.get(&candidate.psm) {
            Some(score) => by_engine
                .entry(score.engine)
                .or_default()
                .push((score.score, *candidate)),
            None => unscored.push(*candidate),
        }
    }
    let mut ordered: Vec<(f64, Candidate)> = Vec::new();
    for mut psms in by_engine.into_values() {
        psms.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.psm.cmp(&b.1.psm)));
        let n = psms.len() as f64;
        ordered.extend(
            psms.into_iter()
                .enumerate()
                .map(|(rank, (_, candidate))| (rank as f64 / n, candidate)),
        );
    }
    ordered.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.psm.cmp(&b.1.psm)));
    let mut seen = HashSet::new();
    let mut psms: Vec<Candidate> = ordered
        .into_iter()
        .map(|(_, candidate)| candidate)
        .filter(|candidate| seen.insert((candidate.peptide, candidate.charge)))
        .take(CALIBRATION_SAMPLE)
        .collect();
    let best = psms.len();
    let room = CALIBRATION_SAMPLE - best;
    if let Some(stride) = unscored.len().checked_div(room) {
        psms.extend(unscored.iter().step_by(stride.max(1)).take(room));
    }
    CalibrationSample { psms, best }
}

impl CalibrationSample {
    /// What the note says the fit was measured on.
    fn describe(&self) -> String {
        let spread = self.psms.len() - self.best;
        match (self.best, spread) {
            (0, spread) => format!("{spread} rank-1 PSMs spread over the project"),
            (best, 0) => format!("the {best} best-scoring rank-1 PSMs, one per peptide and charge"),
            (best, spread) => format!(
                "the {best} best-scoring rank-1 PSMs, one per peptide and charge, and {spread} \
                 without an engine score spread over the project"
            ),
        }
    }
}

pub(crate) struct Deep {
    pub predicted: Option<PredictedFeatures>,
    pub retention: DeepRetention,
    pub notes: Vec<Note>,
}

/// Where the models live and whether they load.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelsStatus {
    pub dir: Option<std::path::PathBuf>,
    pub installed: bool,
    pub problem: Option<String>,
}

pub fn models_status() -> ModelsStatus {
    let dir = Models::default_dir();
    let (installed, problem) = match dir.as_deref().map(Models::load) {
        Some(Ok(models)) => (models.is_some(), None),
        Some(Err(error)) => (false, Some(error.to_string())),
        None => (
            false,
            Some("no home folder to keep them in; set TUFFOLD_MODELS".into()),
        ),
    };
    ModelsStatus {
        dir,
        installed,
        problem,
    }
}

/// Install AlphaPeptDeep's model archive (`pretrained_models_v3.zip`, from
/// [`tuffold_predict::MODELS_URL`]) where validation finds it.
pub fn install_models(archive: &std::path::Path) -> Result<std::path::PathBuf, String> {
    let dir = Models::default_dir().ok_or("no home folder to keep them in; set TUFFOLD_MODELS")?;
    Models::install(archive, &dir).map_err(|error| error.to_string())?;
    Ok(dir)
}

/// The installed models, if any. Missing or broken models are a note, not
/// an error (validation runs without them), so a project validated without
/// deep features says so.
pub(crate) fn installed() -> (Option<Models>, Option<Note>) {
    let missing = |place: String| {
        Note::new(
            NoteCode::DeepLearning,
            format!(
                "No deep-learning models are installed ({place}), so validation ran without \
                 predicted fragment intensities and retention times. Download {MODELS_URL} and \
                 run `tuffold models --install pretrained_models_v3.zip`."
            ),
        )
    };
    let Some(dir) = Models::default_dir() else {
        return (
            None,
            Some(missing(
                "no home folder to keep them in; set TUFFOLD_MODELS".into(),
            )),
        );
    };
    match Models::load(&dir) {
        Ok(Some(models)) => (Some(models), None),
        Ok(None) => (None, Some(missing(format!("looked in {}", dir.display())))),
        Err(error) => (
            None,
            Some(Note::new(
                NoteCode::DeepLearning,
                format!(
                    "The deep-learning models in {} could not be read ({error}); validated \
                     without them. Install them again with `tuffold models --install`.",
                    dir.display()
                ),
            )),
        ),
    }
}

pub(crate) fn prepare(
    project: &Project,
    models: &Models,
    input: &Input,
) -> Result<Deep, ValidateError> {
    let forms = project.rank1_peptides()?;
    let encoded: HashMap<i64, tuffold_predict::Encoded> = forms
        .iter()
        .filter_map(|(&id, form)| Some((id, encode(form)?)))
        .collect();
    let mut notes = Vec::new();
    let unencoded = forms.len() - encoded.len();
    if unencoded > 0 {
        notes.push(Note::new(
            NoteCode::DeepLearning,
            format!(
                "{unencoded} of {} peptides carry a residue or modification the deep-learning \
                 models don't know; their PSMs get the lowest predicted-spectrum scores and the \
                 classical retention-time error.",
                forms.len()
            ),
        ));
    }
    let candidates: Vec<Candidate> = {
        let rows: HashMap<i64, &PsmRow> = input.psms.iter().map(|psm| (psm.id, psm)).collect();
        project
            .rank1_with_peaks()?
            .into_iter()
            .filter_map(|(id, spectrum)| {
                let row = rows.get(&id)?;
                let charge = u8::try_from(row.charge).ok().filter(|&z| z > 0)?;
                encoded.contains_key(&row.peptide_id).then_some(Candidate {
                    psm: id,
                    spectrum,
                    peptide: row.peptide_id,
                    charge,
                })
            })
            .collect()
    };
    let chosen = calibration_sample(&candidates, &engine_scores(input)?);
    let spectrum_ids: Vec<i64> = chosen.psms.iter().map(|psm| psm.spectrum).collect();
    let mut peaks: HashMap<i64, Peaks> = HashMap::new();
    project.for_each_spectrum_peaks(&spectrum_ids, |spectrum, found| {
        if let Some(found) = found {
            peaks.insert(spectrum, found.clone());
        }
    })?;
    let sample: Vec<PsmPeaks> = chosen
        .psms
        .iter()
        .filter_map(|psm| {
            Some(PsmPeaks {
                psm: psm.psm,
                peptide: psm.peptide,
                charge: i32::from(psm.charge),
                peaks: peaks.get(&psm.spectrum)?.clone(),
            })
        })
        .collect();
    // (PSM, peptide, charge, index in the chunk) of the PSMs the models take.
    let spectra = |chunk: &[PsmPeaks]| -> Vec<(i64, i64, u8, usize)> {
        chunk
            .iter()
            .enumerate()
            .filter_map(|(i, psm)| {
                let charge = u8::try_from(psm.charge).ok().filter(|&z| z > 0)?;
                encoded
                    .contains_key(&psm.peptide)
                    .then_some((psm.psm, psm.peptide, charge, i))
            })
            .collect()
    };
    let calibration = spectra(&sample);
    let calibration: Vec<PsmSpectrum<'_>> = calibration
        .iter()
        .map(|&(psm, peptide, charge, i)| PsmSpectrum {
            psm,
            peptide,
            form: &forms[&peptide],
            encoded: &encoded[&peptide],
            charge,
            mz: &sample[i].peaks.mz,
            intensity: &sample[i].peaks.intensity,
        })
        .collect();
    let predicted = match models.calibrate(&calibration) {
        None => {
            notes.push(Note::new(
                NoteCode::DeepLearning,
                format!(
                    "Predicted fragment intensities were not used: AlphaPeptDeep's model doesn't \
                     describe these spectra at any collision energy or instrument it knows (fit \
                     below {MIN_FIT}), as with fragmentation other than HCD or CID."
                ),
            ));
            None
        }
        Some(scorer) => {
            notes.push(Note::new(
                NoteCode::DeepLearning,
                format!(
                    "Fragment intensities predicted by AlphaPeptDeep at NCE {} ({:?}), matched \
                     within {}, fitting the spectra with a Pearson correlation of {:.2} (90th \
                     percentile of {}).",
                    scorer.nce,
                    scorer.instrument,
                    scorer.tolerance,
                    scorer.fit,
                    chosen.describe()
                ),
            ));
            let mut values = HashMap::new();
            let mut predictions = Predictions::new();
            project.for_each_rank1_spectrum(CHUNK, |chunk| {
                let views: Vec<PsmSpectrum<'_>> = spectra(chunk)
                    .into_iter()
                    .map(|(psm, peptide, charge, i)| PsmSpectrum {
                        psm,
                        peptide,
                        form: &forms[&peptide],
                        encoded: &encoded[&peptide],
                        charge,
                        mz: &chunk[i].peaks.mz,
                        intensity: &chunk[i].peaks.intensity,
                    })
                    .collect();
                values.extend(scorer.score(&views, &mut predictions));
                Ok::<(), StoreError>(())
            })?;
            Some(PredictedFeatures {
                names: FEATURE_NAMES.iter().map(|name| name.to_string()).collect(),
                values,
            })
        }
    };
    Ok(Deep {
        predicted,
        retention: models.retention(encoded),
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(psm: i64, peptide: i64, charge: u8) -> Candidate {
        Candidate {
            psm,
            spectrum: psm,
            peptide,
            charge,
        }
    }

    fn scored(engine: usize, score: f64) -> EngineScore {
        EngineScore { engine, score }
    }

    fn ids(sample: &CalibrationSample) -> Vec<i64> {
        sample.psms.iter().map(|psm| psm.psm).collect()
    }

    #[test]
    fn the_best_scoring_psms_calibrate_one_per_peptide_and_charge() {
        // Engine 0: PSMs 1–4, scores 9, 8, 7, 1; PSM 2 repeats PSM 1's
        // peptide and charge, so it is skipped and PSM 3 comes next. Engine
        // 1: PSMs 10, 11 at 0.5, 0.2 (another scale). By rank fraction:
        // PSM 1 (0/4) and 10 (0/2) first, then 2 (1/4, skipped), then 3
        // and 11 (both 2/4, ties by id), then 4 (3/4). PSM 99 has no engine
        // score, so it comes last, filling the room left.
        let candidates = [
            candidate(1, 100, 2),
            candidate(2, 100, 2),
            candidate(3, 100, 3),
            candidate(4, 101, 2),
            candidate(10, 200, 2),
            candidate(11, 201, 2),
            candidate(99, 300, 2),
        ];
        let scores: HashMap<i64, EngineScore> = [
            (1, scored(0, 9.0)),
            (2, scored(0, 8.0)),
            (3, scored(0, 7.0)),
            (4, scored(0, 1.0)),
            (10, scored(1, 0.5)),
            (11, scored(1, 0.2)),
        ]
        .into_iter()
        .collect();
        let sample = calibration_sample(&candidates, &scores);
        assert_eq!(ids(&sample), vec![1, 10, 3, 11, 4, 99]);
        assert_eq!(sample.best, 5);
    }

    #[test]
    fn without_engine_scores_the_sample_spreads_over_the_project() {
        // 5,000 candidates, stride 5,000 / 2,000 = 2: every other PSM.
        let candidates: Vec<Candidate> = (0..5_000).map(|i| candidate(i, i, 2)).collect();
        let sample = calibration_sample(&candidates, &HashMap::new());
        assert_eq!(sample.best, 0);
        assert_eq!(sample.psms.len(), CALIBRATION_SAMPLE);
        assert_eq!(&ids(&sample)[..3], &[0, 2, 4]);
    }

    #[test]
    fn a_small_scored_engine_is_topped_up_from_the_unscored_ones() {
        // 500 scored PSMs (all distinct) leave room for 1,500; 6,000
        // unscored ones give a stride of 6,000 / 1,500 = 4: PSMs 10,000,
        // 10,004, … up to 15,996.
        let mut candidates: Vec<Candidate> = (0..500).map(|i| candidate(i, i, 2)).collect();
        candidates.extend((10_000..16_000).map(|i| candidate(i, i, 2)));
        let scores: HashMap<i64, EngineScore> =
            (0..500).map(|i| (i, scored(0, -(i as f64)))).collect();
        let sample = calibration_sample(&candidates, &scores);
        let ids = ids(&sample);
        assert_eq!(sample.best, 500);
        assert_eq!(ids.len(), CALIBRATION_SAMPLE);
        assert_eq!(&ids[..2], &[0, 1]);
        assert_eq!(&ids[500..503], &[10_000, 10_004, 10_008]);
        assert_eq!(ids[1_999], 15_996);
    }
}
