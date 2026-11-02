//! Validating a project and changing its thresholds.
//!
//! [`validate_project`] runs the whole pipeline (engines, spectra, peptides,
//! proteins) and replaces the stored results. [`apply_thresholds`] is the
//! cheap path for live threshold changes: PSM and protein FDR and the
//! minimum peptide count only change what is counted; a new peptide FDR
//! re-runs protein inference (protein groups are built from passing
//! peptides) but never the rescoring.

use std::collections::HashMap;

use tuffold_core::fragments::Tolerance;
use tuffold_store::{Project, StoreError, ValidationSummary};
use tuffold_validate::{
    LocalizationInput, Note, NoteCode, PeptideEvidence, Settings, Thresholds, ValidationError,
    ambiguous_isoform_count, infer_proteins, localize_spectrum, validate_with,
};

#[derive(Debug, thiserror::Error)]
pub enum ValidateError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error(transparent)]
    Validation(#[from] ValidationError),
    #[error("{0}")]
    Thresholds(String),
    #[error("this project has not been validated since its data last changed; validate it first")]
    NotValidated,
}

/// Whether rescoring may use the installed deep-learning models.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DeepLearning {
    /// When installed, for real (not synthetic) data.
    #[default]
    Auto,
    Off,
}

/// Validate every PSM in the project and store the results.
pub fn validate_project(
    project: &mut Project,
    settings: &Settings,
) -> Result<ValidationSummary, ValidateError> {
    validate_project_with(project, settings, DeepLearning::Auto)
}

/// [`validate_project`], choosing whether deep learning may be used.
pub fn validate_project_with(
    project: &mut Project,
    settings: &Settings,
    deep_learning: DeepLearning,
) -> Result<ValidationSummary, ValidateError> {
    let mut input = project.validation_input()?;
    // Synthetic (demo) spectra aren't real fragmentation for a model to
    // predict, and tests must not change with what a machine has installed.
    let synthetic = project
        .imports()?
        .iter()
        .any(|import| import.format == "synthetic");
    // Without peaks on any rank-1 PSM there is nothing to predict, so the
    // models aren't loaded and missing or broken ones aren't noted (and
    // results don't depend on what the machine has installed).
    let (models, broken) = if settings.rescore
        && deep_learning == DeepLearning::Auto
        && !synthetic
        && project.has_rank1_peaks()?
    {
        crate::deep::installed()
    } else {
        (None, None)
    };
    let deep = models
        .as_ref()
        .map(|models| crate::deep::prepare(project, models, &input))
        .transpose()?;
    if let Some(deep) = &deep {
        input.predicted = deep.predicted.clone();
    }
    let retention = deep
        .as_ref()
        .map(|deep| &deep.retention as &dyn tuffold_validate::RetentionPredictor);
    let mut validation = validate_with(&input, settings, retention)?;
    let (localizations, resolution) = localize_project_spectra(project, &input, &validation)?;
    validation.localizations = localizations;
    validation
        .notes
        .extend(resolution.map(|check| check.note()));
    validation.notes.extend(broken);
    validation
        .notes
        .extend(deep.into_iter().flat_map(|deep| deep.notes));
    project.save_validation(settings, &validation)?;
    current(project)
}

/// Fragment tolerances for PTM localization: high-resolution MS/MS, and the
/// ion-trap fallback (as for the deep features' fragment matching).
pub const HIGH_RESOLUTION: Tolerance = Tolerance::Ppm(20.0);
pub const LOW_RESOLUTION: Tolerance = Tolerance::Da(0.5);
/// Spectra (the first by id) whose matches decide the tolerance.
const RESOLUTION_SAMPLE: usize = 2_000;

/// 20 ppm unless it finds under half the b/y ions 0.5 Da does: on
/// high-resolution spectra nearly every 0.5 Da match is also a 20 ppm one
/// (plus a few chance peaks), on ion-trap spectra few are.
fn localization_tolerance(fine: usize, wide: usize) -> Tolerance {
    if fine * 2 < wide {
        LOW_RESOLUTION
    } else {
        HIGH_RESOLUTION
    }
}

/// The matches that chose localization's fragment tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
struct ResolutionCheck {
    /// Spectra with peaks whose reported isoform's b/y ions were counted.
    spectra: usize,
    /// Ions matched at [`HIGH_RESOLUTION`] and at [`LOW_RESOLUTION`].
    fine: usize,
    wide: usize,
}

impl ResolutionCheck {
    fn tolerance(self) -> Tolerance {
        localization_tolerance(self.fine, self.wide)
    }

    /// Says which tolerance was used and why, in the validation notes.
    fn note(self) -> Note {
        let Self {
            spectra,
            fine,
            wide,
        } = self;
        let message = if wide == 0 {
            // Nothing matched either way: 20 ppm is a default, not a finding.
            format!(
                "PTM sites were localized with fragments matched within 20 ppm, the default: no \
                 reported isoform's b/y ion matched at 20 ppm or at 0.5 Da in {spectra} sampled \
                 spectra."
            )
        } else if self.tolerance() == LOW_RESOLUTION {
            format!(
                "PTM sites were localized with fragments matched within 0.5 Da: at 20 ppm the \
                 reported isoforms' b/y ions matched {fine} times in {spectra} sampled spectra, \
                 against {wide} at 0.5 Da (low-resolution MS/MS)."
            )
        } else {
            format!(
                "PTM sites were localized with fragments matched within 20 ppm: the reported \
                 isoforms' b/y ions matched {fine} times at 20 ppm in {spectra} sampled spectra, \
                 against {wide} at 0.5 Da (high-resolution MS/MS)."
            )
        };
        Note::new(NoteCode::Localization, message)
    }
}

/// Theoretical ions with a peak within `tolerance` (peaks in any order:
/// `.sf3` peaks come in two parts).
fn matched_ions(
    ions: &[tuffold_core::fragments::FragmentIon],
    mz: &[f64],
    tolerance: Tolerance,
) -> usize {
    ions.iter()
        .filter(|ion| {
            let window = tolerance.window(ion.mz);
            mz.iter().any(|&m| (m - ion.mz).abs() <= window)
        })
        .count()
}

/// Spectra localized, and the check that chose the fragment tolerance
/// (`None` when no ambiguous spectrum had peaks to decide from).
type Localized = (
    Vec<(i64, tuffold_validate::SpectrumLocalization)>,
    Option<ResolutionCheck>,
);

fn localize_project_spectra(
    project: &Project,
    input: &tuffold_validate::Input,
    validation: &tuffold_validate::Validation,
) -> Result<Localized, StoreError> {
    if !input
        .peptides
        .iter()
        .any(|peptide| !peptide.is_decoy && peptide.modifications > 0)
    {
        return Ok((Vec::new(), None));
    }
    let forms = project.rank1_peptides()?;
    let ambiguous_peptides: HashMap<i64, &tuffold_core::PeptideForm> = forms
        .iter()
        .filter(|(_, form)| ambiguous_isoform_count(form).is_some())
        .map(|(&id, form)| (id, form))
        .collect();
    if ambiguous_peptides.is_empty() {
        return Ok((Vec::new(), None));
    }
    let mut psm_charges: HashMap<(i64, i64), Option<u8>> = HashMap::new();
    for psm in &input.psms {
        if ambiguous_peptides.contains_key(&psm.peptide_id) {
            let charge = u8::try_from(psm.charge).ok().filter(|&c| c > 0);
            psm_charges
                .entry((psm.spectrum_id, psm.peptide_id))
                .or_insert(charge);
        }
    }
    let mut targets: Vec<(i64, i64)> = validation
        .spectra
        .iter()
        .filter(|s| !s.is_decoy && ambiguous_peptides.contains_key(&s.peptide_id))
        .map(|s| (s.spectrum_id, s.peptide_id))
        .collect();
    targets.sort_by_key(|&(spectrum_id, _)| spectrum_id);
    let spectrum_ids: Vec<i64> = targets.iter().map(|&(id, _)| id).collect();
    let by_spectrum: HashMap<i64, i64> = targets.into_iter().collect();
    let charge_of = |spectrum_id: i64, peptide_id: i64| {
        psm_charges
            .get(&(spectrum_id, peptide_id))
            .copied()
            .flatten()
    };

    // Low-resolution (ion-trap) MS/MS barely matches at 20 ppm; decide from
    // the reported isoforms' b/y ions in a deterministic sample of spectra.
    let sample = &spectrum_ids[..spectrum_ids.len().min(RESOLUTION_SAMPLE)];
    let mut check = ResolutionCheck {
        spectra: 0,
        fine: 0,
        wide: 0,
    };
    project.for_each_spectrum_peaks(sample, |spectrum_id, peaks| {
        let (Some(peaks), Some(&peptide_id)) = (peaks, by_spectrum.get(&spectrum_id)) else {
            return;
        };
        let Some(&form) = ambiguous_peptides.get(&peptide_id) else {
            return;
        };
        let Some(ions) =
            tuffold_core::fragments::fragment_ions(form, charge_of(spectrum_id, peptide_id))
        else {
            return;
        };
        check.spectra += 1;
        check.fine += matched_ions(&ions, &peaks.mz, HIGH_RESOLUTION);
        check.wide += matched_ions(&ions, &peaks.mz, LOW_RESOLUTION);
    })?;
    let tolerance = check.tolerance();

    let mut localizations = Vec::with_capacity(spectrum_ids.len());
    project.for_each_spectrum_peaks(&spectrum_ids, |spectrum_id, peaks| {
        let Some(&peptide_id) = by_spectrum.get(&spectrum_id) else {
            return;
        };
        let Some(&form) = ambiguous_peptides.get(&peptide_id) else {
            return;
        };
        let charge = charge_of(spectrum_id, peptide_id);
        let (mz, intensity) = match peaks {
            Some(peaks) => (peaks.mz.as_slice(), peaks.intensity.as_slice()),
            None => (&[][..], &[][..]),
        };
        let loc = localize_spectrum(
            &LocalizationInput {
                spectrum_id,
                form,
                charge,
                mz,
                intensity,
            },
            tolerance,
        );
        if !matches!(loc, tuffold_validate::SpectrumLocalization::Unambiguous) {
            localizations.push((spectrum_id, loc));
        }
    })?;
    Ok((localizations, (check.spectra > 0).then_some(check)))
}

/// Apply new thresholds to existing results.
pub fn apply_thresholds(
    project: &mut Project,
    thresholds: &Thresholds,
) -> Result<ValidationSummary, ValidateError> {
    thresholds.check().map_err(ValidateError::Thresholds)?;
    let summary = current(project)?;
    // Protein groups are built from the passing peptides, so they only need
    // rebuilding when the peptide cut changes (its value or its mode).
    let before = summary.settings.thresholds;
    if thresholds.mode == before.mode && thresholds.peptide == before.peptide {
        project.save_thresholds(thresholds)?;
    } else {
        let rows = project.inference_rows()?;
        let evidence: Vec<PeptideEvidence<'_>> = rows
            .peptides
            .iter()
            .map(|peptide| PeptideEvidence {
                id: peptide.id,
                sequence: &peptide.sequence,
                is_decoy: peptide.is_decoy,
                score: peptide.score,
                q_value: peptide.q_value,
                pep: peptide.pep,
            })
            .collect();
        let inference = infer_proteins(
            &evidence,
            &rows.peptide_proteins,
            &rows.proteins,
            thresholds,
        );
        project.save_inference(thresholds, &inference)?;
    }
    current(project)
}

/// The stored summary, refusing stale or missing results.
fn current(project: &Project) -> Result<ValidationSummary, ValidateError> {
    match project.validation_summary()? {
        Some(summary) if summary.validated_at.is_some() => Ok(summary),
        _ => Err(ValidateError::NotValidated),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerance_widens_only_when_20_ppm_finds_under_half() {
        // Orbitrap-like: 0.5 Da adds a few chance peaks to 20 ppm's matches.
        assert_eq!(localization_tolerance(900, 1_000), HIGH_RESOLUTION);
        assert_eq!(localization_tolerance(500, 1_000), HIGH_RESOLUTION);
        // Ion trap: 20 ppm catches only what happens to land near.
        assert_eq!(localization_tolerance(499, 1_000), LOW_RESOLUTION);
        assert_eq!(localization_tolerance(40, 1_000), LOW_RESOLUTION);
        // Nothing to go by: stay at 20 ppm.
        assert_eq!(localization_tolerance(0, 0), HIGH_RESOLUTION);
    }

    #[test]
    fn the_note_names_the_tolerance_and_the_counts_behind_it() {
        // 1,405 of 17,519 is under half: the ion trap's 0.5 Da.
        let trap = ResolutionCheck {
            spectra: 2_000,
            fine: 1_405,
            wide: 17_519,
        };
        let note = trap.note();
        assert_eq!(note.code, NoteCode::Localization);
        assert_eq!(
            note.message,
            "PTM sites were localized with fragments matched within 0.5 Da: at 20 ppm the \
             reported isoforms' b/y ions matched 1405 times in 2000 sampled spectra, against \
             17519 at 0.5 Da (low-resolution MS/MS)."
        );
        // 900 of 1,000 is over half: 20 ppm.
        let orbitrap = ResolutionCheck {
            spectra: 37,
            fine: 900,
            wide: 1_000,
        };
        assert_eq!(
            orbitrap.note().message,
            "PTM sites were localized with fragments matched within 20 ppm: the reported \
             isoforms' b/y ions matched 900 times at 20 ppm in 37 sampled spectra, against 1000 \
             at 0.5 Da (high-resolution MS/MS)."
        );
        // No match either way: 20 ppm by default, and the note says so.
        let blank = ResolutionCheck {
            spectra: 3,
            fine: 0,
            wide: 0,
        };
        assert_eq!(
            blank.note().message,
            "PTM sites were localized with fragments matched within 20 ppm, the default: no \
             reported isoform's b/y ion matched at 20 ppm or at 0.5 Da in 3 sampled spectra."
        );
    }

    #[test]
    fn ions_match_within_their_window_in_any_peak_order() {
        use tuffold_core::fragments::{FragmentIon, IonSeries};
        let ion = |mz| FragmentIon {
            series: IonSeries::B,
            ordinal: 2,
            charge: 1,
            mz,
        };
        let ions = [ion(500.0), ion(800.0)];
        // 500.004 is 8 ppm off (window 0.010 Da); 800.3 is 375 ppm off.
        let peaks = [800.3, 120.0, 500.004];
        assert_eq!(matched_ions(&ions, &peaks, HIGH_RESOLUTION), 1);
        assert_eq!(matched_ions(&ions, &peaks, LOW_RESOLUTION), 2);
    }
}
