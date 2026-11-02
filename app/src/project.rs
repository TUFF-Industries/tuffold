//! Project view model: what the app shows for an open `.tuff` file, plus the
//! create / open / add-data / threshold operations behind it. Pure data and
//! engine calls; no GPUI, so everything here is unit-tested and runs off the
//! UI thread.
//!
//! The app validates for the user after every import (with the project's
//! saved settings). Opening a project never writes to it: missing or stale
//! results show as [`ValidationState::Pending`] until the user validates.

use std::path::{Path, PathBuf};

use tuffold_engine::{
    Outcome, OutcomeKind, SearchReport, SearchSettings, apply_thresholds, create_demo_project,
    duplicate_note, import_inputs, search_files, validate_project,
};
use tuffold_store::{
    ImportRecord, LevelCounts, Project, ProjectSummary, StoreError, ValidationSummary,
};
use tuffold_validate::{ThresholdMode, ThresholdRequest, Thresholds, format_rate};

pub const EXTENSION: &str = "tuff";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FindingLevel {
    Error,
    Warning,
}

/// Something the user should know about the last action, tied to an input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Finding {
    pub level: FindingLevel,
    /// The input's file name.
    pub source: String,
    /// Stable issue code, when the finding has one.
    pub code: Option<String>,
    pub message: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusKind {
    Ready,
    Limited,
    Invalid,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Status {
    pub kind: StatusKind,
    pub label: String,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidationState {
    /// No PSMs yet, so nothing to validate.
    NoData,
    /// Data without current results (never validated, or imported since).
    Pending { stale: bool },
    /// Boxed: a summary is far larger than the other states.
    Done(Box<ValidationSummary>),
    /// Validation could not run; the imported data is still saved.
    Failed(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ProjectReport {
    pub path: PathBuf,
    pub summary: ProjectSummary,
    pub imports: Vec<ImportRecord>,
    /// Findings from the most recent create / add-data action.
    pub findings: Vec<Finding>,
    pub validation: ValidationState,
    /// A one-line note about the action just taken (e.g. what was undone).
    pub notice: Option<String>,
}

/// FDR thresholds offered as one-click choices.
pub const FDR_CHOICES: [f64; 3] = [0.001, 0.01, 0.05];
/// Probability thresholds offered as one-click choices.
pub const PROBABILITY_CHOICES: [f64; 3] = [0.9, 0.95, 0.99];

/// The one-click choices for a threshold mode.
pub fn choices(mode: ThresholdMode) -> [f64; 3] {
    match mode {
        ThresholdMode::Fdr => FDR_CHOICES,
        ThresholdMode::Probability => PROBABILITY_CHOICES,
    }
}

/// `thresholds` switched to `mode`, at that mode's defaults (an FDR is no
/// probability), keeping min peptides; unchanged if already in `mode`.
pub fn switched(thresholds: Thresholds, mode: ThresholdMode) -> Thresholds {
    thresholds
        .adjusted(&ThresholdRequest {
            by: Some(mode),
            ..ThresholdRequest::default()
        })
        .unwrap_or(thresholds)
}
/// Minimum peptides per protein group offered as choices.
pub const MIN_PEPTIDE_CHOICES: [u32; 2] = [1, 2];

/// A level validation reports on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Level {
    Spectra,
    Peptides,
    ProteinGroups,
}

impl Level {
    pub const ALL: [Self; 3] = [Self::Spectra, Self::Peptides, Self::ProteinGroups];

    pub fn label(self) -> &'static str {
        match self {
            Self::Spectra => "Spectra",
            Self::Peptides => "Peptides",
            Self::ProteinGroups => "Protein groups",
        }
    }

    fn validated(self) -> tuffold_validate::Level {
        match self {
            Self::Spectra => tuffold_validate::Level::Psm,
            Self::Peptides => tuffold_validate::Level::Peptide,
            Self::ProteinGroups => tuffold_validate::Level::Protein,
        }
    }

    /// This level's threshold: an FDR or a probability, by mode.
    pub fn threshold(self, thresholds: &Thresholds) -> f64 {
        thresholds.value(self.validated())
    }

    /// `thresholds` with this level's threshold replaced.
    pub fn with(self, thresholds: &Thresholds, value: f64) -> Thresholds {
        thresholds.with(self.validated(), value)
    }

    fn counts(self, summary: &ValidationSummary) -> &LevelCounts {
        match self {
            Self::Spectra => &summary.spectra,
            Self::Peptides => &summary.peptides,
            Self::ProteinGroups => &summary.protein_groups,
        }
    }
}

/// One level's headline: how many targets pass, and what that means.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LevelView {
    pub level: Level,
    pub targets: usize,
    pub caption: String,
    /// Why nothing (or less than expected) passes, when the threshold is
    /// below what the target-decoy estimate can resolve.
    pub hint: Option<String>,
}

pub fn level_views(summary: &ValidationSummary) -> Vec<LevelView> {
    let thresholds = &summary.settings.thresholds;
    Level::ALL
        .into_iter()
        .map(|level| {
            let counts = level.counts(summary);
            let threshold = level.threshold(thresholds);
            let mut caption = match thresholds.mode {
                ThresholdMode::Fdr => format!(
                    "at {} FDR · {} passing",
                    format_rate(threshold),
                    plural(counts.decoys, "decoy", "decoys")
                ),
                // A probability says nothing about the error rate of the
                // whole list; the decoys that pass do.
                ThresholdMode::Probability => format!(
                    "at ≥{} probability · est. FDR {} ({})",
                    format_rate(threshold),
                    counts
                        .estimated_fdr()
                        .map_or_else(|| "—".to_string(), format_rate),
                    plural(counts.decoys, "decoy", "decoys")
                ),
            };
            if level == Level::ProteinGroups && thresholds.min_peptides > 1 {
                caption.push_str(&format!(" · ≥{} peptides", thresholds.min_peptides));
            }
            let hint = counts
                .lowest_q
                .filter(|lowest| thresholds.mode == ThresholdMode::Fdr && *lowest > threshold)
                .map(|lowest| {
                    format!(
                        "None can pass {}: with this few {}, the lowest q-value reached is {}.",
                        format_rate(threshold),
                        level.label().to_lowercase(),
                        format_rate(lowest)
                    )
                });
            LevelView {
                level,
                targets: counts.targets,
                caption,
                hint,
            }
        })
        .collect()
}

fn plural(count: usize, one: &str, many: &str) -> String {
    format!("{count} {}", if count == 1 { one } else { many })
}

impl ProjectReport {
    pub fn title(&self) -> String {
        self.path.file_stem().map_or_else(
            || self.path.display().to_string(),
            |stem| stem.to_string_lossy().into_owned(),
        )
    }

    pub fn status(&self) -> Status {
        let count = |level| self.findings.iter().filter(|f| f.level == level).count();
        let (errors, warnings) = (count(FindingLevel::Error), count(FindingLevel::Warning));
        if errors > 0 {
            Status {
                kind: StatusKind::Invalid,
                label: format!("{} not imported", plural(errors, "input", "inputs")),
            }
        } else if warnings > 0 {
            Status {
                kind: StatusKind::Limited,
                label: format!(
                    "Saved with {}",
                    plural(warnings, "limitation", "limitations")
                ),
            }
        } else {
            Status {
                kind: StatusKind::Ready,
                label: "Saved".into(),
            }
        }
    }
}

/// Where a new project goes: `.tuff` is appended when missing.
pub fn project_path(chosen: &Path) -> PathBuf {
    let has_extension = chosen
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case(EXTENSION));
    if has_extension {
        chosen.to_path_buf()
    } else {
        let mut name = chosen.as_os_str().to_owned();
        name.push(".");
        name.push(EXTENSION);
        PathBuf::from(name)
    }
}

/// Create a project at `path` (replacing an existing *project* the user chose
/// to overwrite, never any other file) and import `inputs` into it.
pub fn create(path: &Path, inputs: &[PathBuf]) -> Result<ProjectReport, String> {
    let path = free_project_path(path)?;
    let mut project = Project::create(&path).map_err(|error| error.to_string())?;
    let outcomes = import_inputs(&mut project, inputs, None);
    let validation = ensure_validated(&mut project);
    report(&path, &project, &outcomes, validation)
}

/// Create the synthetic demo project at `path` (same replacement rules as
/// [`create`]) and open it.
pub fn create_demo(path: &Path) -> Result<ProjectReport, String> {
    let path = free_project_path(path)?;
    create_demo_project(&path).map_err(|error| error.to_string())?;
    open(&path)
}

/// `path` with `.tuff`, cleared for a new project: an existing *project* the
/// user chose to overwrite in the save dialog is removed; any other existing
/// file is refused.
fn free_project_path(chosen: &Path) -> Result<PathBuf, String> {
    let path = project_path(chosen);
    if path.exists() {
        Project::open(&path).map_err(|error| match error {
            StoreError::NotAProject(_) => format!(
                "{} exists and is not a Tuffold project; choose another name",
                path.display()
            ),
            other => other.to_string(),
        })?;
        std::fs::remove_file(&path).map_err(|error| format!("{}: {error}", path.display()))?;
    }
    Ok(path)
}

/// Import more inputs into an existing project.
pub fn add(path: &Path, inputs: &[PathBuf]) -> Result<ProjectReport, String> {
    let mut project = Project::open(path).map_err(|error| error.to_string())?;
    let outcomes = import_inputs(&mut project, inputs, None);
    let validation = ensure_validated(&mut project);
    report(path, &project, &outcomes, validation)
}

/// Search spectra files against a FASTA with the built-in search, import
/// the PSMs, and validate. Problems with single files become findings; a
/// FASTA or settings problem fails the whole search.
pub fn search(path: &Path, fasta: &Path, spectra: &[PathBuf]) -> Result<ProjectReport, String> {
    let mut project = Project::open(path).map_err(|error| error.to_string())?;
    let searched = search_files(&mut project, fasta, spectra, SearchSettings::default())
        .map_err(|error| error.to_string())?;
    let validation = ensure_validated(&mut project);
    let mut report = report(path, &project, &[], validation)?;
    report.findings = search_findings(fasta, &searched);
    Ok(report)
}

fn search_findings(fasta: &Path, searched: &SearchReport) -> Vec<Finding> {
    let file_name = |path: &Path| {
        path.file_name().map_or_else(
            || path.display().to_string(),
            |name| name.to_string_lossy().into_owned(),
        )
    };
    let decoys = (searched.dropped_decoys > 0).then(|| Finding {
        level: FindingLevel::Warning,
        source: file_name(fasta),
        code: Some("fasta-decoys-replaced".into()),
        message: format!(
            "{} decoy entries in the FASTA were left out; the search generates its own \
             (rev_ prefix).",
            searched.dropped_decoys
        ),
    });
    let excluded = (!searched.excluded_targets.is_empty()).then(|| Finding {
        level: FindingLevel::Warning,
        source: file_name(fasta),
        code: Some("fasta-entries-left-out".into()),
        message: format!(
            "{} entries were left out (accession contains rev_ or repeats): {}.",
            searched.excluded_targets.len(),
            searched.excluded_targets.join(", ")
        ),
    });
    let runs = searched.runs.iter().flat_map(|run| {
        let source = file_name(&run.input);
        let done = match &run.outcome {
            Err(message) => {
                return vec![Finding {
                    level: FindingLevel::Error,
                    source,
                    code: None,
                    message: message.clone(),
                }];
            }
            Ok(done) => done,
        };
        let skipped = done.skipped.describe().map(|skipped| Finding {
            level: FindingLevel::Warning,
            source: source.clone(),
            code: Some("spectra-not-searched".into()),
            message: format!(
                "{} MS/MS spectra searched; not searched: {skipped}.",
                done.searched
            ),
        });
        // Not a problem, but the user should know their peaks were picked here.
        let centroided = (done.centroided > 0).then(|| Finding {
            level: FindingLevel::Warning,
            source: source.clone(),
            code: Some("profile-centroided".into()),
            message: format!(
                "{} profile-mode spectra were centroided before searching (the stored peaks \
                 stay as recorded).",
                done.centroided
            ),
        });
        let unmatched = done.imported.is_none().then(|| Finding {
            level: FindingLevel::Warning,
            source,
            code: Some("no-matches".into()),
            message: "No spectrum matched a peptide, so nothing was imported (check the FASTA \
                      and the spectra)."
                .into(),
        });
        skipped
            .into_iter()
            .chain(centroided)
            .chain(unmatched)
            .collect()
    });
    decoys.into_iter().chain(excluded).chain(runs).collect()
}

/// Undo the last import, search, or sample-group change (or redo it, when
/// it was just undone), and open the project as it now is.
pub fn undo(path: &Path) -> Result<ProjectReport, String> {
    let note = tuffold_engine::undo(path).map_err(|error| error.to_string())?;
    let notice = if note.redo {
        format!("Redid “{}”. Undo again to undo it.", note.label)
    } else {
        format!("Undid “{}”. Undo again to redo it.", note.label)
    };
    Ok(ProjectReport {
        notice: Some(notice),
        ..open(path)?
    })
}

/// Open a project as it is; nothing is written.
pub fn open(path: &Path) -> Result<ProjectReport, String> {
    let project = Project::open(path).map_err(|error| error.to_string())?;
    let validation = stored_validation(&project).unwrap_or_else(ValidationState::Failed);
    report(path, &project, &[], validation)
}

/// Validate (or re-validate) with the saved settings.
pub fn validate(path: &Path) -> Result<ProjectReport, String> {
    let mut project = Project::open(path).map_err(|error| error.to_string())?;
    let validation = ensure_validated(&mut project);
    report(path, &project, &[], validation)
}

/// Apply new thresholds to a validated project.
/// Star or unstar a protein (by accession) in the project.
pub fn set_starred(path: &Path, accession: &str, starred: bool) -> Result<(), String> {
    tuffold_store::Project::open(path)
        .and_then(|mut project| project.set_starred(accession, starred))
        .map_err(|error| error.to_string())
}

pub fn set_thresholds(path: &Path, thresholds: &Thresholds) -> Result<ProjectReport, String> {
    let mut project = Project::open(path).map_err(|error| error.to_string())?;
    let summary = apply_thresholds(&mut project, thresholds).map_err(|error| error.to_string())?;
    report(
        path,
        &project,
        &[],
        ValidationState::Done(Box::new(summary)),
    )
}

/// Validate with the saved settings when results are missing or stale.
fn ensure_validated(project: &mut Project) -> ValidationState {
    match stored_validation(project) {
        Ok(ValidationState::Pending { .. }) => {
            let settings = project
                .validation_summary()
                .ok()
                .flatten()
                .map(|summary| summary.settings)
                .unwrap_or_default();
            validate_project(project, &settings)
                .map(|summary| ValidationState::Done(Box::new(summary)))
                .unwrap_or_else(|error| ValidationState::Failed(error.to_string()))
        }
        Ok(state) => state,
        Err(message) => ValidationState::Failed(message),
    }
}

/// What the project holds now, without validating.
fn stored_validation(project: &Project) -> Result<ValidationState, String> {
    if project.summary().map_err(|e| e.to_string())?.psms == 0 {
        return Ok(ValidationState::NoData);
    }
    Ok(
        match project.validation_summary().map_err(|e| e.to_string())? {
            Some(summary) if summary.validated_at.is_some() => {
                ValidationState::Done(Box::new(summary))
            }
            Some(_) => ValidationState::Pending { stale: true },
            None => ValidationState::Pending { stale: false },
        },
    )
}

fn report(
    path: &Path,
    project: &Project,
    outcomes: &[Outcome],
    validation: ValidationState,
) -> Result<ProjectReport, String> {
    Ok(ProjectReport {
        path: path.to_path_buf(),
        summary: project.summary().map_err(|error| error.to_string())?,
        imports: project.imports().map_err(|error| error.to_string())?,
        findings: outcomes.iter().flat_map(findings).collect(),
        validation,
        notice: None,
    })
}

fn findings(outcome: &Outcome) -> Vec<Finding> {
    let source = outcome.input.file_name().map_or_else(
        || outcome.input.display().to_string(),
        |name| name.to_string_lossy().into_owned(),
    );
    match &outcome.kind {
        OutcomeKind::Failed(message) => vec![Finding {
            level: FindingLevel::Error,
            source,
            code: None,
            message: message.clone(),
        }],
        OutcomeKind::Imported {
            warnings,
            report,
            notes,
            ..
        } => notes
            .iter()
            .map(|note| Finding {
                level: FindingLevel::Warning,
                source: source.clone(),
                code: Some("import-note".into()),
                message: note.clone(),
            })
            .chain(warnings.iter().map(|issue| Finding {
                level: FindingLevel::Warning,
                source: source.clone(),
                code: Some(issue.code.as_str().to_string()),
                message: issue.message.clone(),
            }))
            .chain(duplicate_note(report).map(|message| Finding {
                level: FindingLevel::Warning,
                source: source.clone(),
                code: Some("duplicates-kept-once".into()),
                message,
            }))
            .collect(),
        OutcomeKind::Attached(report) => {
            let mut found = Vec::new();
            if report.identified_without_peaks > 0 {
                found.push(Finding {
                    level: FindingLevel::Warning,
                    source: source.clone(),
                    code: Some("spectra-missing".into()),
                    message: format!(
                        "{} of run {}'s identified spectra are not in this file, so they have \
                         no peaks; their PSMs stay loaded",
                        plural(report.identified_without_peaks, "spectrum", "spectra"),
                        report.run
                    ),
                });
            }
            if report.ms3_scans > 0 && report.ms3_attached == 0 {
                found.push(Finding {
                    level: FindingLevel::Warning,
                    source,
                    code: Some("ms3-unlinked".into()),
                    message: format!(
                        "run {}'s {} MS3 scans name no identified MS/MS spectrum, so TMT \
                         reporters in them can't be read",
                        report.run, report.ms3_scans
                    ),
                });
            }
            found
        }
    }
}

/// Add GO slim annotations from the picked files; a sentence on what they
/// covered.
pub fn add_go(path: &Path, inputs: &[PathBuf]) -> Result<String, String> {
    let files = tuffold_engine::GoFiles::classify(inputs)?;
    let mut project = Project::open(path).map_err(|e| e.to_string())?;
    let report = tuffold_engine::annotate_go(&mut project, &files).map_err(|e| e.to_string())?;
    Ok(format!(
        "GO: {} of {} proteins with a UniProt accession reach {} slim terms; filter by them \
         in Proteins and Samples.",
        report.annotated, report.uniprot_proteins, report.terms
    ))
}

/// Write the project's results as mzIdentML 1.2 to `target`; a sentence on
/// what it holds.
pub fn export_mzid(source: &Path, target: &Path) -> Result<String, String> {
    let project = Project::open(source).map_err(|e| e.to_string())?;
    let file = std::fs::File::create(target).map_err(|e| format!("{}: {e}", target.display()))?;
    let mut out = std::io::BufWriter::new(file);
    let written = tuffold_engine::write_mzid(&project, &mut out)
        .map_err(|e| e.to_string())
        .and_then(|report| {
            std::io::Write::flush(&mut out).map_err(|e| e.to_string())?;
            Ok(report)
        });
    match written {
        Ok(report) => Ok(format!(
            "Exported {}: {} spectra ({} passing), {} protein groups ({} passing).",
            target.display(),
            report.spectra,
            report.passing_spectra,
            report.protein_groups,
            report.passing_groups
        )),
        Err(message) => {
            drop(out);
            // Best effort: don't leave half a file behind.
            let _ = std::fs::remove_file(target);
            Err(message)
        }
    }
}

/// Write the project's spectral library (.msp) to `target`; a sentence on
/// what was written and left out.
pub fn export_library(source: &Path, target: &Path) -> Result<String, String> {
    let project = Project::open(source).map_err(|e| e.to_string())?;
    let file = std::fs::File::create(target).map_err(|e| format!("{}: {e}", target.display()))?;
    let mut out = std::io::BufWriter::new(file);
    let written = tuffold_engine::write_msp(&project, &mut out).and_then(|report| {
        std::io::Write::flush(&mut out)?;
        Ok(report)
    });
    let report = match written {
        Ok(report) => report,
        Err(error) => {
            drop(out);
            // Best effort: don't leave half a library behind.
            let _ = std::fs::remove_file(target);
            return Err(error.to_string());
        }
    };
    let mut message = format!(
        "Exported {}: {} peptide ions, each its best spectrum.",
        target.display(),
        report.entries
    );
    if report.without_peaks > 0 {
        message.push_str(&format!(
            " {} passing peptide ions have no spectrum with peaks (add the spectra files).",
            report.without_peaks
        ));
    }
    if report.unnamed_modifications > 0 {
        message.push_str(&format!(
            " {} carry a modification without a Unimod name and were left out.",
            report.unnamed_modifications
        ));
    }
    Ok(message)
}

/// Write the project's Excel workbook to `target`; a sentence on what was
/// written (and anything cut at Excel's row limit).
pub fn export_workbook(source: &Path, target: &Path) -> Result<String, String> {
    let project = Project::open(source).map_err(|e| e.to_string())?;
    let written = tuffold_engine::write_workbook(&project, target).map_err(|e| e.to_string())?;
    let sheets: Vec<String> = written
        .sheets
        .iter()
        .map(|(name, rows)| format!("{name} ({rows})"))
        .collect();
    let mut message = format!("Exported {}: {}.", target.display(), sheets.join(", "));
    for (name, total) in &written.cut {
        message.push_str(&format!(
            " {name} has {total} rows but Excel holds {}; the sheet has the first ones.",
            tuffold_engine::MAX_ROWS
        ));
    }
    Ok(message)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixtures() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures")
    }

    fn all_inputs() -> Vec<PathBuf> {
        vec![
            fixtures().join("fragpipe_synthetic"),
            fixtures().join("third_party/pyteomics/test.pep.xml"),
            fixtures().join("third_party/pyteomics/test.mzid"),
            fixtures().join("third_party/pyteomics/test.mgf"),
        ]
    }

    #[test]
    fn undo_reverses_the_last_import_and_says_so() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("study.tuff");
        let created = create(&path, &all_inputs()).unwrap();
        let path = created.path.clone();
        let added = add(&path, &[fixtures().join("fragpipe_synthetic/protein.fas")]).unwrap();
        // A FASTA alone isn't importable: the create stays the undo point.
        assert!(
            added
                .findings
                .iter()
                .any(|f| f.level == FindingLevel::Error)
        );
        let undone = undo(&path).unwrap();
        assert!(undone.imports.is_empty(), "{:?}", undone.imports);
        let notice = undone.notice.unwrap();
        assert!(notice.starts_with("Undid “import 4 inputs”"), "{notice}");
        let redone = undo(&path).unwrap();
        assert_eq!(redone.imports.len(), created.imports.len());
        assert!(redone.notice.unwrap().starts_with("Redid"));
        // Anything else afterwards carries no notice.
        assert_eq!(open(&path).unwrap().notice, None);
    }

    #[test]
    fn search_imports_and_reports_what_it_left_out() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("searched.tuff");
        create(&path, &[]).unwrap();
        let fasta = fixtures().join("fragpipe_synthetic/protein.fas");
        // test.mgf: 2 spectra of 6 peaks (too few to score); protein.fas
        // holds one rev_ decoy.
        let spectra = [fixtures().join("third_party/pyteomics/test.mgf")];
        let report = search(&path, &fasta, &spectra).unwrap();
        let codes: Vec<(&str, Option<&str>)> = report
            .findings
            .iter()
            .map(|f| (f.source.as_str(), f.code.as_deref()))
            .collect();
        assert_eq!(
            codes,
            [
                ("protein.fas", Some("fasta-decoys-replaced")),
                ("test.mgf", Some("spectra-not-searched")),
                ("test.mgf", Some("no-matches")),
            ]
        );
        assert_eq!(
            report.findings[1].message,
            "0 MS/MS spectra searched; not searched: 2 with fewer than 15 peaks."
        );
        // Nothing matched, so nothing was imported (and it can be redone).
        assert!(report.imports.is_empty());

        // A missing spectra file fails on its own.
        let gone = search(&path, &fasta, &[dir.path().join("gone.mgf")]).unwrap();
        assert_eq!(gone.findings[1].level, FindingLevel::Error);
        let missing = search(&path, &dir.path().join("none.fasta"), &spectra).unwrap_err();
        assert!(missing.contains("none.fasta"), "{missing}");
    }

    #[test]
    fn create_imports_everything_and_reports_the_failure() {
        let dir = tempfile::tempdir().unwrap();

        let report = create(&dir.path().join("study"), &all_inputs()).unwrap();

        // Hand-derived totals (see tuffold-cli tests): 33 PSMs, 7 decoys.
        assert_eq!(report.path, dir.path().join("study.tuff"));
        assert_eq!(report.title(), "study");
        assert_eq!((report.summary.psms, report.summary.decoy_psms), (33, 7));
        assert_eq!(report.imports.len(), 3);
        // test.mgf names run "test", which no import created; the FragPipe
        // bundle's psm.tsv is flagged as already filtered.
        let findings: Vec<_> = report
            .findings
            .iter()
            .map(|f| (f.level, f.source.as_str()))
            .collect();
        assert_eq!(
            findings,
            [
                (FindingLevel::Warning, "fragpipe_synthetic"),
                (FindingLevel::Error, "test.mgf")
            ]
        );
        assert_eq!(report.status().label, "1 input not imported");
    }

    #[test]
    fn open_and_add_data_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("study.tuff");
        create(&path, &all_inputs()[..1]).unwrap();

        let opened = open(&path).unwrap();
        assert_eq!(
            (opened.summary.psms, opened.status().kind),
            (5, StatusKind::Ready)
        );

        let added = add(&path, &all_inputs()[2..3]).unwrap();
        assert_eq!(added.summary.psms, 27);
        assert_eq!(added.status().label, "Saved");
    }

    #[test]
    fn replacing_is_limited_to_projects() {
        let dir = tempfile::tempdir().unwrap();
        let precious = dir.path().join("notes.tuff");
        std::fs::write(&precious, b"precious").unwrap();

        let error = create(&precious, &[]).unwrap_err();

        assert!(error.contains("is not a Tuffold project"), "{error}");
        assert_eq!(std::fs::read(&precious).unwrap(), b"precious");
    }

    #[test]
    fn an_existing_project_chosen_in_the_save_dialog_is_replaced() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("study.tuff");
        create(&path, &all_inputs()[..1]).unwrap();

        let fresh = create(&path, &all_inputs()[2..3]).unwrap();

        assert_eq!((fresh.imports.len(), fresh.summary.psms), (1, 22));
    }

    #[test]
    fn importing_validates_and_thresholds_apply() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("study.tuff");
        // The pyteomics mzIdentML has decoys but no score on every rank-1
        // PSM (see tuffold-engine's tests): validated, nothing ranked.
        let report = create(&path, &all_inputs()[2..3]).unwrap();
        let ValidationState::Done(summary) = &report.validation else {
            panic!("expected validation, got {:?}", report.validation);
        };
        assert_eq!(summary.engines[0].method, "unscored");

        let loose = Level::ProteinGroups.with(&Thresholds::default(), 0.05);
        let report = set_thresholds(&path, &loose).unwrap();
        let ValidationState::Done(summary) = &report.validation else {
            panic!("expected validation");
        };
        assert_eq!(summary.settings.thresholds, loose);
    }

    #[test]
    fn opening_never_validates_but_validate_does() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("study.tuff");
        {
            // Imported without the app (as the CLI does): never validated.
            let mut project = Project::create(&path).unwrap();
            import_inputs(&mut project, &all_inputs()[2..3], None);
        }
        let before = std::fs::read(&path).unwrap();
        let report = open(&path).unwrap();
        assert_eq!(report.validation, ValidationState::Pending { stale: false });
        assert_eq!(
            std::fs::read(&path).unwrap(),
            before,
            "open wrote to the file"
        );

        let report = validate(&path).unwrap();
        assert!(matches!(report.validation, ValidationState::Done(_)));

        // Data imported since (outside the app) makes the results stale.
        let mut project = Project::open(&path).unwrap();
        import_inputs(&mut project, &all_inputs()[1..2], None);
        drop(project);
        assert_eq!(
            open(&path).unwrap().validation,
            ValidationState::Pending { stale: true }
        );
        let empty = create(&dir.path().join("empty.tuff"), &[]).unwrap();
        assert_eq!(empty.validation, ValidationState::NoData);
    }

    fn summary_with(protein_groups: LevelCounts, min_peptides: u32) -> ValidationSummary {
        ValidationSummary {
            settings: tuffold_validate::Settings {
                rescore: true,
                thresholds: Thresholds {
                    min_peptides,
                    ..Thresholds::default()
                },
            },
            validated_at: Some(0),
            engines: Vec::new(),
            peptide_scoring: None,
            notes: Vec::new(),
            spectra: LevelCounts {
                targets: 1200,
                decoys: 11,
                lowest_q: Some(0.0008),
            },
            peptides: LevelCounts::default(),
            protein_groups,
            localization: Default::default(),
        }
    }

    #[test]
    fn level_views_explain_thresholds_the_estimate_cannot_reach() {
        let groups = LevelCounts {
            targets: 0,
            decoys: 0,
            lowest_q: Some(0.025),
        };
        let views = level_views(&summary_with(groups, 2));
        assert_eq!(views[0].caption, "at 1% FDR · 11 decoys passing");
        assert_eq!(views[0].hint, None);
        assert_eq!(
            views[2].caption,
            "at 1% FDR · 0 decoys passing · ≥2 peptides"
        );
        assert_eq!(
            views[2].hint.as_deref(),
            Some(
                "None can pass 1%: with this few protein groups, the lowest q-value reached is 2.5%."
            )
        );
        let changed = Level::ProteinGroups.with(&Thresholds::default(), 0.05);
        assert_eq!((changed.protein, changed.psm), (0.05, 0.01));
    }

    #[test]
    fn probability_thresholds_caption_their_estimated_fdr() {
        let groups = LevelCounts {
            targets: 400,
            decoys: 2,
            lowest_q: Some(0.0025),
        };
        let mut summary = summary_with(groups, 1);
        summary.settings.thresholds = Thresholds::probability_default();
        let views = level_views(&summary);
        // (2 + 1) decoys over 400 targets: 0.75%; (11 + 1) over 1,200: 1%.
        assert_eq!(
            views[2].caption,
            "at ≥99% probability · est. FDR 0.75% (2 decoys)"
        );
        assert_eq!(
            views[0].caption,
            "at ≥95% probability · est. FDR 1% (11 decoys)"
        );
        assert!(views.iter().all(|view| view.hint.is_none()));
        // Switching modes starts from that mode's defaults, keeping min peptides.
        let fdr = Thresholds {
            min_peptides: 2,
            ..Thresholds::default()
        };
        let switched = switched(fdr, ThresholdMode::Probability);
        assert_eq!(
            switched,
            Thresholds {
                min_peptides: 2,
                ..Thresholds::probability_default()
            }
        );
        assert_eq!(
            super::switched(switched, ThresholdMode::Probability),
            switched
        );
        assert_eq!(choices(ThresholdMode::Probability), PROBABILITY_CHOICES);
    }

    #[test]
    fn the_demo_opens_validated_with_protein_groups() {
        let dir = tempfile::tempdir().unwrap();
        let report = create_demo(&dir.path().join("demo")).unwrap();
        assert_eq!(report.path, dir.path().join("demo.tuff"));
        let ValidationState::Done(summary) = &report.validation else {
            panic!("the demo is validated: {:?}", report.validation);
        };
        assert!(summary.protein_groups.targets > 100);
        assert_eq!(report.imports[0].format, "synthetic");
    }

    #[test]
    fn extension_is_appended_once() {
        assert_eq!(
            project_path(Path::new("/d/study")),
            PathBuf::from("/d/study.tuff")
        );
        assert_eq!(
            project_path(Path::new("/d/study.TUFF")),
            PathBuf::from("/d/study.TUFF")
        );
        assert_eq!(
            project_path(Path::new("/d/v1.2")),
            PathBuf::from("/d/v1.2.tuff")
        );
    }

    #[test]
    fn the_workbook_exports_and_says_what_it_wrote() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        tuffold_engine::create_demo_project(&path).unwrap();
        let message = export_workbook(&path, &dir.path().join("demo.xlsx")).unwrap();
        assert!(message.starts_with("Exported "), "{message}");
        // 133 of the demo's 150 expressed proteins pass (135 before folds
        // went by peptide sequence; none of either list is unexpressed).
        assert!(message.contains("Proteins (133)"), "{message}");
        let error = export_workbook(&path, &dir.path().join("missing/x.xlsx")).unwrap_err();
        assert!(error.contains("missing"), "{error}");
    }
}
