//! View model for the Quant tab: loading quantities (spectral counts, and
//! MS1 intensity where the project has MS1 data), choosing a comparison, its
//! measure and test, the volcano plot's points, saving sample categories,
//! and exporting. Pure data and engine calls; no GPUI.

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use tuffold_core::reporter::{Impurities, Plex};
use tuffold_engine::{
    IntensityQuant, MbrError, MbrRunReport, ReporterQuant, match_between_runs, quantify_intensity,
    quantify_project, quantify_reporters, write_anova_tsv, write_intensity_tsv, write_qc_tsv,
    write_quant_tsv, write_reporter_tsv,
};
use tuffold_quant::anova::{AnovaComparison, GroupValues};
use tuffold_quant::heatmap::{Heatmap, column_order, heatmap};
use tuffold_quant::intensity::IntensityComparison;
use tuffold_quant::qc::SampleQc;
use tuffold_quant::{Comparison, ComparisonRow, PrimaryTest, Quantities};
use tuffold_store::{PlexDefinition, Project};

/// q-values at or below this are drawn and counted as significant.
pub const SIGNIFICANCE: f64 = 0.05;
/// Points this many on the volcano get a label (the most significant).
pub const LABELLED: usize = 8;
/// Labels closer than this (fraction of the plot, in both axes) are skipped.
const LABEL_GAP: f64 = 0.06;
/// The threshold line's label, as a fraction of the plot (width, height above
/// the line), and the space a dot and a point's label take (offsets from the
/// dot, same units). About 6 mono characters at 10 px on a ~700 × 320 px plot.
const LINE_LABEL: (f64, f64) = (0.08, 0.06);
const DOT: (f64, f64) = (0.008, 0.013);
const POINT_LABEL: [f64; 4] = [0.008, 0.1, 0.012, 0.057];

/// The heat map shows the groups differing at [`SIGNIFICANCE`], up to this
/// many (the most significant).
pub const HEAT_ROWS: usize = 60;
/// Without any significant group, the heat map shows this many of the most
/// significant ones, saying so.
pub const HEAT_FALLBACK: usize = 30;

/// Two categories against each other, or every category at once.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scope {
    Pair,
    All,
    /// Sample QC: PCA and correlations, no comparison.
    Qc,
}

/// What the Quant tab shows: scope, measure, and (for spectral counts in
/// two-category comparisons) the test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QuantChoice {
    pub scope: Scope,
    pub method: Method,
    pub test: tuffold_quant::PrimaryTest,
}

/// What a comparison measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    Counts,
    Intensity,
    Reporters,
}

pub fn plex_title(plex: Plex) -> &'static str {
    match plex {
        Plex::Tmt6 => "TMT6",
        Plex::Tmt6Original => "TMT6 original",
        Plex::Tmt10 => "TMT10",
        Plex::Tmt11 => "TMT11",
        Plex::TmtPro16 => "TMTpro 16",
        Plex::TmtPro18 => "TMTpro 18",
        Plex::Itraq4 => "iTRAQ 4",
        Plex::Itraq8 => "iTRAQ 8",
    }
}

/// All measures of a project. Intensity needs MS1 data and reporters need a
/// defined plex; without them the reason is kept, to show when chosen.
#[derive(Debug, Clone, PartialEq)]
pub struct QuantData {
    pub counts: Quantities,
    pub intensity: Result<IntensityQuant, String>,
    pub reporters: Result<ReporterQuant, String>,
    /// The project's first plex, if any.
    pub plex: Option<PlexDefinition>,
    /// Every plex of the project (reporters are bridged across them).
    pub plexes: Vec<PlexDefinition>,
    /// Every MS run in the project, in name order.
    pub runs: Vec<String>,
}

/// A run in the TMT plex editor: available for the edited plex, or already
/// taken by another plex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlexRun {
    pub name: String,
    /// The other plex this run already belongs to, if any (the run is then
    /// disabled for this plex).
    pub in_other_plex: Option<String>,
}

impl QuantData {
    /// Distinct categories, in first-seen order.
    pub fn categories(&self) -> Vec<String> {
        self.counts.categories()
    }

    /// Distinct categories for `method`: when `Method::Reporters` is active
    /// and a plex is defined, only the categories of the project's plex
    /// channels.
    pub fn categories_for(&self, method: Method) -> Vec<String> {
        if method == Method::Reporters {
            if let Ok(quant) = &self.reporters {
                return quant.quantities.quantities.categories();
            }
            if !self.plexes.is_empty() {
                let mut found = Vec::new();
                for plex in &self.plexes {
                    for (_, sample) in &plex.channels {
                        if let Some(category) = self
                            .counts
                            .samples
                            .iter()
                            .find(|s| &s.name == sample)
                            .and_then(|s| s.category.as_ref())
                            && !found.contains(category)
                        {
                            found.push(category.clone());
                        }
                    }
                }
                return found;
            }
        }
        self.counts.categories()
    }

    /// The comparison pair to pick first for `method`.
    pub fn default_pair_for(&self, method: Method) -> Option<(String, String)> {
        if method == Method::Reporters {
            let categories = self.categories_for(Method::Reporters);
            if let [a, b, ..] = categories.as_slice() {
                return Some((a.clone(), b.clone()));
            }
        }
        default_pair(&self.counts)
    }

    /// The stored plex named `name`, if any.
    pub fn plex_named(&self, name: &str) -> Option<&PlexDefinition> {
        let trimmed = name.trim();
        self.plexes.iter().find(|p| p.name == trimmed)
    }

    /// Every run of the project and whether another plex (anything other than
    /// `editing`, when `Some`) already holds it.
    pub fn runs_for_plex(&self, editing: Option<&str>) -> Vec<PlexRun> {
        plex_runs(&self.runs, &self.plexes, editing)
    }
}

/// Every run in `runs` and which other plex in `plexes` (anything other than
/// `editing`, when `Some`) already holds it.
pub fn plex_runs(
    runs: &[String],
    plexes: &[PlexDefinition],
    editing: Option<&str>,
) -> Vec<PlexRun> {
    let editing = editing.map(str::trim);
    runs.iter()
        .map(|run| {
            let in_other_plex = plexes
                .iter()
                .find(|p| Some(p.name.as_str()) != editing && p.runs.contains(run))
                .map(|p| p.name.clone());
            PlexRun {
                name: run.clone(),
                in_other_plex,
            }
        })
        .collect()
}

/// Whether `plex_name` is the only plex in `plexes` (a new plex when none
/// exist yet, or the single existing plex being edited): only then may an
/// empty run selection mean "every run in the project".
pub fn empty_runs_take_all(plexes: &[PlexDefinition], plex_name: &str) -> bool {
    let name = plex_name.trim();
    plexes.is_empty() || (plexes.len() == 1 && plexes[0].name == name)
}

/// Check the name typed for a new plex against the existing plexes.
pub fn check_new_plex_name(name: &str, plexes: &[PlexDefinition]) -> Result<String, String> {
    let trimmed = name.trim();
    if trimmed.is_empty() {
        return Err("Name the new plex.".into());
    }
    if plexes.iter().any(|p| p.name == trimmed) {
        return Err(format!("A plex named '{trimmed}' already exists."));
    }
    Ok(trimmed.to_string())
}

fn quantify(project: &Project) -> Result<QuantData, String> {
    let plexes = project.plexes().map_err(|e| e.to_string())?;
    let runs = project
        .summary()
        .map_err(|e| e.to_string())?
        .runs
        .into_iter()
        .map(|r| r.name)
        .collect();
    Ok(QuantData {
        counts: quantify_project(project).map_err(|e| e.to_string())?,
        intensity: quantify_intensity(project).map_err(|e| e.to_string()),
        reporters: quantify_reporters(project, None).map_err(|e| e.to_string()),
        plex: plexes.first().cloned(),
        plexes,
        runs,
    })
}

pub fn load(path: &Path) -> Result<QuantData, String> {
    quantify(&Project::open(path).map_err(|e| e.to_string())?)
}

/// What matching between runs did, and the quantities afterwards.
#[derive(Debug, Clone, PartialEq)]
pub struct MbrOutcome {
    /// The project quantified again. Present even when matching failed,
    /// because the engine saves each run as it goes: a failure on a later
    /// file leaves the earlier runs matched, and the tab must show that.
    pub data: Result<QuantData, String>,
    /// The per-run summary, or why matching stopped.
    pub result: Result<String, String>,
}

/// Match between the runs of the given spectra files (mzML with MS1 scans,
/// or timsTOF `.d` folders named like their runs), then quantify again.
pub fn match_runs(path: &Path, files: &[PathBuf]) -> MbrOutcome {
    match Project::open(path) {
        Ok(mut project) => {
            let result = match_between_runs(&mut project, files);
            conclude(&project, result, files.len())
        }
        Err(e) => {
            let message = format!("{MBR_STOPPED}{e}");
            MbrOutcome {
                data: Err(message.clone()),
                result: Err(message),
            }
        }
    }
}

const MBR_STOPPED: &str = "Match between runs stopped: ";

/// Quantify the project again after matching ran with `result` over
/// `files` files, and word what happened. Quantifying always happens: the
/// engine saves run by run, so even a failed matching may have changed what
/// the tab should show.
fn conclude(
    project: &Project,
    result: Result<Vec<MbrRunReport>, MbrError>,
    files: usize,
) -> MbrOutcome {
    let matched = result.is_ok();
    let data = quantify(project).map_err(|e| {
        if matched {
            format!("Matching finished, but quantifying again failed: {e}")
        } else {
            e
        }
    });
    let result = match result {
        Ok(reports) => Ok(mbr_summary(&reports)),
        Err(error) => Err(stopped_message(&error, files)),
    };
    MbrOutcome { data, result }
}

/// Why matching stopped. With several files, runs matched before the one that
/// failed were saved and stay matched, unless nothing had started (an
/// unvalidated project).
fn stopped_message(error: &MbrError, files: usize) -> String {
    let kept = if files > 1 && !matches!(error, MbrError::NotValidated) {
        " Runs matched before it were kept."
    } else {
        ""
    };
    let error = error.to_string();
    format!("{MBR_STOPPED}{}.{kept}", error.trim_end_matches('.'))
}

/// A short report of matching between runs: a headline over every run, then
/// one line per run with how many of the precursors looked for were
/// transferred, and what was off (a ranking by retention-time error alone,
/// damaged MS1 scans left out).
pub fn mbr_summary(reports: &[MbrRunReport]) -> String {
    let looked: usize = reports.iter().map(|r| r.summary.candidates).sum();
    let transferred: usize = reports.iter().map(|r| r.summary.accepted).sum();
    let mut lines = vec![format!(
        "Matched between runs in {} {}: {transferred} of {looked} {} looked for {} transferred \
         at {}% FDR.",
        reports.len(),
        if reports.len() == 1 { "run" } else { "runs" },
        if looked == 1 {
            "precursor"
        } else {
            "precursors"
        },
        if looked == 1 { "was" } else { "were" },
        tuffold_store::TRANSFER_FDR * 100.0
    )];
    for r in reports {
        let mut parts = vec![format!(
            "{}: {} of {} transferred",
            r.run, r.summary.accepted, r.summary.candidates
        )];
        if !r.summary.trained {
            parts.push("ranked by retention-time error (too few known peaks to train)".into());
        }
        if r.ms1_scans_damaged > 0 {
            parts.push(format!(
                "{} damaged MS1 {} left out",
                r.ms1_scans_damaged,
                if r.ms1_scans_damaged == 1 {
                    "scan"
                } else {
                    "scans"
                }
            ));
        }
        lines.push(parts.join(" · "));
    }
    lines.join("\n")
}

/// Save `(biosample, category)` assignments (an empty category clears it)
/// and quantify again.
pub fn save_categories(path: &Path, assignments: &[(String, String)]) -> Result<QuantData, String> {
    let mut project = Project::open(path).map_err(|e| e.to_string())?;
    let pairs: Vec<(&str, Option<&str>)> = assignments
        .iter()
        .map(|(name, category)| (name.as_str(), Some(category.as_str())))
        .collect();
    tuffold_engine::set_categories(&mut project, &pairs).map_err(|e| e.to_string())?;
    quantify(&project)
}

/// What saving a plex does with its reagent lot's impurity certificate.
#[derive(Debug, Clone, PartialEq)]
pub enum Certificate {
    /// The plex's own, if it has one for this reagent set.
    Keep,
    /// Read from this file.
    Load(PathBuf),
    Remove,
}

/// The reference-channel list the plex editor sends on save (`Some`), or
/// `None` when [`build_plex_save`] should keep the stored plex's references
/// via [`tuffold_engine::kept_references`].
///
/// Stored references are kept only when editing an existing plex without
/// changing its reagent set or touching its `ref` toggles; creating a new
/// plex, touching `ref`, or switching reagent sets sends the shown
/// references that sit on a non-empty channel.
pub fn references_for_save(
    edited: Option<&PlexDefinition>,
    reagents: Plex,
    channels: &[(String, String)],
    shown_references: &[String],
    references_touched: bool,
) -> Option<Vec<String>> {
    let keep_stored = !references_touched && edited.is_some_and(|p| p.plex == reagents);
    (!keep_stored).then(|| {
        shown_references
            .iter()
            .filter(|r| {
                channels
                    .iter()
                    .any(|(c, s)| c == *r && !s.trim().is_empty())
            })
            .cloned()
            .collect()
    })
}

/// Build the [`PlexDefinition`] a save sends to [`tuffold_engine::set_plex`],
/// checking whether a new plex's name clashes or whether runs must be picked
/// first because other plexes exist. When `references` is `None`, the stored
/// plex's references for `plex` are kept via [`tuffold_engine::kept_references`].
#[allow(clippy::too_many_arguments)]
pub fn build_plex_save(
    project: &Project,
    name: &str,
    is_new: bool,
    plex: Plex,
    runs: &[String],
    channels: &[(String, String)],
    references: Option<&[String]>,
    certificate: Certificate,
) -> Result<PlexDefinition, String> {
    let plexes = project.plexes().map_err(|e| e.to_string())?;
    let name = if is_new {
        check_new_plex_name(name, &plexes)?
    } else {
        let trimmed = name.trim();
        if trimmed.is_empty() {
            "plex".to_string()
        } else {
            trimmed.to_string()
        }
    };
    if runs.is_empty() && !empty_runs_take_all(&plexes, &name) {
        return Err(format!(
            "Pick at least one run for plex '{name}' (other plexes exist, so an empty run list \
             cannot take every run)."
        ));
    }
    let impurities = match certificate {
        Certificate::Keep => {
            tuffold_engine::kept_impurities(project, &name, plex).map_err(|e| e.to_string())?
        }
        Certificate::Load(file) => {
            let text =
                std::fs::read_to_string(&file).map_err(|e| format!("{}: {e}", file.display()))?;
            Some(Impurities::parse(&text, plex).map_err(|e| format!("{}: {e}", file.display()))?)
        }
        Certificate::Remove => None,
    };
    let references = match references {
        Some(given) => given.to_vec(),
        None => tuffold_engine::kept_references(project, &name, plex).map_err(|e| e.to_string())?,
    };
    Ok(PlexDefinition {
        name,
        plex,
        runs: runs.to_vec(),
        channels: channels.to_vec(),
        impurities,
        references,
    })
}

/// Save a TMT plex definition (reagent set, runs, channel -> biosample name,
/// reference channels bridging it to other plexes, and the lot's impurity
/// certificate) and quantify again.
#[allow(clippy::too_many_arguments)]
pub fn save_plex_named(
    path: &Path,
    name: &str,
    is_new: bool,
    plex: Plex,
    runs: &[String],
    channels: &[(String, String)],
    references: Option<&[String]>,
    certificate: Certificate,
) -> Result<QuantData, String> {
    let mut project = Project::open(path).map_err(|e| e.to_string())?;
    let definition = build_plex_save(
        &project,
        name,
        is_new,
        plex,
        runs,
        channels,
        references,
        certificate,
    )?;
    tuffold_engine::set_plex(&mut project, &definition).map_err(|e| e.to_string())?;
    quantify(&project)
}

/// Save the project's first TMT plex (or create `"plex"` when none exists)
/// and quantify again.
pub fn save_plex(
    path: &Path,
    plex: Plex,
    channels: &[(String, String)],
    references: &[String],
    certificate: Certificate,
) -> Result<QuantData, String> {
    let project = Project::open(path).map_err(|e| e.to_string())?;
    let existing = project
        .plexes()
        .map_err(|e| e.to_string())?
        .into_iter()
        .next();
    drop(project);
    let (name, runs) = match existing {
        Some(existing) => (existing.name, existing.runs),
        None => ("plex".to_string(), Vec::new()),
    };
    save_plex_named(
        path,
        &name,
        false,
        plex,
        &runs,
        channels,
        Some(references),
        certificate,
    )
}

/// Remove a named TMT plex and quantify again.
pub fn remove_plex(path: &Path, name: &str) -> Result<QuantData, String> {
    let mut project = Project::open(path).map_err(|e| e.to_string())?;
    tuffold_engine::remove_plex(&mut project, name).map_err(|e| e.to_string())?;
    quantify(&project)
}

/// What is exported: the table of the measure shown, with its comparison.
pub enum Export {
    Counts(Quantities, Option<Comparison>),
    Intensity(IntensityQuant, Option<IntensityComparison>),
    Reporters(ReporterQuant, Option<IntensityComparison>),
    Anova(AnovaComparison),
    Qc(SampleQc),
}

pub fn export(path: &Path, what: &Export) -> Result<(), String> {
    let file = File::create(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out = BufWriter::new(file);
    match what {
        Export::Counts(quantities, comparison) => {
            write_quant_tsv(&mut out, quantities, comparison.as_ref())
        }
        Export::Intensity(quant, comparison) => {
            write_intensity_tsv(&mut out, &quant.quantities, comparison.as_ref())
        }
        Export::Reporters(quant, comparison) => {
            write_reporter_tsv(&mut out, quant, comparison.as_ref())
        }
        Export::Anova(anova) => write_anova_tsv(&mut out, anova),
        Export::Qc(qc) => write_qc_tsv(&mut out, qc),
    }
    .map_err(|e| format!("{}: {e}", path.display()))
}

/// The comparison to show first: the first two categories, as listed, that
/// have any counted spectra (or any two categories when only channel
/// biosamples are categorized).
pub fn default_pair(quantities: &Quantities) -> Option<(String, String)> {
    let counted: Vec<String> = quantities
        .categories()
        .into_iter()
        .filter(|category| {
            quantities
                .samples
                .iter()
                .zip(&quantities.totals)
                .any(|(sample, &total)| sample.category.as_ref() == Some(category) && total > 0)
        })
        .collect();
    if let [a, b, ..] = counted.as_slice() {
        return Some((a.clone(), b.clone()));
    }
    match quantities.categories().as_slice() {
        [a, b, ..] => Some((a.clone(), b.clone())),
        _ => None,
    }
}

/// The q-value of `test`; `None` when that test could not run for the row
/// (Welch without replicates or variance) — unknown, never shown as 1.
pub fn q_of(row: &ComparisonRow, test: PrimaryTest) -> Option<f64> {
    match test {
        PrimaryTest::Welch => row.welch_q,
        PrimaryTest::Fisher => Some(row.fisher_q),
        PrimaryTest::MannWhitney => row.mann_whitney_q,
    }
}

/// One protein group of a comparison, whichever the measure.
#[derive(Debug, Clone, PartialEq)]
pub struct CompareRow {
    pub group_id: i64,
    pub accession: String,
    /// Mean normalized count, or mean log2 intensity.
    pub mean_a: Option<f64>,
    pub mean_b: Option<f64>,
    pub log2_fold: Option<f64>,
    /// `None` where the test couldn't run.
    pub q: Option<f64>,
}

/// A comparison as the Quant tab shows it.
#[derive(Debug, Clone, PartialEq)]
pub struct Shown {
    pub a: String,
    pub b: String,
    pub samples_a: Vec<String>,
    pub samples_b: Vec<String>,
    pub rows: Vec<CompareRow>,
    /// The test behind the q-values, in words.
    pub test: String,
    /// What the means are, in words.
    pub means: &'static str,
}

/// A spectral-count comparison under `test`, most significant first.
pub fn shown_counts(comparison: &Comparison, test: PrimaryTest) -> Shown {
    let mut rows: Vec<CompareRow> = comparison
        .rows
        .iter()
        .map(|row| CompareRow {
            group_id: row.group_id,
            accession: row.accession.clone(),
            mean_a: Some(row.mean_a),
            mean_b: Some(row.mean_b),
            log2_fold: Some(row.log2_fold),
            q: q_of(row, test),
        })
        .collect();
    rows.sort_by(|x, y| {
        x.q.unwrap_or(f64::INFINITY)
            .total_cmp(&y.q.unwrap_or(f64::INFINITY))
    });
    Shown {
        a: comparison.a.clone(),
        b: comparison.b.clone(),
        samples_a: comparison.samples_a.clone(),
        samples_b: comparison.samples_b.clone(),
        rows,
        test: test.describe().to_string(),
        means: "Mean normalized counts",
    }
}

/// An MS1 intensity comparison (already most significant first).
pub fn shown_intensity(comparison: &IntensityComparison) -> Shown {
    Shown {
        a: comparison.a.clone(),
        b: comparison.b.clone(),
        samples_a: comparison.samples_a.clone(),
        samples_b: comparison.samples_b.clone(),
        rows: comparison
            .rows
            .iter()
            .map(|row| CompareRow {
                group_id: row.group_id,
                accession: row.accession.clone(),
                mean_a: row.mean_a,
                mean_b: row.mean_b,
                log2_fold: row.log2_fold,
                q: row.q,
            })
            .collect(),
        test: "Welch's t-test on log2 median-polish intensities".into(),
        means: "Mean log2 intensity (median polish)",
    }
}

/// A TMT reporter-ion comparison (already most significant first).
pub fn shown_reporters(comparison: &IntensityComparison) -> Shown {
    Shown {
        test: "Welch's t-test on log2 reporter intensities".into(),
        means: "Mean log2 reporter intensity",
        ..shown_intensity(comparison)
    }
}

/// Every category compared at once, with its heat map.
#[derive(Debug, Clone, PartialEq)]
pub struct AllShown {
    pub anova: AnovaComparison,
    pub heatmap: Heatmap,
    /// What the heat map's rows are, in words.
    pub caption: String,
    /// The test, in words.
    pub test: String,
}

/// One-way ANOVA across every category by `method`, and a heat map of the
/// groups that differ (or, without any, of the most significant, labelled
/// as such).
pub fn shown_all(data: &QuantData, method: Method) -> Result<AllShown, String> {
    let (samples, (with_data, groups), anova, values) = match method {
        Method::Counts => (
            &data.counts.samples,
            data.counts.log2_values(),
            data.counts.anova(),
            "log2(normalized counts + 1)",
        ),
        Method::Intensity => {
            let quantities = &data.intensity.as_ref().map_err(Clone::clone)?.quantities;
            (
                &quantities.samples,
                quantities.log2_values(),
                quantities.anova(),
                "log2 median-polish intensities",
            )
        }
        Method::Reporters => {
            let quantities = &data
                .reporters
                .as_ref()
                .map_err(Clone::clone)?
                .quantities
                .quantities;
            (
                &quantities.samples,
                quantities.log2_values(),
                quantities.anova(),
                "log2 reporter intensities",
            )
        }
    };
    let anova = anova.map_err(|e| e.to_string())?;
    let significant: Vec<i64> = anova
        .rows
        .iter()
        .filter(|r| r.q.is_some_and(|q| q <= SIGNIFICANCE))
        .map(|r| r.group_id)
        .collect();
    let (ids, caption): (Vec<i64>, String) = if significant.is_empty() {
        (
            anova
                .rows
                .iter()
                .filter(|r| r.p.is_some())
                .take(HEAT_FALLBACK)
                .map(|r| r.group_id)
                .collect(),
            format!(
                "None differ at {}% FDR: the {HEAT_FALLBACK} most significant, for their \
                 pattern only",
                SIGNIFICANCE * 100.0
            ),
        )
    } else if significant.len() > HEAT_ROWS {
        (
            significant[..HEAT_ROWS].to_vec(),
            format!(
                "The {HEAT_ROWS} most significant of {} groups differing at {}% FDR",
                significant.len(),
                SIGNIFICANCE * 100.0
            ),
        )
    } else {
        (
            significant.clone(),
            format!(
                "The {} groups differing at {}% FDR",
                significant.len(),
                SIGNIFICANCE * 100.0
            ),
        )
    };
    let chosen: Vec<GroupValues> = ids
        .iter()
        .filter_map(|id| groups.iter().find(|g| g.group_id == *id).cloned())
        .collect();
    Ok(AllShown {
        heatmap: heatmap(samples, &column_order(samples, &with_data), &chosen),
        anova,
        caption,
        test: format!("one-way ANOVA on {values}"),
    })
}

/// Sample QC of one measure.
#[derive(Debug, Clone, PartialEq)]
pub struct QcShown {
    pub qc: SampleQc,
    /// The values it's on, in words.
    pub values: &'static str,
}

pub fn shown_qc(data: &QuantData, method: Method) -> Result<QcShown, String> {
    Ok(match method {
        Method::Counts => QcShown {
            qc: data.counts.sample_qc(),
            values: "log2(normalized counts + 1)",
        },
        Method::Intensity => QcShown {
            qc: data
                .intensity
                .as_ref()
                .map_err(Clone::clone)?
                .quantities
                .sample_qc(),
            values: "log2 median-polish intensities",
        },
        Method::Reporters => QcShown {
            qc: data
                .reporters
                .as_ref()
                .map_err(Clone::clone)?
                .quantities
                .quantities
                .sample_qc(),
            values: "log2 reporter intensities",
        },
    })
}

/// A line on the MS1 data behind intensities: scans, precursors measured,
/// and what normalization rested on.
pub fn ms1_status(quant: &IntensityQuant) -> String {
    let with_ms1 = quant
        .runs
        .iter()
        .filter(|r| r.ms1_scans.is_some_and(|n| n > 0))
        .count();
    let traced: usize = quant.runs.iter().map(|r| r.traced).sum();
    let measured: usize = quant.runs.iter().map(|r| r.measured).sum();
    let untraced: usize = quant.runs.iter().map(|r| r.untraced).sum();
    let mut parts = vec![format!(
        "MS1 data in {with_ms1} of {} runs · {measured} of {traced} traced precursors have an \
         elution peak",
        quant.runs.len()
    )];
    if untraced > 0 {
        parts.push(format!(
            "{untraced} identified since their spectra files were added (add them again to \
             trace those)"
        ));
    }
    if !quant.incomplete.is_empty() {
        parts.push(format!(
            "left out (some of their runs have no MS1 data): {}",
            quant.incomplete.join(", ")
        ));
    }
    if quant.transferred > 0 {
        parts.push(format!(
            "{} values matched between runs (1% FDR)",
            quant.transferred
        ));
    }
    let q = &quant.quantities;
    let unnormalized: Vec<&str> = q
        .samples
        .iter()
        .zip(&q.normalized)
        .filter(|(s, done)| !**done && !quant.incomplete.contains(&s.name))
        .map(|(s, _)| s.name.as_str())
        .collect();
    parts.push(if unnormalized.len() == q.samples.len() {
        "not normalized (no two biosamples share enough peptides)".into()
    } else if unnormalized.is_empty() {
        format!(
            "normalized on {} peptides measured in two or more biosamples",
            q.shared_peptides
        )
    } else {
        format!(
            "normalized, except {} (too few peptides shared with the others)",
            unnormalized.join(", ")
        )
    });
    parts.join(" · ")
}

/// A line on the TMT plexes behind reporter intensities: reagent sets, runs,
/// PSMs with every channel, bridging, impurity correction, and channel median
/// normalization.
pub fn reporter_status(quant: &ReporterQuant) -> String {
    let r = &quant.quantities;
    let q = &r.quantities;
    let runs: usize = quant.plexes.iter().map(|p| p.runs.len()).sum();
    let runs = format!("{runs} {}", if runs == 1 { "run" } else { "runs" });
    let mut parts = vec![match quant.plexes.as_slice() {
        [plex] => format!(
            "{} ({}, {runs}) · {} protein groups from {} PSMs with all {} channels",
            plex_title(plex.plex),
            plex.name,
            q.rows.len(),
            r.psms_used,
            r.channels.len(),
        ),
        many => format!(
            "{} plexes ({}; {runs}) · {} protein groups from {} PSMs with every channel",
            many.len(),
            many.iter()
                .map(|p| format!("{} {}", p.name, plex_title(p.plex)))
                .collect::<Vec<_>>()
                .join(", "),
            q.rows.len(),
            r.psms_used,
        ),
    }];
    if r.psms_incomplete > 0 {
        parts.push(format!("{} missing a channel left out", r.psms_incomplete));
    }
    if quant.without_peaks > 0 {
        parts.push(format!(
            "{} passing spectra without peaks (attach the runs' mzML/MGF files)",
            quant.without_peaks
        ));
    }
    if quant.from_ms3 > 0 {
        parts.push(format!(
            "reporters of {} PSMs from their SPS-MS3 scans",
            quant.from_ms3
        ));
    }
    if let (Some(min), true) = (quant.min_purity, quant.impure > 0) {
        parts.push(format!(
            "{} spectra with precursor purity under {:.0}% left out",
            quant.impure,
            min * 100.0
        ));
    }
    if quant.plexes.len() > 1 {
        let mut bridged = "bridged through each plex's reference channels (left out)".to_string();
        if quant.unbridged > 0 {
            bridged.push_str(&format!(
                ", {} values without a reference value left out",
                quant.unbridged
            ));
        }
        parts.push(bridged);
    }
    let corrected: Vec<&str> = quant
        .plexes
        .iter()
        .filter(|p| p.impurities.is_some())
        .map(|p| p.name.as_str())
        .collect();
    parts.push(if corrected.is_empty() {
        "isotope impurities not corrected (no lot certificate)".into()
    } else if corrected.len() == quant.plexes.len() {
        "corrected for the lot's isotope impurities".into()
    } else {
        format!(
            "isotope impurities corrected in {} only",
            corrected.join(", ")
        )
    });
    parts.push("median-normalized across channels".into());
    parts.join(" · ")
}

/// When another reagent set with the same channel names would give more
/// PSMs every channel than the chosen one, say so.
pub fn reporter_reagent_hint(quant: &ReporterQuant) -> Option<String> {
    let (name, other, complete) = quant
        .better_reagents
        .iter()
        .max_by_key(|(_, _, complete)| *complete)?;
    let plex = quant.plexes.iter().find(|p| &p.name == name)?;
    let which = if quant.plexes.len() > 1 {
        format!("In plex {name}, at")
    } else {
        "At".into()
    };
    Some(format!(
        "{which} the {} ({other}) reporter m/z values, {complete} PSMs have every channel \
         (vs {} with {}). If it was labelled with {}, switch the reagent set and save the plex.",
        plex_title(*other),
        quant.quantities.psms_used,
        plex_title(plex.plex),
        plex_title(*other),
    ))
}

#[derive(Debug, Clone, PartialEq)]
pub struct VolcanoPoint {
    pub group_id: i64,
    pub accession: String,
    pub log2_fold: f64,
    /// -log10(q).
    pub height: f64,
    pub significant: bool,
    pub labelled: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Volcano {
    pub points: Vec<VolcanoPoint>,
    /// Symmetric x extent: folds run from -x_max to x_max.
    pub x_max: f64,
    pub y_max: f64,
    /// -log10 of [`SIGNIFICANCE`], where the threshold line goes.
    pub threshold: f64,
    /// Where the threshold line's label starts, as a fraction of the plot's
    /// width: the stretch of the line with the fewest dots and labels on it.
    pub threshold_label_x: f64,
    pub significant: usize,
    /// Rows the test could not score (no point drawn).
    pub untested: usize,
}

/// Points of a comparison: x is the log2 fold change, y is -log10(q). Rows
/// the test could not score are left out (counted in `untested`).
pub fn volcano(shown: &Shown) -> Volcano {
    let mut points: Vec<VolcanoPoint> = shown
        .rows
        .iter()
        .filter_map(|row| {
            let (q, log2_fold) = (row.q?, row.log2_fold?);
            Some(VolcanoPoint {
                group_id: row.group_id,
                accession: row.accession.clone(),
                log2_fold,
                height: -q.max(1e-300).log10(),
                significant: q <= SIGNIFICANCE,
                labelled: false,
            })
        })
        .collect();
    let untested = shown.rows.len() - points.len();
    let threshold = -SIGNIFICANCE.log10();
    let x_max = points
        .iter()
        .map(|p| p.log2_fold.abs())
        .fold(1.0_f64, f64::max)
        * 1.1;
    let y_max = points.iter().map(|p| p.height).fold(threshold, f64::max) * 1.15;
    // Label the most significant points, skipping any that would sit on a
    // label already placed (within LABEL_GAP of the plot in both axes).
    let mut order: Vec<usize> = (0..points.len())
        .filter(|&i| points[i].significant)
        .collect();
    order.sort_by(|&a, &b| points[b].height.total_cmp(&points[a].height));
    let mut placed: Vec<(f64, f64)> = Vec::new();
    for index in order {
        if placed.len() == LABELLED {
            break;
        }
        let at = (
            points[index].log2_fold / (2.0 * x_max),
            points[index].height / y_max,
        );
        let clear = placed
            .iter()
            .all(|(x, y)| (x - at.0).abs() > LABEL_GAP || (y - at.1).abs() > LABEL_GAP);
        if clear {
            points[index].labelled = true;
            placed.push(at);
        }
    }
    let significant = points.iter().filter(|p| p.significant).count();
    let threshold_label_x = line_label_x(&points, x_max, y_max, threshold);
    Volcano {
        points,
        x_max,
        y_max,
        threshold,
        threshold_label_x,
        significant,
        untested,
    }
}

/// The left edge (fraction of the width) of the stretch just above the
/// threshold line that overlaps the fewest dots and point labels. Where the
/// points crowd depends on the measure (reporter ions pile up at fold 0 near
/// the line, MS1 folds near ±2), so no fixed place stays clear. Ties go to
/// the right end, where a reference line is usually named.
fn line_label_x(points: &[VolcanoPoint], x_max: f64, y_max: f64, threshold: f64) -> f64 {
    let line = threshold / y_max;
    let overlaps = |left: f64, rect: [f64; 4]| {
        rect[0] < left + LINE_LABEL.0
            && rect[1] > left
            && rect[2] < line + LINE_LABEL.1
            && rect[3] > line
    };
    let rects: Vec<[f64; 4]> = points
        .iter()
        .flat_map(|p| {
            let (x, y) = ((p.log2_fold + x_max) / (2.0 * x_max), p.height / y_max);
            let dot = [x - DOT.0, x + DOT.0, y - DOT.1, y + DOT.1];
            let label = p.labelled.then(|| {
                [
                    x + POINT_LABEL[0],
                    x + POINT_LABEL[1],
                    y + POINT_LABEL[2],
                    y + POINT_LABEL[3],
                ]
            });
            std::iter::once(dot).chain(label)
        })
        .collect();
    let candidates = (0..=22).map(|i| 0.9 - f64::from(i) * 0.04);
    candidates
        .min_by_key(|&left| rects.iter().filter(|&&r| overlaps(left, r)).count())
        .unwrap_or(0.9)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tuffold_quant::{ProteinCounts, Sample, quantify};

    #[test]
    fn several_plexes_are_described_as_bridged() {
        let plex = |name: &str, impurities: bool| PlexDefinition {
            name: name.into(),
            plex: Plex::Tmt10,
            runs: vec![format!("{name}-1"), format!("{name}-2")],
            channels: vec![("126".into(), format!("{name}-pool"))],
            impurities: impurities.then(|| {
                let text: String = Plex::Tmt10
                    .channels()
                    .iter()
                    .map(|c| format!("{} 0 0 0 0\n", c.name))
                    .collect();
                Impurities::parse(&text, Plex::Tmt10).unwrap()
            }),
            references: vec!["126".into()],
        };
        let quant = ReporterQuant {
            plexes: vec![plex("a", true), plex("b", false)],
            quantities: tuffold_quant::reporter::ReporterQuantities::default(),
            better_reagents: vec![("b".into(), Plex::Tmt11, 40)],
            without_peaks: 0,
            from_ms3: 12,
            min_purity: Some(0.5),
            impure: 7,
            purity_unknown: 0,
            unbridged: 3,
        };
        let status = reporter_status(&quant);
        assert!(
            status.starts_with("2 plexes (a TMT10, b TMT10; 4 runs)"),
            "{status}"
        );
        assert!(
            status.contains("bridged through each plex's reference channels (left out), 3 values")
        );
        assert!(status.contains("isotope impurities corrected in a only"));
        assert!(status.contains("reporters of 12 PSMs from their SPS-MS3 scans"));
        assert!(status.contains("7 spectra with precursor purity under 50% left out"));
        let hint = reporter_reagent_hint(&quant).unwrap();
        assert!(hint.starts_with("In plex b, at the TMT11"), "{hint}");
    }

    fn quantities() -> Quantities {
        let sample = |name: &str, category: &str| Sample {
            name: name.into(),
            category: Some(category.into()),
        };
        let protein = |id, counts: &[usize]| ProteinCounts {
            group_id: id,
            accession: format!("P{id}"),
            length: None,
            counts: counts.to_vec(),
            ..ProteinCounts::default()
        };
        quantify(
            vec![sample("c1", "control"), sample("t1", "treated")],
            vec![
                protein(1, &[200, 900]),
                protein(2, &[500, 500]),
                protein(3, &[300, 30]),
            ],
        )
        .unwrap()
    }

    #[test]
    fn the_first_two_categories_are_compared_first() {
        assert_eq!(
            default_pair(&quantities()),
            Some(("control".into(), "treated".into()))
        );
    }

    #[test]
    fn the_volcano_marks_and_labels_significant_points() {
        let q = quantities();
        let comparison = q.compare("control", "treated").unwrap();
        let volcano = volcano(&shown_counts(&comparison, PrimaryTest::Fisher));
        let point = |id| volcano.points.iter().find(|p| p.group_id == id).unwrap();
        // P1 rises and P3 falls sharply; both are highly significant.
        assert!(point(1).log2_fold > 0.0 && point(1).significant && point(1).labelled);
        assert!(point(3).log2_fold < 0.0 && point(3).significant);
        assert!(volcano.x_max >= point(3).log2_fold.abs());
        assert!(volcano.y_max > point(1).height);
        assert!((volcano.threshold - 1.301_029_995_663_981).abs() < 1e-12);
        assert!(point(3).labelled, "far from P1, so labelled too");
        // Without replicates Welch scores nothing: no points, all untested.
        let welch = super::volcano(&shown_counts(&comparison, PrimaryTest::Welch));
        assert!(welch.points.is_empty());
        assert_eq!(welch.untested, 3);
    }

    #[test]
    fn intensity_comparisons_keep_untested_rows_off_the_plot() {
        use tuffold_quant::intensity::IntensityComparisonRow;
        let row = |id, fold: Option<f64>, q: Option<f64>| IntensityComparisonRow {
            group_id: id,
            accession: format!("P{id}"),
            mean_a: Some(20.0),
            mean_b: fold.map(|f| 20.0 + f),
            measured_a: 2,
            measured_b: usize::from(fold.is_some()) * 2,
            log2_fold: fold,
            p: q,
            q,
        };
        let comparison = IntensityComparison {
            a: "control".into(),
            b: "treated".into(),
            samples_a: vec!["c1".into(), "c2".into()],
            samples_b: vec!["t1".into(), "t2".into()],
            rows: vec![row(1, Some(2.0), Some(0.01)), row(2, None, None)],
        };
        let shown = shown_intensity(&comparison);
        assert_eq!(shown.rows[1].mean_b, None);
        let plot = volcano(&shown);
        assert_eq!(
            (plot.points.len(), plot.untested, plot.significant),
            (1, 1, 1)
        );
    }

    #[test]
    fn the_threshold_label_avoids_points_on_the_line() {
        let point = |log2_fold: f64, height: f64| VolcanoPoint {
            group_id: 0,
            accession: String::new(),
            log2_fold,
            height,
            significant: height >= -SIGNIFICANCE.log10(),
            labelled: false,
        };
        let threshold = -SIGNIFICANCE.log10();
        // Nothing near the line: the right end.
        assert_eq!(line_label_x(&[point(0.0, 10.0)], 2.2, 12.0, threshold), 0.9);
        // A cluster on the line at the right end pushes it left of the cluster.
        let crowded: Vec<VolcanoPoint> = (0..20)
            .map(|i| point(1.8 + f64::from(i) * 0.02, threshold + 0.1))
            .collect();
        let x = line_label_x(&crowded, 2.2, 12.0, threshold);
        assert!(x + LINE_LABEL.0 < (1.8 + 2.2) / 4.4 - DOT.0, "{x}");
        // Reporter ions: a null cluster around fold 0 straddling the line.
        let null: Vec<VolcanoPoint> = (0..40)
            .map(|i| {
                point(
                    -0.1 + f64::from(i) * 0.005,
                    threshold * f64::from(i % 4) / 2.0,
                )
            })
            .collect();
        let x = line_label_x(&null, 2.2, 12.0, threshold);
        assert!(x > 0.55 || x + LINE_LABEL.0 < 0.45, "{x}");
    }

    #[test]
    fn every_category_compares_with_a_heat_map_of_what_differs() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        tuffold_engine::create_demo_project(&path).unwrap();
        let data = load(&path).unwrap();
        let all = shown_all(&data, Method::Intensity).unwrap();
        assert_eq!(all.anova.categories, ["control", "treated"]);
        let significant = all
            .anova
            .rows
            .iter()
            .filter(|r| r.q.is_some_and(|q| q <= SIGNIFICANCE))
            .count();
        assert_eq!(all.heatmap.rows.len(), significant.min(HEAT_ROWS));
        assert!(all.caption.starts_with("The "), "{}", all.caption);
        assert_eq!(
            all.heatmap.columns,
            ["control-A", "control-B", "treated-A", "treated-B"]
        );
        // Spectral counts: the same comparison by counts.
        let counts = shown_all(&data, Method::Counts).unwrap();
        assert!(counts.test.contains("counts"));
        let qc = shown_qc(&data, Method::Counts).unwrap();
        assert_eq!(qc.qc.samples.len(), 4);
        assert!(qc.qc.pca.is_some());
        let qc_table = dir.path().join("qc.tsv");
        export(&qc_table, &Export::Qc(qc.qc)).unwrap();
        let text = std::fs::read_to_string(&qc_table).unwrap();
        assert_eq!(text.lines().count(), 5);
        assert!(text.starts_with("sample\tcategory\tmeasured_groups"));
        assert_eq!(text.lines().nth(1).unwrap().split('\t').count(), 10);
        let table = dir.path().join("anova.tsv");
        export(&table, &Export::Anova(all.anova.clone())).unwrap();
        let text = std::fs::read_to_string(&table).unwrap();
        assert_eq!(text.lines().count(), all.anova.rows.len() + 1);
    }

    #[test]
    fn groups_save_and_the_table_exports() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        tuffold_engine::create_demo_project(&path).unwrap();
        let assignments = vec![
            ("control-A".to_string(), "control".to_string()),
            ("control-B".to_string(), "control".to_string()),
            ("treated-A".to_string(), "treated".to_string()),
            ("treated-B".to_string(), "knockout".to_string()),
        ];
        let data = save_categories(&path, &assignments).unwrap();
        assert_eq!(data.categories(), ["control", "treated", "knockout"]);
        // The demo carries MS1 scans, so intensities come along.
        let intensity = data.intensity.as_ref().unwrap();
        assert!(ms1_status(intensity).starts_with("MS1 data in 4 of 4 runs"));
        // An empty category clears it: treated-B leaves every comparison.
        let mut cleared = assignments.clone();
        cleared[3].1 = String::new();
        let data = save_categories(&path, &cleared).unwrap();
        assert_eq!(data.categories(), ["control", "treated"]);

        let quantities = data.counts.clone();
        let comparison = quantities.compare("control", "treated").unwrap();
        assert_eq!(comparison.samples_b, ["treated-A"]);
        let table = dir.path().join("quant.tsv");
        export(
            &table,
            &Export::Counts(quantities.clone(), Some(comparison)),
        )
        .unwrap();
        let text = std::fs::read_to_string(&table).unwrap();
        assert_eq!(text.lines().count(), quantities.rows.len() + 1);
        let intensity = data.intensity.unwrap();
        let rows = intensity.quantities.rows.len();
        export(&table, &Export::Intensity(intensity, None)).unwrap();
        let text = std::fs::read_to_string(&table).unwrap();
        assert!(text.starts_with("accession\tintensity:"));
        assert_eq!(text.lines().count(), rows + 1);
        let error = export(
            &dir.path().join("missing/x.tsv"),
            &Export::Counts(quantities, None),
        )
        .unwrap_err();
        assert!(error.contains("missing"), "{error}");
    }

    #[test]
    fn plex_saves_hints_better_reagents_and_exports() {
        use tuffold_core::{
            BioSample, Evidence, ImportBatch, MsRun, PeakSpectrum, PeptideForm, Protein, Psm,
            Score, Software, SourceFormat, Spectrum, SpectrumKey,
        };
        use tuffold_store::ImportSource;
        use tuffold_validate::{Settings, ThresholdMode, Thresholds};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("tmt.tuff");
        let mut project = Project::create(&path).unwrap();
        let run = "pool_1";
        let psm = |scan: u32, seq: &str, acc: &str, score: f64| Psm {
            ms_run: run.into(),
            spectrum: SpectrumKey::Scan(scan),
            rank: 1,
            peptide: PeptideForm::unmodified(seq),
            charge: 2,
            calculated_mz: Some(500.0),
            experimental_mz: Some(500.0),
            pass_threshold: None,
            scores: vec![Score {
                name: "hyperscore".into(),
                accession: None,
                value: score,
            }],
            evidence: vec![Evidence {
                accession: acc.into(),
                start: None,
                end: None,
                pre: None,
                post: None,
            }],
        };
        let batch = ImportBatch {
            software: Software {
                name: "Sage".into(),
                version: None,
            },
            biosamples: vec![BioSample {
                name: run.into(),
                category: None,
            }],
            ms_runs: vec![MsRun {
                name: run.into(),
                biosample: run.into(),
                spectra_file: None,
            }],
            spectra: (1..=3)
                .map(|s| Spectrum {
                    ms_run: run.into(),
                    key: SpectrumKey::Scan(s),
                    source_id: None,
                    precursor_mz: Some(500.0),
                    precursor_charge: Some(2),
                    rt_seconds: Some(f64::from(s) * 10.0),
                })
                .collect(),
            psms: vec![
                psm(1, "PEPTIDEK", "P1", 50.0),
                psm(2, "ANOTHERK", "P2", 40.0),
                psm(3, "DECOYPEPK", "rev_P1", 10.0),
            ],
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: None,
                    sequence: None,
                    is_decoy: false,
                    decoy_explicit: false,
                },
                Protein {
                    accession: "P2".into(),
                    description: None,
                    sequence: None,
                    is_decoy: false,
                    decoy_explicit: false,
                },
                Protein {
                    accession: "rev_P1".into(),
                    description: None,
                    sequence: None,
                    is_decoy: true,
                    decoy_explicit: false,
                },
            ],
        };
        project
            .import(
                &ImportSource {
                    path: "/data/tmt.pep.xml".into(),
                    format: SourceFormat::PepXml,
                },
                &batch,
            )
            .unwrap();
        let mz: Vec<f64> = Plex::Tmt6Original
            .channels()
            .into_iter()
            .take(4)
            .map(|c| c.mz)
            .collect();
        let scan = |s: u32, intensities: [f32; 4]| PeakSpectrum {
            key: SpectrumKey::Scan(s),
            native_id: format!("scan={s}"),
            ms_level: 2,
            rt_seconds: None,
            precursor_mz: Some(500.0),
            precursor_charge: Some(2),
            mz: mz.clone(),
            intensity: intensities.to_vec(),
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        };
        let mut attacher = project.attach_spectra(run).unwrap();
        attacher
            .add(&scan(1, [1000.0, 2000.0, 4000.0, 16000.0]))
            .unwrap();
        attacher
            .add(&scan(2, [2000.0, 1000.0, 2000.0, 2000.0]))
            .unwrap();
        attacher.finish().unwrap();
        tuffold_engine::validate_project(
            &mut project,
            &Settings {
                rescore: false,
                thresholds: Thresholds {
                    mode: ThresholdMode::Fdr,
                    psm: 0.9,
                    peptide: 0.9,
                    protein: 0.9,
                    min_peptides: 1,
                },
            },
        )
        .unwrap();
        drop(project);

        assert!(save_plex(&path, Plex::Tmt6, &[], &[], Certificate::Keep).is_err());
        let channels = vec![
            ("126".into(), "wt-1".into()),
            ("127".into(), "wt-2".into()),
            ("128".into(), "ko-1".into()),
            ("129".into(), "ko-2".into()),
        ];
        // Saved with TMT6 first: 127/128 differ from TMT6-original, so the hint fires.
        let wrong = save_plex(&path, Plex::Tmt6, &channels, &[], Certificate::Keep).unwrap();
        let hint = reporter_reagent_hint(wrong.reporters.as_ref().unwrap()).unwrap();
        assert!(
            hint.contains("TMT6 original") && hint.contains("2 PSMs"),
            "{hint}"
        );

        let right =
            save_plex(&path, Plex::Tmt6Original, &channels, &[], Certificate::Keep).unwrap();
        let quant = right.reporters.as_ref().unwrap();
        assert!(reporter_reagent_hint(quant).is_none());
        assert!(
            reporter_status(quant).starts_with("TMT6 original (plex, 1 run) · 2 protein groups"),
            "{}",
            reporter_status(quant)
        );
        assert!(reporter_status(quant).contains("isotope impurities not corrected"));

        // A certificate is loaded, kept when the plex is saved again without
        // one, and removed.
        let certificate = path.with_extension("certificate.txt");
        std::fs::write(
            &certificate,
            "126 0 0 10 0\n127 0 0 0 0\n128 0 0 0 0\n129 0 0 0 0\n130 0 0 0 0\n131 0 0 0 0\n",
        )
        .unwrap();
        let corrected = |certificate: Certificate| {
            let data = save_plex(&path, Plex::Tmt6Original, &channels, &[], certificate).unwrap();
            reporter_status(data.reporters.as_ref().unwrap())
                .contains("corrected for the lot's isotope impurities")
        };
        assert!(corrected(Certificate::Load(certificate.clone())));
        assert!(corrected(Certificate::Keep));
        assert!(!corrected(Certificate::Remove));
        assert!(!corrected(Certificate::Keep));
        let wrong_kit = save_plex(
            &path,
            Plex::Tmt10,
            &channels,
            &[],
            Certificate::Load(certificate),
        );
        assert!(wrong_kit.is_err());

        let data = save_categories(
            &path,
            &[
                ("wt-1".into(), "wt".into()),
                ("wt-2".into(), "wt".into()),
                ("ko-1".into(), "ko".into()),
                ("ko-2".into(), "ko".into()),
            ],
        )
        .unwrap();
        assert_eq!(data.categories_for(Method::Reporters), ["wt", "ko"]);
        assert_eq!(
            data.default_pair_for(Method::Reporters),
            Some(("wt".into(), "ko".into()))
        );
        let quant = data.reporters.as_ref().unwrap();
        let comparison = quant.quantities.quantities.compare("wt", "ko").unwrap();
        let shown = shown_reporters(&comparison);
        assert_eq!(shown.means, "Mean log2 reporter intensity");
        let all = shown_all(&data, Method::Reporters).unwrap();
        assert_eq!(all.test, "one-way ANOVA on log2 reporter intensities");

        let table = dir.path().join("reporters.tsv");
        export(&table, &Export::Reporters(quant.clone(), Some(comparison))).unwrap();
        let text = std::fs::read_to_string(&table).unwrap();
        assert!(
            text.starts_with("accession\tpsms\tintensity:wt-1 (126)"),
            "{text}"
        );
    }

    #[test]
    fn several_plexes_offer_free_runs_refuse_clashes_and_keep_untouched_settings() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        tuffold_engine::create_demo_project(&path).unwrap();
        let channels_a = vec![
            ("126".into(), "pool-a".into()),
            ("127".into(), "c1".into()),
            ("128".into(), "t1".into()),
        ];
        // The first plex can save with an empty run list (takes every run), or
        // with chosen runs.
        let data = save_plex_named(
            &path,
            " pool-a ",
            true,
            Plex::Tmt6,
            &[],
            &channels_a,
            Some(&["126".into()]),
            Certificate::Keep,
        )
        .unwrap();
        assert_eq!(data.plexes.len(), 1);
        assert_eq!(data.plexes[0].name, "pool-a");
        assert_eq!(
            data.plexes[0].runs,
            ["control-A", "control-B", "treated-A", "treated-B"]
        );
        let data = save_plex_named(
            &path,
            " pool-a ",
            false,
            Plex::Tmt6,
            &["control-A".into(), "control-B".into()],
            &channels_a,
            Some(&["126".into()]),
            Certificate::Keep,
        )
        .unwrap();
        assert_eq!(data.plexes.len(), 1);
        assert_eq!(data.plexes[0].name, "pool-a");
        assert_eq!(data.plexes[0].runs, ["control-A", "control-B"]);
        assert!(empty_runs_take_all(&data.plexes, "pool-a"));
        assert!(!empty_runs_take_all(&data.plexes, "pool-b"));

        // Which runs can be offered: editing pool-a leaves all runs free; a
        // new plex sees pool-a's runs marked with pool-a.
        assert_eq!(
            data.runs_for_plex(Some("pool-a")),
            [
                PlexRun {
                    name: "control-A".into(),
                    in_other_plex: None,
                },
                PlexRun {
                    name: "control-B".into(),
                    in_other_plex: None,
                },
                PlexRun {
                    name: "treated-A".into(),
                    in_other_plex: None,
                },
                PlexRun {
                    name: "treated-B".into(),
                    in_other_plex: None,
                },
            ]
        );
        assert_eq!(
            data.runs_for_plex(None),
            [
                PlexRun {
                    name: "control-A".into(),
                    in_other_plex: Some("pool-a".into()),
                },
                PlexRun {
                    name: "control-B".into(),
                    in_other_plex: Some("pool-a".into()),
                },
                PlexRun {
                    name: "treated-A".into(),
                    in_other_plex: None,
                },
                PlexRun {
                    name: "treated-B".into(),
                    in_other_plex: None,
                },
            ]
        );

        // Name clashes and blank names are refused for a new plex.
        assert_eq!(
            check_new_plex_name("  ", &data.plexes).unwrap_err(),
            "Name the new plex."
        );
        assert_eq!(
            check_new_plex_name(" pool-a ", &data.plexes).unwrap_err(),
            "A plex named 'pool-a' already exists."
        );
        assert_eq!(
            save_plex_named(
                &path,
                "pool-a",
                true,
                Plex::Tmt6,
                &["treated-A".into()],
                &channels_a,
                Some(&["126".into()]),
                Certificate::Keep,
            )
            .unwrap_err(),
            "A plex named 'pool-a' already exists."
        );

        // With another plex present, saving with no runs selected is refused.
        let channels_b = vec![
            ("126".into(), "pool-b".into()),
            ("127".into(), "c2".into()),
            ("128".into(), "t2".into()),
        ];
        let no_runs = save_plex_named(
            &path,
            "pool-b",
            true,
            Plex::Tmt6,
            &[],
            &channels_b,
            Some(&["126".into()]),
            Certificate::Keep,
        )
        .unwrap_err();
        assert!(no_runs.contains("Pick at least one run"), "{no_runs}");

        let data = save_plex_named(
            &path,
            "pool-b",
            true,
            Plex::Tmt6,
            &["treated-A".into(), "treated-B".into()],
            &channels_b,
            Some(&["126".into()]),
            Certificate::Keep,
        )
        .unwrap();
        assert_eq!(data.plexes.len(), 2);

        // Saving pool-a again without touching references (`None`) or the
        // certificate (`Certificate::Keep`) keeps both via kept_references and
        // kept_impurities; switching its reagent set sends the shown
        // references so kept_references does not drop them.
        let certificate = path.with_extension("certificate.txt");
        std::fs::write(
            &certificate,
            "126 0 0 10 0\n127 0 0 0 0\n128 0 0 0 0\n129 0 0 0 0\n130 0 0 0 0\n131 0 0 0 0\n",
        )
        .unwrap();
        let data = save_plex_named(
            &path,
            "pool-a",
            false,
            Plex::Tmt6,
            &["control-A".into(), "control-B".into()],
            &channels_a,
            None,
            Certificate::Load(certificate),
        )
        .unwrap();
        let edited_a = data.plex_named("pool-a");
        let shown_refs = vec!["126".to_string(), "129".to_string()];
        assert_eq!(
            references_for_save(edited_a, Plex::Tmt6, &channels_a, &shown_refs, false),
            None
        );
        assert_eq!(
            references_for_save(
                edited_a,
                Plex::Tmt6Original,
                &channels_a,
                &shown_refs,
                false,
            ),
            Some(vec!["126".to_string()])
        );
        assert_eq!(
            references_for_save(None, Plex::Tmt6, &channels_a, &shown_refs, false),
            Some(vec!["126".to_string()])
        );
        let project = Project::open(&path).unwrap();
        let sent = build_plex_save(
            &project,
            "pool-a",
            false,
            Plex::Tmt6,
            &["control-A".into(), "control-B".into()],
            &channels_a,
            references_for_save(edited_a, Plex::Tmt6, &channels_a, &shown_refs, false).as_deref(),
            Certificate::Keep,
        )
        .unwrap();
        assert!(sent.impurities.is_some());
        assert_eq!(sent.references, ["126"]);
        let switched = build_plex_save(
            &project,
            "pool-a",
            false,
            Plex::Tmt6Original,
            &["control-A".into(), "control-B".into()],
            &channels_a,
            references_for_save(
                edited_a,
                Plex::Tmt6Original,
                &channels_a,
                &shown_refs,
                false,
            )
            .as_deref(),
            Certificate::Keep,
        )
        .unwrap();
        assert_eq!(switched.references, ["126"]);
        drop(project);

        // Removing a plex leaves the other and frees its runs.
        assert!(remove_plex(&path, "missing").is_err());
        let after = remove_plex(&path, "pool-a").unwrap();
        assert_eq!(after.plexes.len(), 1);
        assert_eq!(after.plexes[0].name, "pool-b");
        assert_eq!(after.runs_for_plex(Some("pool-b"))[0].in_other_plex, None);
    }

    fn mbr_report(run: &str, candidates: usize, accepted: usize, trained: bool) -> MbrRunReport {
        MbrRunReport {
            run: run.into(),
            summary: tuffold_store::TransferSummary {
                candidates,
                known: 0,
                trained,
                targets_won: accepted,
                decoys_won: 0,
                accepted,
            },
            aligned: 1,
            donors: 1,
            ms1_scans_damaged: 0,
        }
    }

    #[test]
    fn mbr_summary_totals_runs_and_flags_what_was_off() {
        // Totals by hand: 10 + 5 = 15 looked for, 4 + 0 = 4 transferred.
        let mut second = mbr_report("b", 5, 0, false);
        second.ms1_scans_damaged = 3;
        let text = mbr_summary(&[mbr_report("a", 10, 4, true), second]);
        assert_eq!(
            text,
            "Matched between runs in 2 runs: 4 of 15 precursors looked for were transferred at \
             1% FDR.\na: 4 of 10 transferred\nb: 0 of 5 transferred · ranked by retention-time \
             error (too few known peaks to train) · 3 damaged MS1 scans left out"
        );
    }

    #[test]
    fn mbr_summary_uses_singular_nouns() {
        let mut only = mbr_report("solo", 1, 1, true);
        only.ms1_scans_damaged = 1;
        assert_eq!(
            mbr_summary(&[only]),
            "Matched between runs in 1 run: 1 of 1 precursor looked for was transferred at 1% \
             FDR.\nsolo: 1 of 1 transferred · 1 damaged MS1 scan left out"
        );
    }

    #[test]
    fn matching_a_file_that_is_no_run_fails_loudly_and_still_requantifies() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        tuffold_engine::create_demo_project(&path).unwrap();
        let outcome = match_runs(&path, &[PathBuf::from("/nowhere/not-a-run.mzML")]);
        let message = outcome.result.unwrap_err();
        assert!(
            message.starts_with("Match between runs stopped: ")
                && message.contains("no run named 'not-a-run'"),
            "{message}"
        );
        assert!(outcome.data.is_ok());
    }

    #[test]
    fn matching_an_unvalidated_project_says_to_validate() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.tuff");
        Project::create(&path).unwrap();
        let outcome = match_runs(&path, &[PathBuf::from("/nowhere/r1.mzML")]);
        let message = outcome.result.unwrap_err();
        assert!(message.contains("validate it first"), "{message}");
        assert!(outcome.data.is_err());
    }

    #[test]
    fn a_stop_after_several_files_says_earlier_runs_were_kept() {
        let unknown = MbrError::UnknownRun {
            file: "/d/x.mzML".into(),
            run: "x".into(),
            known: vec!["a".into(), "b".into()],
        };
        assert_eq!(
            stopped_message(&unknown, 2),
            "Match between runs stopped: /d/x.mzML: no run named 'x' in this project (runs: a, \
             b). Runs matched before it were kept."
        );
        // One file: nothing came before it. No validation: nothing had begun.
        assert!(!stopped_message(&unknown, 1).contains("kept"));
        assert_eq!(
            stopped_message(&MbrError::NotValidated, 3),
            "Match between runs stopped: this project has no current validation results; \
             validate it first."
        );
        // A message already ending in a full stop does not get a second.
        assert_eq!(
            stopped_message(&MbrError::Read("cannot read.".into()), 1),
            "Match between runs stopped: cannot read."
        );
    }

    #[test]
    fn a_failed_requantification_after_matching_is_worded_as_such() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.tuff");
        Project::create(&path).unwrap();
        let project = Project::open(&path).unwrap();
        // Matching "succeeded" with no reports, but nothing is validated, so
        // quantifying fails: the summary must not stand alone as a success.
        let outcome = conclude(&project, Ok(Vec::new()), 0);
        assert!(outcome.result.is_ok());
        let problem = outcome.data.unwrap_err();
        assert!(
            problem.starts_with("Matching finished, but quantifying again failed: "),
            "{problem}"
        );
    }

    #[test]
    fn a_successful_match_summarizes_and_requantifies_with_the_transfers() {
        use std::collections::HashSet;
        use tuffold_store::{Transfer, TransferSummary};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        tuffold_engine::create_demo_project(&path).unwrap();
        let mut project = Project::open(&path).unwrap();
        let before = super::quantify(&project).unwrap().intensity.unwrap();
        assert_eq!(before.transferred, 0);

        // The engine writes these from the spectra files; here the first run
        // is handed precursors that other runs identified and it did not.
        let runs = project.mbr_inputs().unwrap();
        let acceptor = &runs[0];
        let mut taken: HashSet<(i64, u8)> = acceptor
            .precursors
            .iter()
            .filter(|p| p.passing)
            .map(|p| (p.peptide_id, p.charge))
            .collect();
        let transfers: Vec<Transfer> = runs[1..]
            .iter()
            .flat_map(|run| &run.precursors)
            .filter(|p| p.passing && taken.insert((p.peptide_id, p.charge)))
            .take(10)
            .map(|p| Transfer {
                peptide_id: p.peptide_id,
                charge: p.charge,
                predicted_rt: p.rt,
                area: Some(1.0e6),
                apex_rt: Some(p.rt),
                q_value: Some(0.001),
            })
            .collect();
        let n = transfers.len();
        assert!(n >= 3, "the demo gives the run {n} precursors to receive");
        let summary = TransferSummary {
            candidates: n,
            known: 0,
            trained: false,
            targets_won: n,
            decoys_won: 0,
            accepted: n,
        };
        project
            .save_transfers(acceptor.id, &transfers, &summary)
            .unwrap();
        let report = MbrRunReport {
            run: acceptor.name.clone(),
            summary,
            aligned: 3,
            donors: 3,
            ms1_scans_damaged: 0,
        };
        let run = acceptor.name.clone();

        let outcome = conclude(&project, Ok(vec![report]), 1);
        assert_eq!(
            outcome.result.unwrap(),
            format!(
                "Matched between runs in 1 run: {n} of {n} precursors looked for were \
                 transferred at 1% FDR.\n{run}: {n} of {n} transferred · ranked by \
                 retention-time error (too few known peaks to train)"
            )
        );
        let intensity = outcome.data.unwrap().intensity.unwrap();
        // Some of them reach protein groups that pass; none can exceed what was saved.
        assert!(
            (1..=n).contains(&intensity.transferred),
            "{} of {n}",
            intensity.transferred
        );
        assert!(ms1_status(&intensity).contains("values matched between runs"));
    }
}
