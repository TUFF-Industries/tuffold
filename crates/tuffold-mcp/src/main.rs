//! `tuffold-mcp`: a Model Context Protocol server (stdio) over the Tuffold
//! engine, so an AI assistant can open `.tuff` projects, browse validated
//! proteins, compare sample groups, and run imports, searches, and
//! validation. A thin layer: each tool calls one function in [`tools`] on a
//! blocking thread and returns its JSON; failures come back as tool errors
//! the assistant can read.

mod privacy;
mod tools;

use std::path::PathBuf;
use std::sync::Arc;

use rmcp::handler::server::router::tool::ToolRouter;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::{CallToolResult, ContentBlock, ServerCapabilities, ServerConfig};
use rmcp::{ErrorData, ServerHandler, ServiceExt, tool, tool_handler, tool_router};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::Value;
use tuffold_engine::{SearchOptions, SearchSettings};
use tuffold_validate::{ThresholdMode, ThresholdRequest};

use crate::privacy::{Level, Privacy};

const INSTRUCTIONS: &str = "Tuffold works on .tuff project files (SQLite). Typical flow: \
    import_data or search_spectra, then validate, then list_protein_groups / protein_group / \
    peptide_spectra to browse, set_categories and compare_categories for quantitation (or \
    set_plex and quantify_reporters for TMT and iTRAQ). Paths are local file paths. Lists are paged \
    (total is always reported). q-values are target-decoy estimates; null means unknown, never \
    zero.";

#[derive(Deserialize, JsonSchema)]
struct ProjectArgs {
    /// Absolute path to the .tuff project file.
    project: PathBuf,
}

#[derive(Deserialize, JsonSchema)]
struct ListArgs {
    project: PathBuf,
    /// Keep groups whose accession or description contains this text.
    filter: Option<String>,
    #[serde(default)]
    offset: usize,
    /// Rows to return (default 50, at most 500).
    limit: Option<usize>,
}

#[derive(Deserialize, JsonSchema)]
struct GroupArgs {
    project: PathBuf,
    group_id: i64,
}

#[derive(Deserialize, JsonSchema)]
struct SpectraArgs {
    project: PathBuf,
    peptide_id: i64,
    /// Spectra to return, best first (default 20, at most 500).
    limit: Option<usize>,
}

#[derive(Deserialize, JsonSchema)]
struct CompareArgs {
    project: PathBuf,
    /// The reference category.
    baseline: String,
    /// The category compared against it (fold changes are versus/baseline).
    versus: String,
    /// `counts` (spectral counts, the default), `intensity` (MS1 precursor
    /// intensity, median polish), or `reporters` (reporter ions; needs a plex).
    method: Option<String>,
    /// With `method` `reporters`: one plex to quantify alone; when omitted,
    /// every plex, bridged through their reference channels when several.
    plex: Option<String>,
    #[serde(default)]
    offset: usize,
    /// Rows to return, most significant first (default 25, at most 500).
    limit: Option<usize>,
}

#[derive(Deserialize, JsonSchema)]
struct ExportArgs {
    project: PathBuf,
    /// The file to write: an absolute path that doesn't exist yet (.xlsx for
    /// the workbook, .msp for the library, .mzid for mzIdentML).
    to: PathBuf,
}

#[derive(Deserialize, JsonSchema)]
struct QcArgs {
    project: PathBuf,
    /// `counts` (spectral counts, the default), `intensity` (MS1, median polish), or
    /// `reporters` (reporter ions; needs a plex).
    method: Option<String>,
    /// With `method` `reporters`: one plex alone; when omitted, every plex,
    /// bridged through their reference channels when several.
    plex: Option<String>,
}

#[derive(Deserialize, JsonSchema)]
struct AllCategoriesArgs {
    project: PathBuf,
    /// `counts` (spectral counts, the default), `intensity` (MS1, median polish), or
    /// `reporters` (reporter ions; needs a plex).
    method: Option<String>,
    /// With `method` `reporters`: one plex to quantify alone; when omitted,
    /// every plex, bridged through their reference channels when several.
    plex: Option<String>,
    #[serde(default)]
    offset: usize,
    /// Rows to return, most significant first (default 25, at most 500).
    limit: Option<usize>,
}

#[derive(Deserialize, JsonSchema)]
struct CategoriesArgs {
    project: PathBuf,
    /// Biosample name -> category; an empty category clears it.
    assignments: std::collections::BTreeMap<String, String>,
}

#[derive(Deserialize, JsonSchema)]
struct PlexArgs {
    project: PathBuf,
    /// Name for the plex (default "plex").
    name: Option<String>,
    /// Reagent set: "tmt6", "tmt6-original", "tmt10", "tmt11", "tmtpro16",
    /// "tmtpro18", "itraq4", or "itraq8".
    label: String,
    /// MS runs in this plex; when omitted or empty, every run in the project
    /// not in another plex.
    #[serde(default)]
    runs: Vec<String>,
    /// Reporter channel -> biosample name (e.g. "126": "control-1", "127N":
    /// "treated-1"). Channel biosamples are created when new.
    channels: std::collections::BTreeMap<String, String>,
    /// Local path to the reagent lot's certificate of isotope impurities
    /// (text: one line per channel of the kit, its name then 4 percentages
    /// for -2/-1/+1/+2 Da or 8 for TMTpro's -2C/-N-C/-C/-N/+N/+C/+N+C/+2C),
    /// or "none" to remove it. When omitted the plex keeps its own.
    impurities: Option<String>,
    /// Channels holding the pooled reference sample that bridges this plex
    /// to the project's others (e.g. ["126"]); [] clears them. When omitted
    /// the plex keeps its own. With several plexes, quantify_reporters
    /// without a plex puts them on one scale through these.
    references: Option<Vec<String>>,
}

#[derive(Deserialize, JsonSchema)]
struct RemovePlexArgs {
    project: PathBuf,
    /// Name of the plex to remove.
    name: String,
}

#[derive(Deserialize, JsonSchema)]
struct ReportersArgs {
    project: PathBuf,
    /// One plex to quantify alone; when omitted, every plex, bridged through
    /// their reference channels when there are several.
    plex: Option<String>,
    /// Optional reference category to compare against `versus`.
    baseline: Option<String>,
    /// Optional category compared against `baseline` (fold changes are
    /// versus/baseline).
    versus: Option<String>,
    #[serde(default)]
    offset: usize,
    /// Rows to return (default 50, or 25 when comparing, at most 500).
    limit: Option<usize>,
}

#[derive(Deserialize, JsonSchema)]
struct MbrArgs {
    project: PathBuf,
    /// The runs' spectra files (.mzML with MS1 scans, or timsTOF .d folders),
    /// named like their runs.
    files: Vec<PathBuf>,
}

#[derive(Deserialize, JsonSchema)]
struct ImportArgs {
    /// The project file; created when it does not exist.
    project: PathBuf,
    /// FragPipe or MaxQuant txt folders, .pep.xml, .mzid, Proteome Discoverer
    /// .pdResult/.msf, Mascot .dat, Scaffold .sf3 projects and .sfdb files, and
    /// .mzML/.mgf spectra files or Bruker timsTOF .d folders.
    inputs: Vec<PathBuf>,
}

#[derive(Deserialize, JsonSchema)]
struct SearchArgs {
    /// The project file; created when it does not exist.
    project: PathBuf,
    /// Protein FASTA (targets; decoys are generated).
    fasta: PathBuf,
    /// .mzML or .mgf files, or Bruker timsTOF .d folders; each becomes its own run and biosample.
    spectra: Vec<PathBuf>,
    /// trypsin (default), trypsin/p, lys-c, lys-n, arg-c, glu-c, asp-n, chymotrypsin.
    enzyme: Option<String>,
    /// Allow one non-enzymatic end.
    #[serde(default)]
    semi_specific: bool,
    /// Default 2.
    missed_cleavages: Option<u8>,
    /// ± ppm, default 20.
    precursor_ppm: Option<f32>,
    /// ± ppm, default 20.
    fragment_ppm: Option<f32>,
    /// ± Da instead of fragment_ppm (low-resolution MS/MS, e.g. 0.5 for an
    /// ion trap); giving both is an error.
    fragment_da: Option<f32>,
    /// Extra fixed modifications as SITE:MOD, e.g. "K:TMT6plex", "C:57.021464".
    /// Sites: a residue letter, ^ $ for peptide termini, [ ] for protein termini.
    #[serde(default)]
    fixed: Vec<String>,
    /// Extra variable modifications as SITE:MOD, e.g. "S:Phospho", "[:Acetyl".
    #[serde(default)]
    variable: Vec<String>,
    /// Isobaric label fixed on K and peptide N-termini: "tmt", "tmtpro",
    /// "itraq4", or "itraq8".
    label: Option<String>,
    /// Drop the default fixed carbamidomethyl C.
    #[serde(default)]
    no_carbamidomethyl: bool,
    /// Drop the default variable oxidized M.
    #[serde(default)]
    no_oxidation: bool,
}

impl SearchArgs {
    fn options(&self) -> SearchOptions {
        SearchOptions {
            precursor_ppm: self.precursor_ppm,
            fragment_ppm: self.fragment_ppm,
            fragment_da: self.fragment_da,
            enzyme: self.enzyme.clone(),
            semi_specific: self.semi_specific,
            missed_cleavages: self.missed_cleavages,
            fixed: self.fixed.clone(),
            variable: self.variable.clone(),
            label: self.label.clone(),
            max_variable_mods: None,
            no_carbamidomethyl: self.no_carbamidomethyl,
            no_oxidation: self.no_oxidation,
            // Raising the budget past what is free is the person's call.
            max_memory_gb: None,
        }
    }
}

#[derive(Deserialize, JsonSchema)]
struct ThresholdArgs {
    project: PathBuf,
    /// Threshold by "fdr" (q-values) or "probability" (1 − PEP). Implied by
    /// which values are given; defaults to the project's current mode.
    by: Option<String>,
    /// PSM FDR, e.g. 0.01 (fdr mode).
    psm_fdr: Option<f64>,
    /// Peptide FDR (fdr mode).
    peptide_fdr: Option<f64>,
    /// Protein group FDR (fdr mode).
    protein_fdr: Option<f64>,
    /// PSM probability, e.g. 0.95 (probability mode).
    psm_probability: Option<f64>,
    /// Peptide probability (probability mode).
    peptide_probability: Option<f64>,
    /// Protein group probability, e.g. 0.99 (probability mode).
    protein_probability: Option<f64>,
    /// Unique or razor peptides a group needs.
    min_peptides: Option<u32>,
}

#[derive(Deserialize, JsonSchema)]
struct ValidateArgs {
    #[serde(flatten)]
    thresholds: ThresholdArgs,
    /// Rescore each engine (default true); false ranks by its best score.
    rescore: Option<bool>,
}

impl ThresholdArgs {
    /// The change asked for; unnamed levels keep the project's values.
    fn request(&self) -> Result<ThresholdRequest, String> {
        let by = match self.by.as_deref() {
            None => None,
            Some(text) => Some(
                ThresholdMode::parse(text)
                    .ok_or_else(|| format!("by takes \"fdr\" or \"probability\", not {text:?}"))?,
            ),
        };
        Ok(ThresholdRequest {
            by,
            fdr: [self.psm_fdr, self.peptide_fdr, self.protein_fdr],
            probability: [
                self.psm_probability,
                self.peptide_probability,
                self.protein_probability,
            ],
            min_peptides: self.min_peptides,
        })
    }
}

#[derive(Clone)]
struct Tuffold {
    tool_router: ToolRouter<Self>,
    privacy: Arc<Privacy>,
}

impl Tuffold {
    /// Run a tool on a blocking thread (they touch SQLite and can take a
    /// while) and turn its outcome into a tool result, with the project's
    /// names hidden in privacy mode. `work` gets the aliases to map names
    /// the assistant sent back to real ones.
    async fn run(
        &self,
        project: PathBuf,
        work: impl FnOnce(&Privacy) -> Result<Value, String> + Send + 'static,
    ) -> CallToolResult {
        let privacy = Arc::clone(&self.privacy);
        let outcome = tokio::task::spawn_blocking(move || match work(&privacy) {
            Ok(value) => Ok(privacy.hide(&project, value)),
            Err(message) => Err(privacy.hide_text(&project, &message)),
        })
        .await
        .unwrap_or_else(|error| Err(format!("the tool stopped unexpectedly: {error}")));
        match outcome.and_then(|value| ContentBlock::json(value).map_err(|e| e.message.to_string()))
        {
            Ok(content) => CallToolResult::success(vec![content]),
            Err(message) => CallToolResult::error(vec![ContentBlock::text(message)]),
        }
    }
}

#[tool_router]
impl Tuffold {
    fn new(privacy: Privacy) -> Self {
        Self {
            tool_router: Self::tool_router(),
            privacy: Arc::new(privacy),
        }
    }

    #[tool(
        description = "Summarize a project: counts, runs, biosamples and their categories, \
                       plexes, imports, and validation (engines, thresholds, passing counts, \
                       notes).",
        annotations(read_only_hint = true)
    )]
    async fn project_summary(
        &self,
        Parameters(args): Parameters<ProjectArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| tools::project_summary(&args.project))
            .await)
    }

    #[tool(
        description = "List protein groups passing the thresholds, best first, with spectral \
                       counts per biosample. Paged; optional text filter.",
        annotations(read_only_hint = true)
    )]
    async fn list_protein_groups(
        &self,
        Parameters(args): Parameters<ListArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                // Hidden accessions are matched exactly, and descriptions,
                // which the assistant never saw, not at all.
                let hides_proteins = privacy.level() == Level::Proteins;
                let filter = args
                    .filter
                    .as_deref()
                    .map(|filter| privacy.reveal(&args.project, filter));
                tools::list_protein_groups(
                    &args.project,
                    filter.as_deref(),
                    !hides_proteins,
                    args.offset,
                    args.limit.unwrap_or(50),
                )
            })
            .await)
    }

    #[tool(
        description = "One protein group: member proteins, its peptides (unique, razor, or \
                       shared, with q-values and spectra), and sequence coverage.",
        annotations(read_only_hint = true)
    )]
    async fn protein_group(
        &self,
        Parameters(args): Parameters<GroupArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::protein_group(&args.project, args.group_id)
            })
            .await)
    }

    #[tool(
        description = "A peptide's validated spectra, best first (run, scan, charge, m/z, \
                       q-value, whether it passes).",
        annotations(read_only_hint = true)
    )]
    async fn peptide_spectra(
        &self,
        Parameters(args): Parameters<SpectraArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::peptide_spectra(&args.project, args.peptide_id, args.limit.unwrap_or(20))
            })
            .await)
    }

    #[tool(
        description = "Compare two biosample categories: fold changes and q-values, most \
                       significant first. By spectral counts (Welch's t-test with replicates, \
                       else Fisher's exact test), by MS1 intensity (method \"intensity\": \
                       median-polish protein intensities, Welch's t-test), or by TMT reporter ions \
                       (method \"reporters\": median-normalized channel intensities, Welch's \
                       t-test).",
        annotations(read_only_hint = true)
    )]
    async fn compare_categories(
        &self,
        Parameters(args): Parameters<CompareArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let limit = args.limit.unwrap_or(25);
                let baseline = privacy.reveal(&args.project, &args.baseline);
                let versus = privacy.reveal(&args.project, &args.versus);
                let plex = args
                    .plex
                    .as_deref()
                    .map(|plex| privacy.reveal(&args.project, plex));
                match args.method.as_deref().unwrap_or("counts") {
                    "counts" => tools::compare_categories(
                        &args.project,
                        &baseline,
                        &versus,
                        args.offset,
                        limit,
                    ),
                    "intensity" => tools::compare_intensities(
                        &args.project,
                        &baseline,
                        &versus,
                        args.offset,
                        limit,
                    ),
                    "reporters" => tools::compare_reporters(
                        &args.project,
                        plex.as_deref(),
                        &baseline,
                        &versus,
                        args.offset,
                        limit,
                    ),
                    other => Err(format!(
                        "method is counts, intensity, or reporters, not '{other}'"
                    )),
                }
            })
            .await)
    }

    #[tool(
        description = "Compare every biosample category at once with a one-way ANOVA per \
                       protein group (spectral counts, MS1 intensity with method \"intensity\", \
                       or reporter ions with method \"reporters\"): per-category means, F, \
                       and q-values, most significant first.",
        annotations(read_only_hint = true)
    )]
    async fn compare_all_categories(
        &self,
        Parameters(args): Parameters<AllCategoriesArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let plex = args
                    .plex
                    .as_deref()
                    .map(|plex| privacy.reveal(&args.project, plex));
                tools::compare_all_categories(
                    &args.project,
                    args.method.as_deref().unwrap_or("counts"),
                    plex.as_deref(),
                    args.offset,
                    args.limit.unwrap_or(25),
                )
            })
            .await)
    }

    #[tool(
        description = "Sample QC by counts, MS1 \"intensity\", or \"reporters\": groups measured \
                       and median per sample, Pearson correlation of every pair, and a PCA over \
                       the groups measured in every sample. For spotting failed runs, swapped \
                       labels, and batch effects.",
        annotations(read_only_hint = true)
    )]
    async fn sample_qc(
        &self,
        Parameters(args): Parameters<QcArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let plex = args
                    .plex
                    .as_deref()
                    .map(|plex| privacy.reveal(&args.project, plex));
                tools::sample_qc(
                    &args.project,
                    args.method.as_deref().unwrap_or("counts"),
                    plex.as_deref(),
                )
            })
            .await)
    }

    #[tool(
        description = "Assign biosamples to categories (e.g. control, treated) for comparisons."
    )]
    async fn set_categories(
        &self,
        Parameters(args): Parameters<CategoriesArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let pairs: Vec<(String, String)> = args
                    .assignments
                    .iter()
                    .map(|(sample, category)| {
                        (
                            privacy.reveal(&args.project, sample),
                            privacy.reveal(&args.project, category),
                        )
                    })
                    .collect();
                tools::set_categories(&args.project, &pairs)
            })
            .await)
    }

    #[tool(
        description = "Define (or redefine) an isobaric plex: its reagent set (tmt6, \
                       tmt6-original, tmt10, tmt11, tmtpro16, tmtpro18, itraq4, itraq8), its MS runs (all \
                       unassigned runs when omitted), which biosample each reporter \
                       channel holds, and optionally the lot's isotope impurity certificate \
                       (reporters are then corrected for it) and the reference channels that \
                       bridge it to other plexes."
    )]
    async fn set_plex(
        &self,
        Parameters(args): Parameters<PlexArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let name = privacy.reveal(&args.project, args.name.as_deref().unwrap_or("plex"));
                let runs: Vec<String> = args
                    .runs
                    .iter()
                    .map(|run| privacy.reveal(&args.project, run))
                    .collect();
                let channels: Vec<(String, String)> = args
                    .channels
                    .iter()
                    .map(|(channel, sample)| {
                        (channel.clone(), privacy.reveal(&args.project, sample))
                    })
                    .collect();
                tools::set_plex(
                    &args.project,
                    &name,
                    &args.label,
                    &runs,
                    &channels,
                    args.impurities.as_deref(),
                    args.references.as_deref(),
                )
            })
            .await)
    }

    #[tool(
        description = "Remove a named isobaric plex (its runs are free for another plex \
                       again; its channel biosamples stay)."
    )]
    async fn remove_plex(
        &self,
        Parameters(args): Parameters<RemovePlexArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let name = privacy.reveal(&args.project, &args.name);
                tools::remove_plex(&args.project, &name)
            })
            .await)
    }

    #[tool(
        description = "Quantify isobaric plexes by reporter ions: median-normalized channel \
                       intensities and log2 values per passing protein group, channel loading \
                       shifts, and any reagent set that would match more PSMs. Several plexes \
                       are put on one scale through their reference channels (internal \
                       reference scaling), which are then left out. Optionally compare two \
                       categories (baseline and versus) by Welch's t-test. Paged.",
        annotations(read_only_hint = true)
    )]
    async fn quantify_reporters(
        &self,
        Parameters(args): Parameters<ReportersArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |privacy| {
                let plex = args
                    .plex
                    .as_deref()
                    .map(|plex| privacy.reveal(&args.project, plex));
                let baseline = args
                    .baseline
                    .as_deref()
                    .map(|b| privacy.reveal(&args.project, b));
                let versus = args
                    .versus
                    .as_deref()
                    .map(|v| privacy.reveal(&args.project, v));
                let default_limit = if baseline.is_some() || versus.is_some() {
                    25
                } else {
                    50
                };
                tools::quantify_reporters(
                    &args.project,
                    plex.as_deref(),
                    baseline.as_deref(),
                    versus.as_deref(),
                    args.offset,
                    args.limit.unwrap_or(default_limit),
                )
            })
            .await)
    }

    #[tool(
        description = "Match between runs: in each given run, measure the MS1 peak of precursors \
                       identified only in other runs, where they should elute (aligned retention \
                       time), keeping transfers at 1% FDR by target-decoy competition (decoys \
                       looked for at the wrong time). Needs a validated project and the runs' \
                       mzML files; MS1 intensity quantitation then uses the accepted transfers."
    )]
    async fn match_between_runs(
        &self,
        Parameters(args): Parameters<MbrArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::match_between_runs(&args.project, &args.files)
            })
            .await)
    }

    #[tool(
        description = "Import search results (FragPipe folders, pepXML, mzIdentML) and spectra \
                       (.mzML, .mgf, Bruker .d) into a project, creating it if needed. Importing clears \
                       validation results: validate afterwards. An error when nothing imports."
    )]
    async fn import_data(
        &self,
        Parameters(args): Parameters<ImportArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::import_data(&args.project, &args.inputs)
            })
            .await)
    }

    #[tool(
        description = "Search .mzML/.mgf spectra or Bruker .d folders against a FASTA with the built-in search \
                       (Sage's core) and import the matches. Defaults: trypsin, 2 missed \
                       cleavages, carbamidomethyl C fixed, oxidized M variable, 20 ppm; enzyme, \
                       TMT or iTRAQ labels, other modifications, and a fragment tolerance in Da \
                       (low-resolution MS/MS) can be set. Validate afterwards. An \
                       error when nothing matches."
    )]
    async fn search_spectra(
        &self,
        Parameters(args): Parameters<SearchArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                let settings = SearchSettings::from_options(&args.options())?;
                tools::search_spectra(&args.project, &args.fasta, &args.spectra, settings)
            })
            .await)
    }

    #[tool(
        description = "Validate a project: rank each engine's PSMs (rescored by default), let \
                       engines vote, and estimate FDR for spectra, peptides, and protein groups."
    )]
    async fn validate(
        &self,
        Parameters(args): Parameters<ValidateArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.thresholds.project.clone();
        Ok(self
            .run(project, move |_| {
                let request = args.thresholds.request()?;
                tools::validate(
                    &args.thresholds.project,
                    &request,
                    args.rescore.unwrap_or(true),
                )
            })
            .await)
    }

    #[tool(
        description = "Apply new thresholds to a validated project without validating again: by \
                       FDR (psm_fdr, ...) or by probability (psm_probability, ...; 1 − PEP, e.g. \
                       0.95 peptide, 0.99 protein). Unnamed levels keep their value; switching \
                       mode starts from that mode's defaults. Probability mode reports the \
                       decoy-estimated FDR of what passes."
    )]
    async fn set_thresholds(
        &self,
        Parameters(args): Parameters<ThresholdArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::set_thresholds(&args.project, &args.request()?)
            })
            .await)
    }

    #[tool(
        description = "Undo the last import, search, or sample-group change (the project is \
                       copied before each). Calling it again redoes the change."
    )]
    async fn undo_last_change(
        &self,
        Parameters(args): Parameters<ProjectArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| tools::undo_last_change(&args.project))
            .await)
    }

    #[tool(
        description = "Export a validated project as an Excel workbook (.xlsx; absolute path, \
                       must not exist): summary, methods, proteins, peptides, spectra, and ANOVA \
                       tables at the saved thresholds."
    )]
    async fn export_workbook(
        &self,
        Parameters(args): Parameters<ExportArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::export_workbook(&args.project, &args.to)
            })
            .await)
    }

    #[tool(
        description = "Export a validated project's results as mzIdentML 1.2 (.mzid; absolute \
                       path, must not exist): every validated spectrum and protein group, target \
                       and decoy, with q-values, PEPs, and pass flags, e.g. for PRIDE."
    )]
    async fn export_mzid(
        &self,
        Parameters(args): Parameters<ExportArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::export_mzid(&args.project, &args.to)
            })
            .await)
    }

    #[tool(
        description = "Export a validated project's spectral library (NIST .msp; absolute path, \
                       must not exist): the best passing spectrum of every peptide ion, with b/y \
                       annotations, readable by SpectraST, Skyline, and NIST MS Search."
    )]
    async fn export_library(
        &self,
        Parameters(args): Parameters<ExportArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| {
                tools::export_library(&args.project, &args.to)
            })
            .await)
    }

    #[tool(
        description = "A methods paragraph for publication, written from what the project did \
                       (engines, rescoring, thresholds, results, quantitation), with references.",
        annotations(read_only_hint = true)
    )]
    async fn methods_text(
        &self,
        Parameters(args): Parameters<ProjectArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| tools::methods_text(&args.project))
            .await)
    }

    #[tool(description = "Create the synthetic demo project (labelled synthetic) to explore with.")]
    async fn create_demo(
        &self,
        Parameters(args): Parameters<ProjectArgs>,
    ) -> Result<CallToolResult, ErrorData> {
        let project = args.project.clone();
        Ok(self
            .run(project, move |_| tools::create_demo(&args.project))
            .await)
    }
}

#[tool_handler(router = self.tool_router)]
impl ServerHandler for Tuffold {
    fn get_info(&self) -> ServerConfig {
        let instructions = match self.privacy.level().instructions() {
            Some(privacy) => format!("{INSTRUCTIONS} {privacy}"),
            None => INSTRUCTIONS.to_string(),
        };
        ServerConfig::new(ServerCapabilities::builder().enable_tools().build())
            .with_instructions(instructions)
    }
}

#[tokio::main]
async fn main() {
    tuffold_engine::install_panic_log("tuffold-mcp");
    let args: Vec<String> = std::env::args().skip(1).collect();
    let level = match Level::parse(&args, std::env::var(privacy::ENV).ok().as_deref()) {
        Ok(level) => level,
        Err(message) => {
            // stdout carries the protocol; diagnostics go to stderr.
            eprintln!("tuffold-mcp: {message}");
            std::process::exit(2);
        }
    };
    let served = Tuffold::new(Privacy::new(level))
        .serve(rmcp::transport::stdio())
        .await;
    let result = match served {
        Ok(service) => service
            .waiting()
            .await
            .map(|_| ())
            .map_err(|e| e.to_string()),
        Err(error) => Err(error.to_string()),
    };
    if let Err(message) = result {
        // stdout carries the protocol; diagnostics go to stderr.
        eprintln!("tuffold-mcp: {message}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(by: Option<&str>, protein_probability: Option<f64>) -> ThresholdArgs {
        ThresholdArgs {
            project: PathBuf::from("/p.tuff"),
            by: by.map(String::from),
            psm_fdr: None,
            peptide_fdr: None,
            protein_fdr: None,
            psm_probability: None,
            peptide_probability: None,
            protein_probability,
            min_peptides: None,
        }
    }

    #[test]
    fn threshold_arguments_become_a_request() {
        let request = args(Some("probability"), Some(0.9)).request().unwrap();
        assert_eq!(request.by, Some(ThresholdMode::Probability));
        assert_eq!(request.probability, [None, None, Some(0.9)]);
        let error = args(Some("vibes"), None).request().unwrap_err();
        assert!(error.contains("probability"), "{error}");
    }
}
