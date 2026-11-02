//! `tuffold` — command-line access to the Tuffold engine.

/// `println!` that stops quietly when the reader has gone (`tuffold info x |
/// head`), as command-line tools do, instead of panicking; any other write
/// error is reported and fails the command.
macro_rules! out {
    ($($arg:tt)*) => {
        $crate::write_out(format_args!("{}\n", format_args!($($arg)*)))
    };
}

/// `print!` likewise.
macro_rules! out_raw {
    ($($arg:tt)*) => {
        $crate::write_out(format_args!($($arg)*))
    };
}

fn write_out(text: std::fmt::Arguments<'_>) {
    use std::io::Write as _;
    if let Err(error) = std::io::stdout().lock().write_fmt(text) {
        if error.kind() == std::io::ErrorKind::BrokenPipe {
            std::process::exit(0);
        }
        eprintln!("error: writing output: {error}");
        std::process::exit(1);
    }
}

mod quant;
mod search;
mod validation;

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};
use tuffold_engine::demo::DEMO_SEED;
use tuffold_engine::{
    InputKind, OutcomeKind, classify, create_demo_project, duplicate_note, import_inputs,
};
use tuffold_io::{BundleInspection, inspect_bundle};
use tuffold_store::{Project, ProjectSummary};

#[derive(Parser)]
#[command(
    name = "tuffold",
    version,
    about = "Tuffold: open proteomics validation and viewing"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Check a FragPipe results folder and report every problem found.
    Inspect {
        bundle: PathBuf,
        /// Print the inspection as JSON.
        #[arg(long)]
        json: bool,
    },
    /// Import search results into a project file, creating it if needed.
    ///
    /// Accepts FragPipe results folders, MaxQuant txt folders, pepXML (.pep.xml/.pepXML),
    /// mzIdentML (.mzid), optionally gzip-compressed, Proteome Discoverer
    /// results (.pdResult, .msf) and Scaffold files (.sf3, .sfdb), one import per
    /// search engine, plus spectra files
    /// (.mzML, .mgf, Bruker timsTOF .d folders) whose peaks attach to the identified spectra.
    /// Identifications are imported first, then spectra.
    Import {
        /// Files or folders to import.
        #[arg(required = true)]
        inputs: Vec<PathBuf>,
        /// The project file (.tuff) to import into.
        #[arg(long, short = 'o')]
        into: PathBuf,
        /// The MS run a single spectra file belongs to, when its file name
        /// differs from the run name.
        #[arg(long)]
        run: Option<String>,
    },
    /// Search spectra files (.mzML, .mgf, Bruker timsTOF .d) against a FASTA with the built-in
    /// search (Sage's search core), import the PSMs into a project (creating
    /// it if needed), and attach the peaks. Validate the project afterwards.
    Search {
        /// Spectra files to search; each is its own run and biosample.
        #[arg(required = true)]
        spectra: Vec<PathBuf>,
        /// The protein database (targets; decoys are generated).
        #[arg(long)]
        fasta: PathBuf,
        /// The project file (.tuff) to import into.
        #[arg(long, short = 'o')]
        into: PathBuf,
        #[command(flatten)]
        settings: search::SearchArgs,
    },
    /// Summarize a project file.
    Info { project: PathBuf },
    /// Undo the last import, search, or sample-group change (run it again to
    /// redo). The copy it restores is kept next to the project as
    /// PROJECT.undo.
    Undo { project: PathBuf },
    /// Validate a project: rank each engine's PSMs (rescoring them by
    /// default), let engines vote per spectrum, estimate FDR for spectra,
    /// peptides, and protein groups, and store the results.
    Validate {
        project: PathBuf,
        /// Rank each engine by its best single score instead of rescoring.
        #[arg(long)]
        no_rescore: bool,
        /// Rescore without the deep-learning models even when installed.
        #[arg(long)]
        no_deep: bool,
        #[command(flatten)]
        thresholds: validation::ThresholdArgs,
    },
    /// Create a demo project from synthetic data (labelled as such), to
    /// explore Tuffold without lab data.
    Demo {
        /// The project file (.tuff) to create.
        #[arg(long, short = 'o')]
        into: PathBuf,
    },
    /// List biosamples and their categories; assign with NAME=CATEGORY.
    Samples {
        project: PathBuf,
        /// `NAME=CATEGORY` assignments (`NAME=` clears a category).
        #[arg(long = "set")]
        assignments: Vec<String>,
    },
    /// Define an isobaric plex (TMT, TMTpro, iTRAQ): its reagent set, its runs (one
    /// run, or the fractions of one labelled pool), and the sample in each
    /// channel; `--remove` removes a named plex; without `--label`, list the
    /// project's plexes.
    Plex {
        project: PathBuf,
        /// The plex's name (one plex per labelled pool).
        name: Option<String>,
        /// tmt6, tmt6-original (127 and 129 heavier, kits until about 2012), tmt10, tmt11, tmtpro16, tmtpro18, itraq4, or itraq8.
        #[arg(long, conflicts_with = "remove")]
        label: Option<String>,
        /// Remove the named plex (its channel biosamples stay).
        #[arg(long)]
        remove: bool,
        /// A run of the plex (repeatable); default: every run not in another
        /// plex.
        #[arg(long = "run", conflicts_with = "remove")]
        runs: Vec<String>,
        /// The reagent lot's certificate of isotope impurities: a text file,
        /// one line per channel of the kit with its percentages, 4 (-2, -1,
        /// +1, +2 Da) or 8 (TMTpro: -2C, -N-C, -C, -N, +N, +C, +N+C, +2C);
        /// `none` removes it. Default: the plex's own, if it has one.
        #[arg(long, conflicts_with = "remove")]
        impurities: Option<String>,
        /// A channel holding the pooled reference sample that bridges this
        /// plex to the others (repeatable; `none` clears). Default: the
        /// plex's own.
        #[arg(long = "reference", conflicts_with = "remove")]
        references: Vec<String>,
        /// `CHANNEL=SAMPLE` assignments (e.g. 126=control-1); channels left
        /// out are not quantified.
        #[arg(conflicts_with = "remove")]
        channels: Vec<String>,
    },
    /// Quantitation, by spectral counts (normalized counts and NSAF; Fisher's
    /// exact test, and Welch's t-test with replicates) or by MS1 precursor
    /// intensity (median polish; Welch's t-test), per biosample, optionally
    /// comparing two categories (Benjamini-Hochberg q-values).
    Quant {
        project: PathBuf,
        /// `counts` (spectral counts), `intensity` (MS1 precursor
        /// intensity; needs spectra files with MS1 scans), or `reporters`
        /// (reporter ions of a plex; see `tuffold plex`).
        #[arg(long, default_value = "counts", value_name = "METHOD")]
        by: String,
        /// With --by reporters: the plex to quantify (default: every plex,
        /// bridged by their reference channels when there are several).
        #[arg(long)]
        plex: Option<String>,
        /// Compare category B against A.
        #[arg(long, num_args = 2, value_names = ["A", "B"], conflicts_with_all = ["anova", "qc"])]
        compare: Option<Vec<String>>,
        /// Compare every category at once (one-way ANOVA); `--export` then
        /// writes the ANOVA table.
        #[arg(long, conflicts_with = "qc")]
        anova: bool,
        /// Sample QC instead of a comparison: groups measured, median, PCA,
        /// and Pearson correlation of every pair of samples; `--export` then
        /// writes them.
        #[arg(long)]
        qc: bool,
        /// Write a tab-separated table.
        #[arg(long)]
        export: Option<PathBuf>,
        /// With --by reporters: leave out spectra whose precursor made up
        /// less than this share of the isolation window (0–1; `none` keeps
        /// every spectrum). Default 0.5.
        #[arg(long, value_name = "FRACTION")]
        min_purity: Option<String>,
    },
    /// Match between runs: measure precursors identified in other runs where
    /// they should elute in each given run, keeping transfers at 1%
    /// target-decoy FDR (decoys looked for at the wrong time). Needs a
    /// validated project and the runs' spectra files (mzML with MS1 scans, or
    /// timsTOF .d folders); MS1 intensities then use the accepted transfers.
    Mbr {
        project: PathBuf,
        /// The runs' spectra files (named like their runs).
        #[arg(required = true)]
        files: Vec<PathBuf>,
    },
    /// Change the thresholds of a validated project without validating again.
    Thresholds {
        project: PathBuf,
        #[command(flatten)]
        thresholds: validation::ThresholdArgs,
    },
    /// Summarize PTM site localization (probabilities, Ascores, and decoy-AA
    /// false localization rate) for passing PSMs with ambiguous variable
    /// modifications, and optionally export their site scores as a TSV table.
    Localize {
        project: PathBuf,
        /// Site probability cut, as a percentage ("95%") or fraction ("0.95").
        #[arg(long, value_parser = validation::parse_rate)]
        site_prob: Option<f64>,
        /// Write a tab-separated table of localized PSMs.
        #[arg(long)]
        export: Option<PathBuf>,
    },
    /// Export a validated project, by the file's extension: .xlsx (Excel
    /// workbook: summary, methods, proteins, peptides, spectra, ANOVA), .mzid
    /// (mzIdentML 1.2 results, e.g. for PRIDE), or .msp (spectral library for
    /// SpectraST, Skyline, NIST MS Search).
    Export {
        project: PathBuf,
        /// The file to write (.xlsx, .mzid, or .msp); it must not exist.
        #[arg(long = "to", value_name = "FILE")]
        to: PathBuf,
    },
    /// Add GO annotations: the proteins' annotations rolled up to a GO
    /// slim's terms, for filtering (replaces an earlier GO import).
    Go {
        project: PathBuf,
        /// The GO ontology (go-basic.obo).
        #[arg(long, value_name = "OBO")]
        ontology: PathBuf,
        /// A GO slim (goslim_generic.obo).
        #[arg(long, value_name = "OBO")]
        slim: PathBuf,
        /// Annotation files keyed by UniProt accession (goa_human.gaf.gz).
        /// Repeatable.
        #[arg(long = "annotations", value_name = "GAF", required = true)]
        annotations: Vec<PathBuf>,
    },
    /// Print a methods paragraph for publication, written from what the
    /// project actually did (engines, rescoring, thresholds, results), with
    /// references.
    Methods { project: PathBuf },
    /// Deep-learning models for rescoring (AlphaPeptDeep's fragment-intensity
    /// and retention-time models): show whether they are installed, or
    /// install them from AlphaPeptDeep's model archive. Once installed,
    /// validation uses them whenever it rescores.
    Models {
        /// AlphaPeptDeep's model archive (pretrained_models_v3.zip), downloaded
        /// from its GitHub releases.
        #[arg(long, value_name = "ZIP")]
        install: Option<PathBuf>,
    },
    /// A report to send with a bug: version, system, models, a project's
    /// counts, and the last crash logs. Holds no sample, run, or protein names.
    Diagnostics {
        /// A project to include the counts of.
        project: Option<PathBuf>,
        /// Write it here instead of printing it.
        #[arg(long, value_name = "FILE")]
        to: Option<PathBuf>,
    },
}

fn main() -> ExitCode {
    tuffold_engine::install_panic_log("tuffold");
    let result = match Cli::parse().command {
        Command::Inspect { bundle, json } => return inspect(&bundle, json),
        Command::Import { inputs, into, run } => import(&inputs, &into, run.as_deref()),
        Command::Search {
            spectra,
            fasta,
            into,
            settings,
        } => search::search(&spectra, &fasta, &into, &settings),
        Command::Info { project } => info(&project),
        Command::Undo { project } => undo(&project),
        Command::Validate {
            project,
            no_rescore,
            no_deep,
            thresholds,
        } => validation::validate(&project, no_rescore, no_deep, &thresholds),
        Command::Thresholds {
            project,
            thresholds,
        } => validation::thresholds(&project, &thresholds),
        Command::Localize {
            project,
            site_prob,
            export,
        } => validation::localize(&project, site_prob, export.as_deref()),
        Command::Demo { into } => demo(&into),
        Command::Export { project, to } => export(&project, &to),
        Command::Go {
            project,
            ontology,
            slim,
            annotations,
        } => go(
            &project,
            &tuffold_engine::GoFiles {
                ontology,
                slim,
                annotations,
            },
        ),
        Command::Models { install } => models(install.as_deref()),
        Command::Diagnostics { project, to } => diagnostics(project.as_deref(), to.as_deref()),
        Command::Methods { project } => Project::open(&project)
            .map_err(|e| e.to_string())
            .and_then(|p| tuffold_engine::methods(&p).map_err(|e| e.to_string()))
            .map(|m| out_raw!("{}", m.text())),
        Command::Samples {
            project,
            assignments,
        } => quant::samples(&project, &assignments),
        Command::Plex {
            project,
            name,
            label,
            remove,
            runs,
            impurities,
            references,
            channels,
        } => quant::plex(
            &project,
            name.as_deref(),
            label.as_deref(),
            remove,
            &runs,
            impurities.as_deref(),
            &references,
            &channels,
        ),
        Command::Mbr { project, files } => mbr(&project, &files),
        Command::Quant {
            project,
            by,
            plex,
            compare,
            anova,
            qc,
            export,
            min_purity,
        } => {
            let given = min_purity.is_some();
            quant::min_purity(min_purity.as_deref()).and_then(|purity| match (by.as_str(), anova) {
                _ if given && by != "reporters" => {
                    Err("--min-purity applies to --by reporters only".into())
                }
                ("counts" | "intensity" | "reporters", _) if qc => {
                    quant::qc(&project, &by, plex.as_deref(), purity, export.as_deref())
                }
                ("counts" | "intensity" | "reporters", true) => {
                    quant::anova(&project, &by, plex.as_deref(), purity, export.as_deref())
                }
                ("reporters", false) => quant::reporters(
                    &project,
                    plex.as_deref(),
                    purity,
                    compare.as_deref(),
                    export.as_deref(),
                ),
                ("counts", false) => quant::quant(&project, compare.as_deref(), export.as_deref()),
                ("intensity", false) => {
                    quant::intensity(&project, compare.as_deref(), export.as_deref())
                }
                (other, _) => Err(format!(
                    "--by takes counts, intensity, or reporters, not '{other}'"
                )),
            })
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("error: {message}");
            ExitCode::FAILURE
        }
    }
}

fn undo(project: &Path) -> Result<(), String> {
    let note = tuffold_engine::undo(project).map_err(|error| error.to_string())?;
    let (done, again) = if note.redo {
        ("redid", "undo")
    } else {
        ("undid", "redo")
    };
    out!(
        "{done}: {} (run tuffold undo again to {again} it)",
        note.label
    );
    Ok(())
}

fn inspect(bundle: &Path, json: bool) -> ExitCode {
    let inspection = inspect_bundle(bundle);
    let printed = if json {
        match serde_json::to_string_pretty(&inspection) {
            Ok(rendered) => {
                out!("{rendered}");
                true
            }
            Err(error) => {
                eprintln!("error: could not serialize output: {error}");
                false
            }
        }
    } else {
        print_inspection(&inspection);
        true
    };
    if printed && inspection.ok() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

fn import(inputs: &[PathBuf], into: &Path, run: Option<&str>) -> Result<(), String> {
    // Reject unrecognized inputs before creating or touching the project.
    if let Some(unknown) = inputs.iter().find(|input| classify(input).is_none()) {
        if tuffold_io::is_thermo_raw(unknown) {
            return Err(format!("{}: {}", unknown.display(), tuffold_io::THERMO_RAW));
        }
        return Err(format!(
            "{}: unrecognized input; expected a FragPipe results folder, a .pep.xml/.pepXML \
             file, an .mzid file, a MaxQuant txt folder, a .dat, .pdResult/.msf, .sf3 or .sfdb file, or a \
             .mzML/.mgf spectra file or Bruker .d folder",
            unknown.display()
        ));
    }
    let spectra_inputs = inputs
        .iter()
        .filter(|input| classify(input) == Some(InputKind::Spectra))
        .count();
    if run.is_some() && spectra_inputs != 1 {
        return Err("--run names the run of exactly one spectra file".into());
    }
    let mut project = if into.exists() {
        Project::open(into)
    } else {
        Project::create(into)
    }
    .map_err(|error| error.to_string())?;

    let outcomes = import_inputs(&mut project, inputs, run);
    for outcome in &outcomes {
        let input = outcome.input.display();
        match &outcome.kind {
            OutcomeKind::Imported {
                format,
                software,
                report,
                warnings,
                notes,
            } => {
                for issue in warnings {
                    out!("warning [{}] {}", issue.code, issue.message);
                }
                for note in notes {
                    out!("note: {input}: {note}");
                }
                if let Some(note) = duplicate_note(report) {
                    out!("warning [duplicates-kept-once] {input}: {note}");
                }
                out!(
                    "imported {input} ({}, {software}): {} PSMs; new: {} runs, {} spectra, \
                     {} peptides, {} proteins",
                    format.as_str(),
                    report.psms,
                    report.new_ms_runs,
                    report.new_spectra,
                    report.new_peptides,
                    report.new_proteins,
                );
            }
            OutcomeKind::Attached(report) => {
                out!(
                    "attached {input} to run {}: {} of {} MS/MS spectra matched identified \
                     spectra; {} identified spectra still without peaks; {}",
                    report.run,
                    report.attached,
                    report.read,
                    report.identified_without_peaks,
                    tuffold_engine::ms1_note(report),
                );
                if let Some(note) = tuffold_engine::ms3_note(report) {
                    out!("  {note}");
                }
            }
            OutcomeKind::Failed(message) => eprintln!("error: {input}: {message}"),
        }
    }
    print_summary(&project.summary().map_err(|error| error.to_string())?);
    let failures = outcomes
        .iter()
        .filter(|outcome| outcome.is_failure())
        .count();
    if failures > 0 {
        return Err(format!(
            "{failures} of {} inputs failed to import",
            outcomes.len()
        ));
    }
    Ok(())
}

fn info(path: &Path) -> Result<(), String> {
    let project = Project::open(path).map_err(|error| error.to_string())?;
    for record in project.imports().map_err(|error| error.to_string())? {
        let version = record
            .software_version
            .as_deref()
            .unwrap_or("version unknown");
        out!(
            "import {}: {} [{}] {} ({}), {} PSMs",
            record.id,
            record.source_path,
            record.format,
            record.software,
            version,
            record.psms
        );
    }
    print_summary(&project.summary().map_err(|error| error.to_string())?);
    match project
        .validation_summary()
        .map_err(|error| error.to_string())?
    {
        Some(summary) => validation::print_validation(&summary),
        None => out!("validation: not yet run; run `tuffold validate`"),
    }
    Ok(())
}

fn demo(into: &Path) -> Result<(), String> {
    let summary = create_demo_project(into).map_err(|error| error.to_string())?;
    out!(
        "created {} from synthetic data (seed {DEMO_SEED}); nothing in it is measured",
        into.display()
    );
    validation::print_validation(&summary);
    Ok(())
}

fn print_inspection(inspection: &BundleInspection) {
    out!("bundle: {}", inspection.path.display());
    for (file, rows) in &inspection.row_counts {
        out!("  {file:<24} {rows} rows");
    }
    out!(
        "  experiment folders: {}",
        inspection.experiment_folders.join(", ")
    );
    match &inspection.selected_fasta {
        Some(fasta) => out!(
            "  fasta: {} ({} records)",
            fasta.display(),
            inspection.fasta_record_count.unwrap_or_default()
        ),
        None => out!("  fasta: none"),
    }
    for issue in &inspection.warnings {
        out!("warning [{}] {}", issue.code, issue.message);
    }
    for issue in &inspection.errors {
        out!("error [{}] {}", issue.code, issue.message);
    }
    out!("{}", if inspection.ok() { "OK" } else { "INVALID" });
}

fn print_summary(summary: &ProjectSummary) {
    out!(
        "project: {} imports, {} biosamples, {} runs, {} spectra, {} PSMs ({} decoy), \
         {} peptides, {} proteins ({} decoy, {} with sequence)",
        summary.imports,
        summary.biosamples,
        summary.ms_runs,
        summary.spectra,
        summary.psms,
        summary.decoy_psms,
        summary.peptides,
        summary.proteins,
        summary.decoy_proteins,
        summary.proteins_with_sequence,
    );
    for run in &summary.runs {
        out!(
            "  run {} (biosample {}): {} spectra, {} PSMs",
            run.name,
            run.biosample,
            run.spectra,
            run.psms
        );
    }
}

fn diagnostics(project: Option<&Path>, to: Option<&Path>) -> Result<(), String> {
    let report = tuffold_engine::diagnostics_report(project);
    match to {
        None => out_raw!("{report}"),
        Some(to) => {
            std::fs::write(to, report).map_err(|e| format!("{}: {e}", to.display()))?;
            out!("wrote {}; read it before sending it", to.display());
        }
    }
    Ok(())
}

fn models(install: Option<&Path>) -> Result<(), String> {
    if let Some(archive) = install {
        let dir = tuffold_engine::install_models(archive)?;
        out!(
            "installed AlphaPeptDeep's MS2 and RT models in {}",
            dir.display()
        );
        return Ok(());
    }
    let status = tuffold_engine::models_status();
    let place = status
        .dir
        .as_ref()
        .map_or_else(|| "nowhere".to_string(), |dir| dir.display().to_string());
    if status.installed {
        out!("deep-learning models installed in {place}; validation uses them when rescoring");
    } else if let Some(problem) = status.problem {
        return Err(format!("the models in {place} can't be used: {problem}"));
    } else {
        out!(
            "no deep-learning models installed (looked in {place}). Download {} and run \
             `tuffold models --install pretrained_models_v3.zip`.",
            tuffold_engine::MODELS_URL
        );
    }
    Ok(())
}

fn go(project: &Path, files: &tuffold_engine::GoFiles) -> Result<(), String> {
    let mut project = Project::open(project).map_err(|e| e.to_string())?;
    let report = tuffold_engine::annotate_go(&mut project, files).map_err(|e| e.to_string())?;
    out!(
        "GO: {} of {} proteins with a UniProt accession reach {} slim terms",
        report.annotated,
        report.uniprot_proteins,
        report.terms
    );
    if report.unknown_terms > 0 {
        eprintln!(
            "note: {} annotations name terms the ontology doesn't have (obsolete or newer); \
             they were left out",
            report.unknown_terms
        );
    }
    Ok(())
}

/// Write `to` with `write`, removing it again if writing fails.
fn write_new<T, E: std::fmt::Display>(
    to: &Path,
    write: impl FnOnce(&mut std::io::BufWriter<std::fs::File>) -> Result<T, E>,
) -> Result<T, String> {
    let file = std::fs::File::create(to).map_err(|e| format!("{}: {e}", to.display()))?;
    let mut out = std::io::BufWriter::new(file);
    let written = write(&mut out)
        .map_err(|e| e.to_string())
        .and_then(|value| {
            std::io::Write::flush(&mut out).map_err(|e| e.to_string())?;
            Ok(value)
        });
    if written.is_err() {
        drop(out);
        // Best effort: don't leave half a file behind.
        let _ = std::fs::remove_file(to);
    }
    written
}

fn export_mzid(project: &Project, to: &Path) -> Result<(), String> {
    let report = write_new(to, |out| tuffold_engine::write_mzid(project, out))?;
    out!(
        "wrote {}: {} spectra ({} passing), {} protein groups ({} passing)",
        to.display(),
        report.spectra,
        report.passing_spectra,
        report.protein_groups,
        report.passing_groups
    );
    Ok(())
}

fn export_library(project: &Project, to: &Path) -> Result<(), String> {
    let report = write_new(to, |out| tuffold_engine::write_msp(project, out))?;
    out!("wrote {}: {} peptide ions", to.display(), report.entries);
    if report.without_peaks > 0 {
        eprintln!(
            "note: {} passing peptide ions have no spectrum with peaks; add the runs' spectra \
             files to include them",
            report.without_peaks
        );
    }
    if report.unnamed_modifications > 0 {
        eprintln!(
            "note: {} peptide ions carry a modification without a Unimod name and were left \
             out (.msp names every modification)",
            report.unnamed_modifications
        );
    }
    Ok(())
}

fn export(project: &Path, to: &Path) -> Result<(), String> {
    let extension = to
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase());
    if !matches!(extension.as_deref(), Some("xlsx" | "mzid" | "msp")) {
        return Err(format!(
            "{}: export to .xlsx (workbook), .mzid (mzIdentML), or .msp (spectral library)",
            to.display()
        ));
    }
    if to.exists() {
        return Err(format!(
            "{} already exists; choose a new file name",
            to.display()
        ));
    }
    let project = Project::open(project).map_err(|e| e.to_string())?;
    match extension.as_deref() {
        Some("mzid") => return export_mzid(&project, to),
        Some("msp") => return export_library(&project, to),
        _ => {}
    }
    let written = tuffold_engine::write_workbook(&project, to).map_err(|e| e.to_string())?;
    for (sheet, rows) in &written.sheets {
        out!("  {sheet}: {rows} rows");
    }
    for (sheet, total) in &written.cut {
        eprintln!(
            "warning: {sheet} has {total} rows but Excel holds {}; the sheet has the first \
             ones (the Summary sheet says so)",
            tuffold_engine::MAX_ROWS
        );
    }
    out!("wrote {}", to.display());
    Ok(())
}

fn mbr(path: &std::path::Path, files: &[PathBuf]) -> Result<(), String> {
    let mut project = tuffold_store::Project::open(path).map_err(|e| e.to_string())?;
    let reports =
        tuffold_engine::match_between_runs(&mut project, files).map_err(|e| e.to_string())?;
    for r in reports {
        let s = &r.summary;
        out!(
            "{}: {} of {} other runs aligned; {} precursors looked for, {} transferred at {}% \
             FDR (pairs won by targets {}, by decoys {}; {} known peaks{})",
            r.run,
            r.aligned,
            r.donors,
            s.candidates,
            s.accepted,
            tuffold_store::TRANSFER_FDR * 100.0,
            s.targets_won,
            s.decoys_won,
            s.known,
            if s.trained {
                ", model trained on them"
            } else {
                ", too few to train: ranked by retention-time error"
            }
        );
        if r.ms1_scans_damaged > 0 {
            out!(
                "  {} damaged MS1 scans in {}'s spectra file left out",
                r.ms1_scans_damaged,
                r.run
            );
        }
    }
    Ok(())
}
