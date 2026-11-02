//! `tuffold search`: the built-in search (Sage's search core).

use std::path::{Path, PathBuf};

use tuffold_engine::{SearchOptions, SearchSettings, prepare_search, search_prepared};
use tuffold_io::is_spectra_file;
use tuffold_store::Project;

#[derive(clap::Args)]
pub struct SearchArgs {
    /// Precursor tolerance, ± ppm [default: 20].
    #[arg(long)]
    precursor_ppm: Option<f32>,
    /// Fragment tolerance, ± ppm [default: 20].
    #[arg(long)]
    fragment_ppm: Option<f32>,
    /// Fragment tolerance in ± Da instead (low-resolution MS/MS, e.g. 0.5
    /// for an ion trap).
    #[arg(long, conflicts_with = "fragment_ppm")]
    fragment_da: Option<f32>,
    /// trypsin, trypsin/p, lys-c, lys-n, arg-c, glu-c, asp-n, chymotrypsin
    /// [default: trypsin].
    #[arg(long)]
    enzyme: Option<String>,
    /// Allow one non-enzymatic end.
    #[arg(long)]
    semi: bool,
    /// Missed cleavages allowed [default: 2].
    #[arg(long)]
    missed_cleavages: Option<u8>,
    /// A fixed modification, SITE:MOD (e.g. K:TMT6plex, C:57.021464; sites
    /// ^ $ are peptide termini, [ ] protein termini). Repeatable.
    #[arg(long = "fixed", value_name = "SITE:MOD")]
    fixed: Vec<String>,
    /// A variable modification, SITE:MOD (e.g. [:Acetyl, S:Phospho).
    /// Repeatable.
    #[arg(long = "variable", value_name = "SITE:MOD")]
    variable: Vec<String>,
    /// Isobaric label fixed on K and peptide N-termini: tmt, tmtpro, itraq4,
    /// or itraq8.
    #[arg(long)]
    label: Option<String>,
    /// Variable modifications per peptide [default: 2].
    #[arg(long)]
    max_variable_mods: Option<usize>,
    /// Do not fix carbamidomethyl on cysteine.
    #[arg(long)]
    no_carbamidomethyl: bool,
    /// Do not search oxidized methionine.
    #[arg(long)]
    no_oxidation: bool,
    /// Most memory the search index may take, in GB [default: three quarters
    /// of the memory free when the search starts]. A search that would need
    /// more is refused before it starts.
    #[arg(long, value_name = "GB")]
    max_memory: Option<f64>,
}

impl SearchArgs {
    fn settings(&self) -> Result<SearchSettings, String> {
        SearchSettings::from_options(&SearchOptions {
            precursor_ppm: self.precursor_ppm,
            fragment_ppm: self.fragment_ppm,
            fragment_da: self.fragment_da,
            enzyme: self.enzyme.clone(),
            semi_specific: self.semi,
            missed_cleavages: self.missed_cleavages,
            fixed: self.fixed.clone(),
            variable: self.variable.clone(),
            label: self.label.clone(),
            max_variable_mods: self.max_variable_mods,
            no_carbamidomethyl: self.no_carbamidomethyl,
            no_oxidation: self.no_oxidation,
            max_memory_gb: self.max_memory,
        })
    }
}

pub fn search(
    spectra: &[PathBuf],
    fasta: &Path,
    into: &Path,
    args: &SearchArgs,
) -> Result<(), String> {
    // Check the inputs before creating or touching the project.
    if let Some(other) = spectra.iter().find(|input| !is_spectra_file(input)) {
        if tuffold_io::is_thermo_raw(other) {
            return Err(format!("{}: {}", other.display(), tuffold_io::THERMO_RAW));
        }
        return Err(format!(
            "{}: not a spectra file; the search reads .mzML and .mgf files and Bruker timsTOF .d folders",
            other.display()
        ));
    }
    let settings = args.settings()?;
    if !fasta.is_file() {
        return Err(format!("{}: no such FASTA file", fasta.display()));
    }
    out!("search: {}", settings.describe());
    // The index first: a search too big for memory is refused before a
    // project is created or touched.
    let searcher = prepare_search(fasta, settings).map_err(|error| error.to_string())?;
    let mut project = if into.exists() {
        Project::open(into)
    } else {
        Project::create(into)
    }
    .map_err(|error| error.to_string())?;
    let report =
        search_prepared(&mut project, &searcher, spectra).map_err(|error| error.to_string())?;
    out!(
        "digested {} proteins: {} peptide forms with generated decoys (index estimated at \
         {:.1} GB of a {:.1} GB memory budget)",
        report.targets,
        report.peptides,
        report.index_bytes as f64 / 1e9,
        report.memory_budget as f64 / 1e9
    );
    if report.dropped_decoys > 0 {
        out!(
            "note: {} decoy entries in the FASTA were left out; the search generates its own \
             (rev_ prefix)",
            report.dropped_decoys
        );
    }
    if !report.excluded_targets.is_empty() {
        out!(
            "warning: {} FASTA entries were left out (accession contains rev_ or repeats): {}",
            report.excluded_targets.len(),
            report.excluded_targets.join(", ")
        );
    }
    for run in &report.runs {
        let input = run.input.display();
        match &run.outcome {
            Ok(searched) => {
                let found = searched.imported.as_ref().map_or_else(
                    || "no matches, nothing imported".to_string(),
                    |done| {
                        format!(
                            "{} PSMs imported, peaks attached to {} spectra; {}",
                            done.import.psms,
                            done.attach.attached,
                            tuffold_engine::ms1_note(&done.attach)
                        )
                    },
                );
                out!(
                    "searched {input} (run {}): {} spectra; {found}",
                    run.run,
                    searched.searched
                );
                if searched.centroided > 0 {
                    out!(
                        "note: {input}: {} profile-mode spectra centroided before searching",
                        searched.centroided
                    );
                }
                if let Some(skipped) = searched.skipped.describe() {
                    out!("note: {input}: not searched: {skipped}");
                }
            }
            Err(message) => out!("failed {input}: {message}"),
        }
    }
    out!("next: tuffold validate {}", into.display());
    match report.failures() {
        0 => Ok(()),
        failed => Err(format!("{failed} spectra file(s) could not be searched")),
    }
}
