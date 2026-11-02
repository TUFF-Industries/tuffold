//! Built-in search: digest a FASTA once, then for each spectra file search
//! its MS/MS spectra with Sage, import the PSMs (one import per file,
//! software "Sage"), then attach the peaks and trace MS1 intensity by
//! reading the file once more, as an import of it would. One file's MS/MS
//! spectra are held in memory at a time. Validation is a separate step, as
//! after any import.

use std::path::{Path, PathBuf};

use tuffold_core::{MsRun, PeakSpectrum, SourceFormat};
use tuffold_io::fasta::parse_fasta;
use tuffold_io::{read_spectra, spectra_run_name};
use tuffold_search::{SearchError, SearchSettings, Searcher, Skipped};

use crate::import::attach_file;
use crate::undo::{Checkpoint, UndoError, import_count};
use tuffold_store::{AttachReport, ImportReport, ImportSource, Project};

#[derive(Debug, thiserror::Error)]
pub enum SearchFailure {
    #[error("{path}: {source}")]
    Fasta {
        path: String,
        source: std::io::Error,
    },
    #[error(transparent)]
    Search(#[from] SearchError),
    #[error("nothing was searched: {0}")]
    Undo(#[from] UndoError),
}

/// The digested database, and what happened to each spectra file.
#[derive(Debug)]
pub struct SearchReport {
    pub targets: usize,
    /// Target and decoy peptide forms searched.
    pub peptides: usize,
    /// Peak bytes building the index was estimated to take, and the budget
    /// it was checked against.
    pub index_bytes: u64,
    pub memory_budget: u64,
    /// Decoys in the FASTA, replaced by generated ones.
    pub dropped_decoys: usize,
    /// Target accessions left out (they contain the decoy tag, or repeat).
    pub excluded_targets: Vec<String>,
    /// Set when results were imported but the undo copy couldn't be kept.
    pub undo_problem: Option<String>,
    pub runs: Vec<RunReport>,
}

impl SearchReport {
    pub fn failures(&self) -> usize {
        self.runs.iter().filter(|r| r.outcome.is_err()).count()
    }
}

#[derive(Debug)]
pub struct RunReport {
    pub input: PathBuf,
    pub run: String,
    pub outcome: Result<RunSearched, String>,
}

#[derive(Debug)]
pub struct RunSearched {
    pub searched: usize,
    /// Profile spectra centroided before searching.
    pub centroided: usize,
    pub skipped: Skipped,
    /// `None` when no spectrum matched: nothing was imported, so the file
    /// can be searched again (another FASTA, other settings).
    pub imported: Option<Imported>,
}

#[derive(Debug)]
pub struct Imported {
    pub import: ImportReport,
    pub attach: AttachReport,
}

/// Search `spectra` files against `fasta` and import the results. A FASTA
/// or settings problem, or an index too big for the memory budget, stops
/// everything before the project is touched; a spectra file that cannot be
/// read or imported is reported and the others still run.
pub fn search_files(
    project: &mut Project,
    fasta: &Path,
    spectra: &[PathBuf],
    settings: SearchSettings,
) -> Result<SearchReport, SearchFailure> {
    let searcher = prepare_search(fasta, settings)?;
    search_prepared(project, &searcher, spectra)
}

/// Digest `fasta` into the search index, refusing an index that would not
/// fit the memory budget. Needs no project, so a front end that would create
/// one can do this first.
pub fn prepare_search(fasta: &Path, settings: SearchSettings) -> Result<Searcher, SearchFailure> {
    let records = parse_fasta(fasta).map_err(|source| SearchFailure::Fasta {
        path: fasta.display().to_string(),
        source,
    })?;
    Ok(Searcher::new(&records, settings)?)
}

/// [`search_files`] with the index already built by [`prepare_search`].
pub fn search_prepared(
    project: &mut Project,
    searcher: &Searcher,
    spectra: &[PathBuf],
) -> Result<SearchReport, SearchFailure> {
    let label = match spectra {
        [one] => format!(
            "search {}",
            one.file_name().map_or_else(
                || one.display().to_string(),
                |n| n.to_string_lossy().into_owned()
            )
        ),
        many => format!("search {} spectra files", many.len()),
    };
    let checkpoint = Checkpoint::take(project, &label)?;
    let before = import_count(project);
    let runs: Vec<RunReport> = spectra
        .iter()
        .map(|input| {
            let absolute = std::fs::canonicalize(input).unwrap_or_else(|_| input.clone());
            let run = spectra_run_name(&absolute);
            let outcome = read_all(&absolute)
                .and_then(|peaks| search_run(project, searcher, &absolute, &run, &peaks));
            RunReport {
                input: absolute,
                run,
                outcome,
            }
        })
        .collect();
    // Counted from the project, so a run whose PSMs imported but whose peaks
    // didn't attach is still undoable.
    let imported = import_count(project) > before;
    let undo_problem = match (imported, checkpoint) {
        (true, Some(checkpoint)) => checkpoint.keep().err().map(|error| {
            format!("the undo copy could not be kept, so this can't be undone ({error})")
        }),
        _ => None,
    };
    let (index, memory_budget) = searcher.index_estimate();
    Ok(SearchReport {
        targets: searcher.targets(),
        peptides: searcher.peptides(),
        index_bytes: index.peak_bytes,
        memory_budget,
        dropped_decoys: searcher.dropped_decoys(),
        excluded_targets: searcher.excluded_targets().to_vec(),
        undo_problem,
        runs,
    })
}

fn read_all(path: &Path) -> Result<Vec<PeakSpectrum>, String> {
    read_spectra(path)
        .and_then(Iterator::collect)
        .map_err(|error| error.to_string())
}

/// Search one run's spectra, import the PSMs, and attach the file. The run
/// is its own biosample, as with any spectra file imported on its own.
fn search_run(
    project: &mut Project,
    searcher: &Searcher,
    source: &Path,
    run: &str,
    peaks: &[PeakSpectrum],
) -> Result<RunSearched, String> {
    let ms_run = MsRun {
        name: run.to_string(),
        biosample: run.to_string(),
        spectra_file: source
            .file_name()
            .map(|name| name.to_string_lossy().into_owned()),
    };
    // Sage asserts on inputs it cannot handle; one bad file must not take
    // the other files (or the app) down with it.
    let found = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        searcher.search_run(&ms_run, peaks)
    }))
    .map_err(|_| "the search engine failed on this file; nothing was imported".to_string())?;
    if found.batch.psms.is_empty() {
        return Ok(RunSearched {
            searched: found.searched,
            centroided: found.centroided,
            skipped: found.skipped,
            imported: None,
        });
    }
    let source_path = source;
    let source = ImportSource {
        path: source.display().to_string(),
        format: SourceFormat::SageSearch,
    };
    let import = project
        .import(&source, &found.batch)
        .map_err(|error| error.to_string())?;
    let attach = attach_file(project, source_path, run).map_err(|error| {
        format!(
            "{} PSMs were imported, but the peaks could not be attached: {error}",
            import.psms
        )
    })?;
    Ok(RunSearched {
        searched: found.searched,
        centroided: found.centroided,
        skipped: found.skipped,
        imported: Some(Imported { import, attach }),
    })
}
