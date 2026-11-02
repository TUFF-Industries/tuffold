//! MSFragger's own results inside a FragPipe folder.
//!
//! FragPipe's `psm.tsv` is Philosopher's report: PeptideProphet's or
//! Percolator's verdicts, already cut at FragPipe's FDR, decoys gone. Tuffold
//! cannot validate a list that has been filtered, so when the folder still
//! holds the search engine's own pepXML files (one per run, every rank-1
//! match, decoys included) those are imported instead, through the pepXML
//! importer.
//!
//! - A search file is a `.pepXML`/`.pep.xml` (optionally gzipped) in the
//!   folder or one of its experiment subfolders, except PeptideProphet's
//!   `interact-*` output and MSFragger's `_rankN` files past rank 1 (written
//!   when more than one hit per spectrum is reported for Percolator).
//! - Runs in an experiment subfolder belong to that experiment (the
//!   biosample, its category the name's prefix, as for `psm.tsv`); runs at
//!   the top level are their own biosample.
//! - A run found in two files, or files from different search engines, are
//!   refused rather than merged.
//! - Protein sequences come from the folder's FASTA (top level first, else the
//!   first experiment subfolder holding one), and only fill unknown ones.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use tuffold_core::ImportBatch;

use super::biosample;
use crate::fasta::{build_identifier_index, find_fasta, parse_fasta};
use crate::pepxml::import_pepxml;
use crate::xml::XmlImportError;

const SEARCH_EXTENSIONS: &[&str] = &[".pepxml", ".pep.xml", ".pepxml.gz", ".pep.xml.gz"];
const POST_PROCESSOR_PREFIX: &str = "interact-";

/// One search-engine pepXML in a FragPipe folder.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SearchFile {
    pub path: PathBuf,
    /// The experiment subfolder it is in, `None` at the top level.
    pub experiment: Option<String>,
}

/// The search engine's own pepXML files in `folder` and its experiment
/// subfolders, sorted by experiment, then file name. Hidden entries are
/// skipped.
pub fn find_search_files(folder: &Path) -> Result<Vec<SearchFile>, XmlImportError> {
    let mut found = Vec::new();
    for (path, is_dir) in sorted_entries(folder)? {
        if is_dir {
            let experiment = crate::tsv::file_name(&path);
            for (inner, inner_is_dir) in sorted_entries(&path)? {
                if !inner_is_dir && is_search_file(&crate::tsv::file_name(&inner)) {
                    found.push(SearchFile {
                        path: inner,
                        experiment: Some(experiment.clone()),
                    });
                }
            }
        } else if is_search_file(&crate::tsv::file_name(&path)) {
            found.push(SearchFile {
                path,
                experiment: None,
            });
        }
    }
    Ok(found)
}

/// Non-hidden entries of `dir`, sorted by name, with whether each is a
/// directory.
fn sorted_entries(dir: &Path) -> Result<Vec<(PathBuf, bool)>, XmlImportError> {
    let io = |source| XmlImportError::Io {
        file: dir.display().to_string(),
        source,
    };
    let mut entries = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(io)? {
        let entry = entry.map_err(io)?;
        if entry.file_name().to_string_lossy().starts_with('.') {
            continue;
        }
        let is_dir = entry.file_type().map_err(io)?.is_dir();
        entries.push((entry.path(), is_dir));
    }
    entries.sort();
    Ok(entries)
}

fn is_search_file(name: &str) -> bool {
    let lowered = name.to_ascii_lowercase();
    let Some(stem) = SEARCH_EXTENSIONS
        .iter()
        .filter(|extension| lowered.ends_with(*extension))
        .max_by_key(|extension| extension.len())
        .map(|extension| &lowered[..lowered.len() - extension.len()])
    else {
        return false;
    };
    if lowered.starts_with(POST_PROCESSOR_PREFIX) {
        return false;
    }
    // `run_rank2.pepXML` and beyond hold lower-ranked hits only.
    match stem.rsplit_once("_rank") {
        Some((_, rank)) if !rank.is_empty() && rank.bytes().all(|b| b.is_ascii_digit()) => {
            rank == "1"
        }
        _ => true,
    }
}

/// Import every search file into one batch.
pub fn import_search_files(
    folder: &Path,
    files: &[SearchFile],
) -> Result<ImportBatch, XmlImportError> {
    let mut merged = ImportBatch::default();
    let mut run_files: HashMap<String, PathBuf> = HashMap::new();
    let mut proteins: HashSet<String> = HashSet::new();
    for file in files {
        let mut batch = import_pepxml(&file.path)?;
        if merged.software.name.is_empty() {
            merged.software = batch.software.clone();
        } else if batch.software.name != merged.software.name {
            return Err(unsupported(
                &file.path,
                format!(
                    "is a {} search, but {} holds {} searches; import them separately",
                    batch.software.name,
                    folder.display(),
                    merged.software.name
                ),
            ));
        }
        for run in &batch.ms_runs {
            if let Some(first) = run_files.insert(run.name.clone(), file.path.clone()) {
                return Err(unsupported(
                    &file.path,
                    format!(
                        "run {} is also in {}; remove one copy",
                        run.name,
                        first.display()
                    ),
                ));
            }
        }
        if let Some(experiment) = &file.experiment {
            for run in &mut batch.ms_runs {
                run.biosample = experiment.clone();
            }
            batch.biosamples = vec![biosample(experiment)];
        }
        for sample in batch.biosamples {
            if !merged
                .biosamples
                .iter()
                .any(|known| known.name == sample.name)
            {
                merged.biosamples.push(sample);
            }
        }
        merged.ms_runs.extend(batch.ms_runs);
        merged.spectra.extend(batch.spectra);
        merged.psms.extend(batch.psms);
        merged.proteins.extend(
            batch
                .proteins
                .into_iter()
                .filter(|protein| proteins.insert(protein.accession.clone())),
        );
    }
    fill_sequences(folder, files, &mut merged)?;
    Ok(merged)
}

/// Sequences from the folder's FASTA for proteins that have none.
fn fill_sequences(
    folder: &Path,
    files: &[SearchFile],
    batch: &mut ImportBatch,
) -> Result<(), XmlImportError> {
    let experiments = files
        .iter()
        .filter_map(|file| file.experiment.as_ref())
        .map(|experiment| folder.join(experiment));
    let mut fasta = None;
    for dir in std::iter::once(folder.to_path_buf()).chain(experiments) {
        let io = |source| XmlImportError::Io {
            file: dir.display().to_string(),
            source,
        };
        if let Some(path) = find_fasta(&dir).map_err(io)? {
            fasta = Some(parse_fasta(&path).map_err(|source| XmlImportError::Io {
                file: path.display().to_string(),
                source,
            })?);
            break;
        }
    }
    let Some(records) = fasta else {
        return Ok(());
    };
    let index = build_identifier_index(&records);
    for protein in batch
        .proteins
        .iter_mut()
        .filter(|protein| protein.sequence.is_none())
    {
        if let Some(&at) = index.get(&protein.accession) {
            protein.sequence = Some(records[at].sequence.clone());
        }
    }
    Ok(())
}

fn unsupported(path: &Path, problem: String) -> XmlImportError {
    XmlImportError::Unsupported {
        file: path.display().to_string(),
        problem,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn search_files_are_the_engines_not_the_post_processors() {
        for (name, expected) in [
            ("A-01.pepXML", true),
            ("A-01.pep.xml.gz", true),
            ("A-01_rank1.pepXML", true),
            ("A-01_rank2.pepXML", false),
            ("A-01_rank12.pepXML", false),
            ("interact-A-01.pep.xml", false),
            ("Interact-A-01.pepXML", false),
            ("my_ranked_run.pepXML", true),
            ("A-01_rank.pepXML", true),
            ("A-01.pin", false),
            ("psm.tsv", false),
        ] {
            assert_eq!(is_search_file(name), expected, "{name}");
        }
    }
}
