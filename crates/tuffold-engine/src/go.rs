//! GO annotations for a project: read the ontology, a GO slim, and GAF
//! annotation files, roll the project's proteins' annotations up to the
//! slim's terms, and store them (replacing any earlier GO import). Only the
//! project's own UniProt accessions are kept, so a 15 MB GAF leaves a few
//! kilobytes behind.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use tuffold_io::go::{GoError, read_gaf, read_obo, slim_map};
use tuffold_store::{Project, StoreError};

#[derive(Debug, thiserror::Error)]
pub enum GoImportError {
    #[error(transparent)]
    Go(#[from] GoError),
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error("{0}: the GO slim has no terms the ontology knows; is it a GO slim .obo?")]
    EmptySlim(String),
}

/// The files a GO import reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GoFiles {
    /// The full ontology (`go-basic.obo`).
    pub ontology: PathBuf,
    /// A GO slim (`goslim_generic.obo`).
    pub slim: PathBuf,
    /// GAF annotation files (`goa_human.gaf.gz`), keyed by UniProt.
    pub annotations: Vec<PathBuf>,
}

impl GoFiles {
    /// Sort picked files into ontology, slim, and annotations: `.gaf` /
    /// `.gaf.gz` are annotations, an `.obo` named like a slim is the slim,
    /// another `.obo` the ontology. Says what is missing or doubled.
    pub fn classify(paths: &[PathBuf]) -> Result<Self, String> {
        let name = |path: &Path| {
            path.file_name()
                .map(|n| n.to_string_lossy().to_ascii_lowercase())
                .unwrap_or_default()
        };
        let (mut ontology, mut slim, mut annotations) = (Vec::new(), Vec::new(), Vec::new());
        for path in paths {
            let lower = name(path);
            if lower.ends_with(".gaf") || lower.ends_with(".gaf.gz") {
                annotations.push(path.clone());
            } else if lower.ends_with(".obo") && lower.contains("slim") {
                slim.push(path.clone());
            } else if lower.ends_with(".obo") {
                ontology.push(path.clone());
            } else {
                return Err(format!(
                    "{}: not a GO file (expected .obo or .gaf/.gaf.gz)",
                    path.display()
                ));
            }
        }
        let one = |found: Vec<PathBuf>, what: &str, example: &str| match found.as_slice() {
            [only] => Ok(only.clone()),
            [] => Err(format!("pick {what} too (e.g. {example})")),
            _ => Err(format!("pick one {what}, not {}", found.len())),
        };
        if annotations.is_empty() {
            return Err("pick annotation files too (e.g. goa_human.gaf.gz)".into());
        }
        Ok(Self {
            ontology: one(ontology, "GO ontology", "go-basic.obo")?,
            slim: one(slim, "GO slim", "goslim_generic.obo")?,
            annotations,
        })
    }
}

/// What a GO import found.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GoReport {
    /// The project's proteins with a UniProt accession.
    pub uniprot_proteins: usize,
    /// Of those, how many reach at least one slim term.
    pub annotated: usize,
    /// Slim terms some protein reaches.
    pub terms: usize,
    /// Annotations (of the project's proteins) to terms the ontology lacks.
    pub unknown_terms: usize,
}

/// Read the GO files and store the project's slim annotations.
pub fn annotate_go(project: &mut Project, files: &GoFiles) -> Result<GoReport, GoImportError> {
    let ontology = read_obo(&files.ontology)?;
    let slim: HashSet<String> = read_obo(&files.slim)?
        .terms
        .into_keys()
        .filter(|id| ontology.terms.contains_key(id))
        .collect();
    if slim.is_empty() {
        return Err(GoImportError::EmptySlim(files.slim.display().to_string()));
    }
    let wanted: HashSet<String> = project.uniprot_accessions()?.into_iter().collect();
    let mut annotations = Vec::new();
    for path in &files.annotations {
        annotations.extend(
            read_gaf(path)?
                .into_iter()
                .filter(|(accession, _)| wanted.contains(accession)),
        );
    }
    let map = slim_map(&ontology, &slim, &annotations);
    let used: HashSet<&String> = map.by_accession.values().flatten().collect();
    let mut terms: Vec<(String, String, String)> = used
        .iter()
        .map(|id| {
            let term = &ontology.terms[id.as_str()];
            (
                term.id.clone(),
                term.name.clone(),
                term.aspect.label().to_string(),
            )
        })
        .collect();
    terms.sort();
    let mut pairs: Vec<(String, String)> = map
        .by_accession
        .iter()
        .flat_map(|(accession, ids)| ids.iter().map(move |id| (accession.clone(), id.clone())))
        .collect();
    pairs.sort();
    project.set_go(&terms, &pairs)?;
    Ok(GoReport {
        uniprot_proteins: wanted.len(),
        annotated: map.by_accession.len(),
        terms: terms.len(),
        unknown_terms: map.unknown_terms,
    })
}
