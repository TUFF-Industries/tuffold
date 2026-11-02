//! Gene Ontology: the ontology (OBO), a GO slim (OBO), and annotations
//! (GAF 2.x, gzip-aware), rolled up to slim terms per UniProt accession.
//!
//! - Terms keep their name, aspect (`namespace`), and parents through
//!   `is_a` and `relationship: part_of`, the relations GO slims are mapped
//!   through. Obsolete terms are dropped (annotations to them don't count).
//! - GAF rows with a `NOT` qualifier are skipped; every other row maps its
//!   DB object ID (column 2, a UniProt accession in GOA files) to its GO ID
//!   (column 5).
//! - An accession gets every slim term its annotations reach, the
//!   annotated term itself included. Annotations reaching no slim term are
//!   counted, not shown.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
pub use tuffold_core::uniprot::uniprot_accession;

#[derive(Debug, thiserror::Error)]
pub enum GoError {
    #[error("{file}: {source}")]
    Io {
        file: String,
        #[source]
        source: io::Error,
    },
    #[error("{file}: no [Term] stanzas with GO ids; is it a GO .obo file?")]
    NotAnOntology { file: String },
    #[error("{file}: line {line}: a GAF row needs at least 9 tab-separated columns, found {found}")]
    BadGafRow {
        file: String,
        line: usize,
        found: usize,
    },
}

/// Which of GO's three ontologies a term belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Aspect {
    Process,
    Function,
    Component,
}

impl Aspect {
    fn from_namespace(namespace: &str) -> Option<Self> {
        match namespace {
            "biological_process" => Some(Self::Process),
            "molecular_function" => Some(Self::Function),
            "cellular_component" => Some(Self::Component),
            _ => None,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Process => "process",
            Self::Function => "function",
            Self::Component => "component",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GoTerm {
    pub id: String,
    pub name: String,
    pub aspect: Aspect,
    /// `is_a` and `part_of` parents.
    pub parents: Vec<String>,
}

/// GO terms by id.
#[derive(Debug, Clone, Default)]
pub struct Ontology {
    pub terms: HashMap<String, GoTerm>,
}

fn open(path: &Path) -> Result<Box<dyn BufRead>, GoError> {
    let io = |source| GoError::Io {
        file: path.display().to_string(),
        source,
    };
    let file = File::open(path).map_err(io)?;
    let gzip = path
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"));
    Ok(if gzip {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

/// Read an OBO file's GO terms.
pub fn read_obo(path: &Path) -> Result<Ontology, GoError> {
    let file = path.display().to_string();
    let ontology = parse_obo(open(path)?).map_err(|source| GoError::Io {
        file: file.clone(),
        source,
    })?;
    if ontology.terms.is_empty() {
        return Err(GoError::NotAnOntology { file });
    }
    Ok(ontology)
}

#[derive(Default)]
struct Stanza {
    id: Option<String>,
    name: Option<String>,
    aspect: Option<Aspect>,
    parents: Vec<String>,
    obsolete: bool,
}

impl Stanza {
    fn finish(self, terms: &mut HashMap<String, GoTerm>) {
        if let (Some(id), Some(name), Some(aspect), false) =
            (self.id, self.name, self.aspect, self.obsolete)
        {
            terms.insert(
                id.clone(),
                GoTerm {
                    id,
                    name,
                    aspect,
                    parents: self.parents,
                },
            );
        }
    }
}

/// Parse OBO text: `[Term]` stanzas with GO ids.
pub fn parse_obo(reader: impl BufRead) -> io::Result<Ontology> {
    let mut terms = HashMap::new();
    let mut stanza: Option<Stanza> = None;
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.starts_with('[') {
            if let Some(done) = stanza.take() {
                done.finish(&mut terms);
            }
            stanza = (line == "[Term]").then(Stanza::default);
            continue;
        }
        let Some(current) = stanza.as_mut() else {
            continue;
        };
        let Some((tag, value)) = line.split_once(':') else {
            continue;
        };
        // Drop trailing comments (`! name`) and qualifiers (`{...}`).
        let value = value.split(" ! ").next().unwrap_or("").trim();
        let first = value.split_whitespace().next().unwrap_or("");
        match tag.trim() {
            "id" if first.starts_with("GO:") => current.id = Some(first.to_string()),
            "name" => current.name = Some(value.to_string()),
            "namespace" => current.aspect = Aspect::from_namespace(first),
            "is_a" if first.starts_with("GO:") => current.parents.push(first.to_string()),
            "relationship" => {
                let mut parts = value.split_whitespace();
                if parts.next() == Some("part_of")
                    && let Some(parent) = parts.next().filter(|p| p.starts_with("GO:"))
                {
                    current.parents.push(parent.to_string());
                }
            }
            "is_obsolete" => current.obsolete = first == "true",
            _ => {}
        }
    }
    if let Some(done) = stanza {
        done.finish(&mut terms);
    }
    Ok(Ontology { terms })
}

/// Accession → GO ids from a GAF file (`NOT` rows skipped).
pub fn read_gaf(path: &Path) -> Result<Vec<(String, String)>, GoError> {
    let file = path.display().to_string();
    let mut annotations = Vec::new();
    for (index, line) in open(path)?.lines().enumerate() {
        let line = line.map_err(|source| GoError::Io {
            file: file.clone(),
            source,
        })?;
        if line.starts_with('!') || line.trim().is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 9 {
            return Err(GoError::BadGafRow {
                file,
                line: index + 1,
                found: columns.len(),
            });
        }
        if columns[3].split('|').any(|q| q.eq_ignore_ascii_case("NOT")) {
            continue;
        }
        annotations.push((columns[1].to_string(), columns[4].to_string()));
    }
    Ok(annotations)
}

/// What the roll-up found.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SlimMap {
    /// Accession → slim term ids its annotations reach.
    pub by_accession: HashMap<String, BTreeSet<String>>,
    /// Annotations to terms the ontology doesn't have (obsolete or newer).
    pub unknown_terms: usize,
    /// Annotations reaching no slim term.
    pub outside_slim: usize,
}

/// Roll annotations up to the slim terms each reaches.
pub fn slim_map(
    ontology: &Ontology,
    slim: &HashSet<String>,
    annotations: &[(String, String)],
) -> SlimMap {
    let mut reach: HashMap<&str, BTreeSet<String>> = HashMap::new();
    let mut map = SlimMap::default();
    for (accession, term) in annotations {
        if !ontology.terms.contains_key(term) {
            map.unknown_terms += 1;
            continue;
        }
        if !reach.contains_key(term.as_str()) {
            let found = slim_terms(ontology, slim, term);
            reach.insert(term.as_str(), found);
        }
        let found = &reach[term.as_str()];
        if found.is_empty() {
            map.outside_slim += 1;
            continue;
        }
        map.by_accession
            .entry(accession.clone())
            .or_default()
            .extend(found.iter().cloned());
    }
    map
}

/// The slim terms among `term` and its ancestors.
fn slim_terms(ontology: &Ontology, slim: &HashSet<String>, term: &str) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    let mut seen: HashSet<&str> = HashSet::new();
    let mut stack = vec![term];
    while let Some(current) = stack.pop() {
        if !seen.insert(current) {
            continue;
        }
        if slim.contains(current) {
            found.insert(current.to_string());
        }
        if let Some(node) = ontology.terms.get(current) {
            stack.extend(node.parents.iter().map(String::as_str));
        }
    }
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    const OBO: &str = "format-version: 1.2
ontology: go

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0007049
name: cell cycle
namespace: biological_process
is_a: GO:0009987 ! cellular process

[Term]
id: GO:0009987
name: cellular process
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000082
name: G1/S transition of mitotic cell cycle
namespace: biological_process
relationship: part_of GO:0007049 ! cell cycle

[Term]
id: GO:0005737
name: cytoplasm
namespace: cellular_component

[Term]
id: GO:0000001
name: old term
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
name: part of
";

    #[test]
    fn obo_terms_keep_names_aspects_and_both_relations() {
        let ontology = parse_obo(OBO.as_bytes()).unwrap();
        // Five live terms; the obsolete one and the Typedef are dropped.
        assert_eq!(ontology.terms.len(), 5);
        let transition = &ontology.terms["GO:0000082"];
        assert_eq!(transition.name, "G1/S transition of mitotic cell cycle");
        assert_eq!(transition.aspect, Aspect::Process);
        assert_eq!(transition.parents, ["GO:0007049"]);
        assert_eq!(ontology.terms["GO:0005737"].aspect, Aspect::Component);
    }

    #[test]
    fn annotations_roll_up_to_the_slim_terms_they_reach() {
        let ontology = parse_obo(OBO.as_bytes()).unwrap();
        let slim: HashSet<String> = ["GO:0007049", "GO:0005737"].map(String::from).into();
        let annotations: Vec<(String, String)> = [
            // part_of then is_a reaches cell cycle.
            ("P12345", "GO:0000082"),
            ("P12345", "GO:0005737"),
            // cellular process is above the slim term: no slim term.
            ("Q99999", "GO:0009987"),
            // Obsolete (dropped from the ontology): unknown.
            ("Q99999", "GO:0000001"),
        ]
        .map(|(a, t)| (a.to_string(), t.to_string()))
        .to_vec();
        let map = slim_map(&ontology, &slim, &annotations);
        assert_eq!(
            map.by_accession["P12345"],
            BTreeSet::from(["GO:0005737".to_string(), "GO:0007049".to_string()])
        );
        assert!(!map.by_accession.contains_key("Q99999"));
        assert_eq!((map.outside_slim, map.unknown_terms), (1, 1));
    }

    #[test]
    fn gaf_rows_map_accessions_and_skip_not() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.gaf");
        std::fs::write(
            &path,
            "!gaf-version: 2.2\n\
             UniProtKB\tP12345\tCDK1\tinvolved_in\tGO:0000082\tPMID:1\tIDA\t\tP\n\
             UniProtKB\tP12345\tCDK1\tNOT|located_in\tGO:0005737\tPMID:1\tIDA\t\tC\n",
        )
        .unwrap();
        assert_eq!(
            read_gaf(&path).unwrap(),
            [("P12345".to_string(), "GO:0000082".to_string())]
        );
        std::fs::write(&path, "UniProtKB\tP1\tX\n").unwrap();
        assert!(matches!(
            read_gaf(&path),
            Err(GoError::BadGafRow {
                line: 1,
                found: 3,
                ..
            })
        ));
    }
}
