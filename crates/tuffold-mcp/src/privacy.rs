//! Privacy mode: the assistant sees aliases, never the project's names.
//!
//! Started with `--private` (or `TUFFOLD_MCP_PRIVATE=names`), every reply
//! has its biosample, category, run, spectra-file, and import names replaced
//! by aliases (`sample-1`, `group-2`, `run-3`, `file-1`, `import-2`); with
//! `--private=proteins` protein accessions become `protein-N` and their
//! descriptions are left out too. Names the assistant sends back (a
//! category to compare, a sample to group) are mapped to the real ones
//! here, so the tools work unchanged and nothing real leaves the computer
//! except what the assistant was itself given (the paths it is asked to
//! open). The level is set by whoever starts the server; no tool can change
//! it.
//!
//! Aliases are numbered per project, in the order the project stores the
//! names, the first time this server sees them, and kept for the server's
//! lifetime, so `sample-2` means the same sample in every reply. Every
//! string in a reply is rewritten, keys included: a whole value that is a
//! name becomes its alias, and names of [`MIN_TEXT_NAME`] characters or more
//! are replaced inside text too (messages, notes, paths), at word
//! boundaries, longest first. Shorter names only as whole values: a
//! category "A" would otherwise hit every sentence.
//!
//! Not hidden: peptide sequences (which can identify a protein, so
//! `proteins` hides names, not identity), numbers, software, and settings.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, PoisonError};

use serde_json::Value;
use tuffold_store::{Project, ProjectNames};

/// Names shorter than this are aliased only as whole values.
const MIN_TEXT_NAME: usize = 3;
/// The environment variable that turns privacy mode on.
pub const ENV: &str = "TUFFOLD_MCP_PRIVATE";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Level {
    Off,
    /// Sample, category, run, file, and import names.
    Names,
    /// Names, plus protein accessions and descriptions.
    Proteins,
}

impl Level {
    /// From the command line (`--private`, `--private=names`,
    /// `--private=proteins`, `--private=off`), else the [`ENV`] variable.
    pub fn parse(args: &[String], env: Option<&str>) -> Result<Self, String> {
        let mut level = match env.map(str::trim) {
            None | Some("") => Self::Off,
            Some(value) => Self::named(value)
                .ok_or_else(|| format!("{ENV}={value:?}: use off, names, or proteins"))?,
        };
        for arg in args {
            level = match arg.as_str() {
                "--private" => Self::Names,
                other => other
                    .strip_prefix("--private=")
                    .and_then(Self::named)
                    .ok_or_else(|| {
                        format!(
                            "unknown argument {other:?}; tuffold-mcp takes \
                             --private[=names|proteins|off]"
                        )
                    })?,
            };
        }
        Ok(level)
    }

    fn named(value: &str) -> Option<Self> {
        match value.to_ascii_lowercase().as_str() {
            "off" | "0" | "false" => Some(Self::Off),
            "names" | "on" | "1" | "true" => Some(Self::Names),
            "proteins" => Some(Self::Proteins),
            _ => None,
        }
    }

    /// What the assistant is told, appended to the server's instructions.
    pub fn instructions(self) -> Option<&'static str> {
        match self {
            Self::Off => None,
            Self::Names => Some(
                "Privacy mode: sample, category, run, and file names are aliases (sample-1, \
                 group-2, run-3, file-1, import-2). Use them as given; the server maps them \
                 back. Don't guess real names.",
            ),
            Self::Proteins => Some(
                "Privacy mode: sample, category, run, and file names are aliases (sample-1, \
                 group-2, run-3, file-1, import-2), and so are protein accessions (protein-7; \
                 descriptions are withheld). Use them as given; the server maps them back. \
                 Don't guess real names or identify proteins from their peptides.",
            ),
        }
    }
}

/// One project's aliases.
#[derive(Default)]
struct Book {
    /// Real name → alias.
    alias_of: HashMap<String, String>,
    /// Alias → real name.
    real_of: HashMap<String, String>,
    /// Aliases handed out per kind.
    issued: HashMap<&'static str, usize>,
    /// Names replaced inside text, by first character, longest first.
    in_text: HashMap<char, Vec<String>>,
}

impl Book {
    /// Give every name not seen before its alias.
    fn learn(&mut self, names: &ProjectNames, level: Level) {
        let mut added = false;
        for source in &names.sources {
            added |= self.learn_one("import", source);
            // The file name alone shows up in messages too.
            let base = source.rsplit(['/', '\\']).next().unwrap_or(source);
            if let Some(alias) = self.alias_of.get(source).cloned() {
                added |= self.learn_as(base, alias);
            }
        }
        for file in &names.spectra_files {
            added |= self.learn_one("file", file);
        }
        for sample in &names.biosamples {
            added |= self.learn_one("sample", sample);
        }
        for category in &names.categories {
            added |= self.learn_one("group", category);
        }
        for plex in &names.plexes {
            added |= self.learn_one("plex", plex);
        }
        for run in &names.runs {
            added |= self.learn_one("run", run);
        }
        if level == Level::Proteins {
            for accession in &names.accessions {
                added |= self.learn_one("protein", accession);
            }
        }
        if added {
            for list in self.in_text.values_mut() {
                list.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));
            }
        }
    }

    fn learn_one(&mut self, kind: &'static str, real: &str) -> bool {
        if real.is_empty() || self.alias_of.contains_key(real) {
            return false;
        }
        let issued = self.issued.entry(kind).or_default();
        *issued += 1;
        let alias = format!("{kind}-{issued}");
        self.real_of.insert(alias.clone(), real.to_string());
        self.learn_as(real, alias)
    }

    /// Map `real` to an existing alias (no new number).
    fn learn_as(&mut self, real: &str, alias: String) -> bool {
        if real.is_empty() || self.alias_of.contains_key(real) {
            return false;
        }
        self.alias_of.insert(real.to_string(), alias);
        if real.chars().count() >= MIN_TEXT_NAME
            && let Some(first) = real.chars().next()
        {
            self.in_text
                .entry(first)
                .or_default()
                .push(real.to_string());
        }
        true
    }

    fn reveal(&self, name: &str) -> String {
        self.real_of
            .get(name)
            .cloned()
            .unwrap_or_else(|| name.to_string())
    }

    /// `text` with every name replaced, in one pass (an alias is never
    /// rewritten again).
    fn hide_text(&self, text: &str) -> String {
        if let Some(alias) = self.alias_of.get(text) {
            return alias.clone();
        }
        let mut hidden = String::with_capacity(text.len());
        let mut at = 0;
        let mut previous: Option<char> = None;
        while let Some(current) = text[at..].chars().next() {
            let rest = &text[at..];
            let starts_word = previous.is_none_or(|c| !c.is_alphanumeric());
            let name = starts_word
                .then(|| self.in_text.get(&current))
                .flatten()
                .and_then(|names| {
                    names.iter().find(|name| {
                        rest.starts_with(name.as_str())
                            && rest[name.len()..]
                                .chars()
                                .next()
                                .is_none_or(|c| !c.is_alphanumeric())
                    })
                });
            match name {
                Some(name) => {
                    hidden.push_str(&self.alias_of[name]);
                    at += name.len();
                    previous = name.chars().last();
                }
                None => {
                    hidden.push(current);
                    at += current.len_utf8();
                    previous = Some(current);
                }
            }
        }
        hidden
    }

    fn hide(&self, value: Value, level: Level) -> Value {
        match value {
            Value::String(text) => Value::String(self.hide_text(&text)),
            Value::Array(items) => items
                .into_iter()
                .map(|item| self.hide(item, level))
                .collect(),
            Value::Object(fields) => fields
                .into_iter()
                .map(|(key, value)| {
                    let value = if level == Level::Proteins && key == "description" {
                        Value::Null
                    } else {
                        self.hide(value, level)
                    };
                    (self.hide_text(&key), value)
                })
                .collect(),
            other => other,
        }
    }
}

/// The server's aliases, one book per project.
pub struct Privacy {
    level: Level,
    books: Mutex<HashMap<PathBuf, Book>>,
}

impl Privacy {
    pub fn new(level: Level) -> Self {
        Self {
            level,
            books: Mutex::default(),
        }
    }

    pub fn level(&self) -> Level {
        self.level
    }

    /// Run `work` on the project's book after learning its current names.
    /// A project that can't be read (not created yet, say) adds none.
    fn with_book<T>(&self, project: &Path, work: impl FnOnce(&Book) -> T) -> T {
        let names = project
            .is_absolute()
            .then(|| {
                Project::open(project)
                    .and_then(|opened| opened.names())
                    .ok()
            })
            .flatten();
        let mut books = self.books.lock().unwrap_or_else(PoisonError::into_inner);
        let book = books.entry(project.to_path_buf()).or_default();
        if let Some(names) = names {
            book.learn(&names, self.level);
        }
        work(book)
    }

    /// The real name behind an alias the assistant sent; anything else is
    /// passed on as given (a new category it names, say).
    pub fn reveal(&self, project: &Path, name: &str) -> String {
        if self.level == Level::Off {
            return name.to_string();
        }
        self.with_book(project, |book| book.reveal(name))
    }

    /// A reply with every name replaced.
    pub fn hide(&self, project: &Path, value: Value) -> Value {
        if self.level == Level::Off {
            return value;
        }
        self.with_book(project, |book| book.hide(value, self.level))
    }

    /// An error message with every name replaced.
    pub fn hide_text(&self, project: &Path, text: &str) -> String {
        if self.level == Level::Off {
            return text.to_string();
        }
        self.with_book(project, |book| book.hide_text(text))
    }
}

#[cfg(test)]
mod tests {
    use serde_json::json;

    use super::*;

    fn names() -> ProjectNames {
        ProjectNames {
            biosamples: vec!["Patient 7 tumor".into(), "Patient 7".into(), "B1".into()],
            categories: vec!["cancer".into(), "A".into()],
            plexes: vec!["cohort_tmt".into()],
            runs: vec!["Patient 7 tumor".into(), "run_B1".into()],
            spectra_files: vec!["/data/lab/run_B1.mzML".into()],
            sources: vec!["/data/lab/secret/cohort.pep.xml".into()],
            accessions: vec!["sp|P01116|RASK_HUMAN".into()],
        }
    }

    fn book(level: Level) -> Book {
        let mut book = Book::default();
        book.learn(&names(), level);
        book
    }

    #[test]
    fn aliases_are_numbered_per_kind_in_stored_order() {
        let book = book(Level::Names);
        // Sources first (and their file names), then files, samples,
        // categories, plexes, runs; a run named like its sample keeps the
        // sample's alias.
        for (real, alias) in [
            ("/data/lab/secret/cohort.pep.xml", "import-1"),
            ("cohort.pep.xml", "import-1"),
            ("/data/lab/run_B1.mzML", "file-1"),
            ("Patient 7 tumor", "sample-1"),
            ("Patient 7", "sample-2"),
            ("B1", "sample-3"),
            ("cancer", "group-1"),
            ("A", "group-2"),
            ("cohort_tmt", "plex-1"),
            ("run_B1", "run-1"),
        ] {
            assert_eq!(book.alias_of[real], alias, "{real}");
            assert_eq!(
                book.reveal(alias),
                if real == "cohort.pep.xml" {
                    "/data/lab/secret/cohort.pep.xml"
                } else {
                    real
                }
            );
        }
        assert!(!book.alias_of.contains_key("sp|P01116|RASK_HUMAN"));
        assert_eq!(book.reveal("treated"), "treated");
    }

    #[test]
    fn names_inside_text_are_replaced_at_word_boundaries_longest_first() {
        let book = book(Level::Names);
        assert_eq!(
            book.hide_text("Patient 7 tumor vs Patient 7 (cancer); Patient 70 stays"),
            "sample-1 vs sample-2 (group-1); Patient 70 stays"
        );
        // Short names only as whole values.
        assert_eq!(book.hide_text("A"), "group-2");
        assert_eq!(book.hide_text("A merged result, B1"), "A merged result, B1");
        assert_eq!(
            book.hide_text("could not read /data/lab/secret/cohort.pep.xml: line 3"),
            "could not read import-1: line 3"
        );
        assert_eq!(book.hide_text("cancerous"), "cancerous");
        // An alias that looks like a name is not rewritten again.
        assert_eq!(book.hide_text("run_B1 and run_B1x"), "run-1 and run_B1x");
    }

    #[test]
    fn replies_are_rewritten_keys_included() {
        let reply = json!({
            "biosamples": [{"name": "Patient 7", "category": "cancer"}, {"name": "B1", "category": "A"}],
            "spectra_by_biosample": {"Patient 7 tumor": 12},
            "groups": [{"accession": "sp|P01116|RASK_HUMAN", "description": "GTPase KRas", "q_value": 0.001}],
        });
        let names_only = book(Level::Names).hide(reply.clone(), Level::Names);
        assert_eq!(
            names_only,
            json!({
                "biosamples": [{"name": "sample-2", "category": "group-1"}, {"name": "sample-3", "category": "group-2"}],
                "spectra_by_biosample": {"sample-1": 12},
                "groups": [{"accession": "sp|P01116|RASK_HUMAN", "description": "GTPase KRas", "q_value": 0.001}],
            })
        );
        let proteins = book(Level::Proteins).hide(reply, Level::Proteins);
        assert_eq!(
            proteins["groups"],
            json!([{"accession": "protein-1", "description": null, "q_value": 0.001}])
        );
    }

    #[test]
    fn the_level_comes_from_whoever_starts_the_server() {
        let args = |list: &[&str]| list.iter().map(|a| a.to_string()).collect::<Vec<_>>();
        assert_eq!(Level::parse(&[], None), Ok(Level::Off));
        assert_eq!(Level::parse(&args(&["--private"]), None), Ok(Level::Names));
        assert_eq!(
            Level::parse(&args(&["--private=proteins"]), None),
            Ok(Level::Proteins)
        );
        assert_eq!(Level::parse(&[], Some("names")), Ok(Level::Names));
        // The command line wins over the environment.
        assert_eq!(
            Level::parse(&args(&["--private=off"]), Some("proteins")),
            Ok(Level::Off)
        );
        assert!(Level::parse(&args(&["--privat"]), None).is_err());
        assert!(Level::parse(&[], Some("maybe")).is_err());
    }
}
