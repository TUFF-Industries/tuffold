//! Shared plumbing for the streaming XML importers (pepXML, mzIdentML).
//!
//! Files are read as a stream (gzip-compressed when the name ends in `.gz`),
//! so memory stays flat however large the search result is. Errors carry the
//! line number, recovered by rescanning the file only when an error occurs.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use quick_xml::events::BytesStart;
use quick_xml::{Reader, XmlVersion};

#[derive(Debug, thiserror::Error)]
pub enum XmlImportError {
    #[error("{file}: {source}")]
    Io {
        file: String,
        #[source]
        source: io::Error,
    },
    #[error("{file}: line {line}: malformed XML: {message}")]
    Malformed {
        file: String,
        line: usize,
        message: String,
    },
    #[error("{file}: line {line}: <{element}> {problem}")]
    Invalid {
        file: String,
        line: usize,
        element: String,
        problem: String,
    },
    #[error("{file}: {problem}")]
    Unsupported { file: String, problem: String },
}

pub(crate) fn open(path: &Path) -> Result<Box<dyn BufRead>, XmlImportError> {
    let file = File::open(path).map_err(|source| io_error(path, source))?;
    let is_gzip = path
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"));
    Ok(if is_gzip {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

pub(crate) fn reader(path: &Path) -> Result<Reader<Box<dyn BufRead>>, XmlImportError> {
    let mut reader = Reader::from_reader(open(path)?);
    reader.config_mut().trim_text(true);
    Ok(reader)
}

pub(crate) fn io_error(path: &Path, source: io::Error) -> XmlImportError {
    XmlImportError::Io {
        file: crate::tsv::file_name(path),
        source,
    }
}

/// 1-based line of the first non-whitespace byte at or after `offset` in the
/// (decompressed) document. Event positions sit just after the previous
/// event, before any trimmed whitespace, so the skip lands on the element.
pub(crate) fn line_at(path: &Path, offset: u64) -> usize {
    let Ok(stream) = open(path) else {
        return 0;
    };
    let offset = usize::try_from(offset).unwrap_or(usize::MAX);
    let mut line = 1;
    for (position, byte) in stream.bytes().enumerate() {
        let Ok(byte) = byte else { break };
        if position >= offset && !byte.is_ascii_whitespace() {
            break;
        }
        if byte == b'\n' {
            line += 1;
        }
    }
    line
}

pub(crate) fn malformed(path: &Path, offset: u64, error: impl std::fmt::Display) -> XmlImportError {
    XmlImportError::Malformed {
        file: crate::tsv::file_name(path),
        line: line_at(path, offset),
        message: error.to_string(),
    }
}

pub(crate) fn invalid(
    path: &Path,
    offset: u64,
    element: &str,
    problem: impl Into<String>,
) -> XmlImportError {
    XmlImportError::Invalid {
        file: crate::tsv::file_name(path),
        line: line_at(path, offset),
        element: element.to_string(),
        problem: problem.into(),
    }
}

/// An element's attributes by local name (namespace prefixes dropped),
/// with entities unescaped.
pub(crate) struct Attrs {
    pub element: String,
    values: HashMap<String, String>,
}

impl Attrs {
    pub fn read(start: &BytesStart<'_>) -> Result<Self, String> {
        let element = start.local_name().as_ref().to_string();
        let mut values = HashMap::new();
        for attribute in start.attributes() {
            let attribute = attribute.map_err(|error| error.to_string())?;
            let key = attribute.key.local_name().as_ref().to_string();
            let value = attribute
                .normalized_value(XmlVersion::Implicit1_0)
                .map_err(|error| error.to_string())?;
            values.insert(key, value.into_owned());
        }
        Ok(Self { element, values })
    }

    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str).map(str::trim)
    }

    /// A present, non-empty attribute or a description of what is missing.
    pub fn required(&self, name: &str) -> Result<&str, String> {
        self.get(name)
            .filter(|value| !value.is_empty())
            .ok_or_else(|| format!("is missing its required '{name}' attribute"))
    }

    pub fn parse<T: std::str::FromStr>(&self, name: &str) -> Result<Option<T>, String> {
        match self.get(name).filter(|value| !value.is_empty()) {
            None => Ok(None),
            Some(value) => value
                .parse()
                .map(Some)
                .map_err(|_| format!("has '{name}' = {value:?}, which is not a valid number")),
        }
    }

    /// A finite floating-point attribute; NaN and infinities are rejected.
    pub fn number(&self, name: &str) -> Result<Option<f64>, String> {
        match self.parse::<f64>(name)? {
            Some(value) if !value.is_finite() => Err(format!("has non-finite '{name}' = {value}")),
            other => Ok(other),
        }
    }

    pub fn required_number(&self, name: &str) -> Result<f64, String> {
        self.number(name)?
            .ok_or_else(|| format!("is missing its required '{name}' attribute"))
    }
}

/// The file name of a path-like attribute value (`/a/b/run.mzML` -> `run.mzML`).
pub(crate) fn base_name(value: &str) -> &str {
    value.rsplit(['/', '\\']).next().unwrap_or(value)
}

const SPECTRA_EXTENSIONS: &[&str] = &[
    ".mzml", ".mzxml", ".mgf", ".raw", ".d", ".wiff", ".ms2", ".mzml.gz", ".mgf.gz",
];

/// The MS run named by a path-like value: its file name without a spectra
/// extension (`/data/run1.mzML` -> `run1`, `proteinscape://x/Fraction_X` ->
/// `Fraction_X`). Every importer names runs this way so results meet.
pub(crate) fn run_name(location: &str) -> String {
    let name = base_name(location);
    let lowered = name.to_ascii_lowercase();
    SPECTRA_EXTENSIONS
        .iter()
        .filter(|extension| lowered.ends_with(*extension))
        .max_by_key(|extension| extension.len())
        .map_or(name, |extension| &name[..name.len() - extension.len()])
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_names_drop_folders_schemes_and_spectra_extensions() {
        assert_eq!(run_name("/data/run1.mzML"), "run1");
        assert_eq!(run_name("C:\\raw\\HeLa_01.RAW"), "HeLa_01");
        assert_eq!(run_name("proteinscape://server/Fraction_X"), "Fraction_X");
        assert_eq!(run_name("/data/runA"), "runA");
        assert_eq!(run_name("x.mzML.gz"), "x");
    }

    #[test]
    fn base_names_strip_both_separators() {
        assert_eq!(base_name("/data/run1"), "run1");
        assert_eq!(base_name("C:\\data\\run1.mzML"), "run1.mzML");
        assert_eq!(base_name("run1"), "run1");
    }

    #[test]
    fn line_at_counts_newlines_before_the_offset() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.xml");
        std::fs::write(&path, "a\nb\nc").unwrap();
        assert_eq!(line_at(&path, 0), 1);
        assert_eq!(line_at(&path, 2), 2);
        assert_eq!(line_at(&path, 4), 3);
        // Offset 1 is the newline after "a": the next element is on line 2.
        assert_eq!(line_at(&path, 1), 2);
    }
}
