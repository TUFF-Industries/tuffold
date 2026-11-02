//! Unquoted tab-separated report primitives.
//!
//! Philosopher/FragPipe write TSVs without CSV quoting, so fields are split on
//! raw tabs and every value is read as a trimmed string; typed coercion happens
//! in the mappers, where a bad value can name its file, line, and column.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::diagnostics::{InspectionIssue, IssueCode, quote_list};

const BOM: char = '\u{feff}';

/// A logical model field and the normalized headers that may carry it.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ColumnSpec {
    pub logical: &'static str,
    pub variants: &'static [&'static str],
    pub required: bool,
}

impl ColumnSpec {
    pub const fn required(logical: &'static str, variants: &'static [&'static str]) -> Self {
        Self {
            logical,
            variants,
            required: true,
        }
    }

    pub const fn optional(logical: &'static str, variants: &'static [&'static str]) -> Self {
        Self {
            logical,
            variants,
            required: false,
        }
    }
}

/// Logical field name to header column index.
pub(crate) type ResolvedColumns = BTreeMap<&'static str, usize>;

/// Collapse case and separator drift so `Q-value`, `Q_Value`, and `q value` match.
pub(crate) fn normalize_header(header: &str) -> String {
    let mut normalized = String::with_capacity(header.len());
    for ch in header.trim().to_lowercase().chars() {
        if ch.is_ascii_alphanumeric() {
            normalized.push(ch);
        } else if !normalized.ends_with(' ') {
            normalized.push(' ');
        }
    }
    normalized.trim().to_string()
}

fn strip_line_ending(line: &str) -> &str {
    let line = line.strip_suffix('\n').unwrap_or(line);
    line.strip_suffix('\r').unwrap_or(line)
}

/// The header names of a report, with a UTF-8 byte-order mark removed.
pub(crate) fn read_header(path: &Path) -> io::Result<Vec<String>> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut line = String::new();
    reader.read_line(&mut line)?;
    let line = strip_line_ending(&line);
    let line = line.strip_prefix(BOM).unwrap_or(line);
    Ok(line.split('\t').map(str::to_string).collect())
}

/// Map logical fields to header indices, pushing one problem per required miss.
///
/// When two headers normalize to the same text, the first one wins.
pub(crate) fn resolve_columns(
    file_name: &str,
    headers: &[String],
    specs: &[ColumnSpec],
    problems: &mut Vec<String>,
) -> ResolvedColumns {
    let mut lookup: HashMap<String, usize> = HashMap::new();
    for (index, header) in headers.iter().enumerate() {
        lookup.entry(normalize_header(header)).or_insert(index);
    }
    let mut resolved = ResolvedColumns::new();
    for spec in specs {
        match spec
            .variants
            .iter()
            .find_map(|variant| lookup.get(*variant))
        {
            Some(&index) => {
                resolved.insert(spec.logical, index);
            }
            None if spec.required => problems.push(format!(
                "{file_name}: no required column for '{}' (accepted headers: {}); found: {}",
                spec.logical,
                quote_list(spec.variants),
                quote_list(headers),
            )),
            None => {}
        }
    }
    resolved
}

/// Column-oriented values of a report: one trimmed string per data row.
#[derive(Debug, Default)]
pub(crate) struct Columns {
    values: HashMap<&'static str, Vec<String>>,
    rows: usize,
}

impl Columns {
    pub fn get(&self, logical: &str) -> Option<&[String]> {
        self.values.get(logical).map(Vec::as_slice)
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    /// The value at `row`, or the empty string for an absent optional column.
    pub fn value(&self, logical: &str, row: usize) -> &str {
        self.get(logical).map_or("", |column| column[row].as_str())
    }
}

/// Read the requested columns. Rows must be structurally sound (see
/// [`scan_report_structure`]); a short row yields empty trailing values.
pub(crate) fn read_columns(path: &Path, resolved: &ResolvedColumns) -> io::Result<Columns> {
    let mut reader = BufReader::new(File::open(path)?);
    let mut line = String::new();
    reader.read_line(&mut line)?; // header
    let mut values: HashMap<&'static str, Vec<String>> = resolved
        .keys()
        .map(|&logical| (logical, Vec::new()))
        .collect();
    let mut rows = 0;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let fields: Vec<&str> = strip_line_ending(&line).split('\t').collect();
        for (&logical, &index) in resolved {
            let value = fields.get(index).map_or("", |field| field.trim());
            values
                .get_mut(logical)
                .expect("every resolved column has a value vector")
                .push(value.to_string());
        }
        rows += 1;
    }
    Ok(Columns { values, rows })
}

/// Stream one report: validate UTF-8 and row shape, count data rows.
///
/// One pass, constant memory. Stops at the first structural problem: one
/// precise, line-numbered error is more actionable than a flood, and parsing
/// would stop there anyway.
pub(crate) fn scan_report_structure(path: &Path) -> (usize, Vec<InspectionIssue>) {
    let name = file_name(path);
    let file = match File::open(path) {
        Ok(file) => file,
        Err(error) => {
            return (
                0,
                vec![InspectionIssue::new(
                    IssueCode::UnreadableFile,
                    format!("{name}: file cannot be opened: {error}"),
                )],
            );
        }
    };
    let mut reader = BufReader::new(file);
    let mut raw = Vec::new();
    let fail = |rows, code, message: String| (rows, vec![InspectionIssue::new(code, message)]);
    let read = |reader: &mut BufReader<File>, raw: &mut Vec<u8>| {
        raw.clear();
        reader.read_until(b'\n', raw)
    };

    match read(&mut reader, &mut raw) {
        Ok(0) => {
            return fail(
                0,
                IssueCode::EmptyFile,
                format!("{name}: file is empty (0 bytes); expected a TSV report with a header row"),
            );
        }
        Ok(_) => {}
        Err(error) => {
            return fail(0, IssueCode::UnreadableFile, format!("{name}: {error}"));
        }
    }
    if raw.trim_ascii().is_empty() {
        return fail(
            0,
            IssueCode::EmptyFile,
            format!("{name}: first line is blank; expected a TSV header row"),
        );
    }
    if let Err(error) = std::str::from_utf8(&raw) {
        return fail(
            0,
            IssueCode::InvalidEncoding,
            format!("{name}: header line is not valid UTF-8: {error}"),
        );
    }
    let expected_fields = count_tabs(&raw) + 1;

    let mut rows = 0;
    let mut line_number = 1;
    loop {
        line_number += 1;
        match read(&mut reader, &mut raw) {
            Ok(0) => break,
            Ok(_) => {}
            Err(error) => {
                return fail(
                    rows,
                    IssueCode::UnreadableFile,
                    format!("{name}: line {line_number} cannot be read: {error}"),
                );
            }
        }
        if let Err(error) = std::str::from_utf8(&raw) {
            return fail(
                rows,
                IssueCode::InvalidEncoding,
                format!(
                    "{name}: line {line_number} is not valid UTF-8: {error}; \
                     reports must be UTF-8 encoded"
                ),
            );
        }
        if raw.trim_ascii().is_empty() {
            return fail(
                rows,
                IssueCode::MalformedTsv,
                format!("{name}: line {line_number} is blank; TSV rows must be non-empty"),
            );
        }
        let fields = count_tabs(&raw) + 1;
        if fields != expected_fields {
            return fail(
                rows,
                IssueCode::MalformedTsv,
                format!(
                    "{name}: line {line_number} has {fields} fields but the header has \
                     {expected_fields}; check for a missing or extra tab"
                ),
            );
        }
        rows += 1;
    }
    (rows, Vec::new())
}

fn count_tabs(raw: &[u8]) -> usize {
    raw.iter().filter(|&&byte| byte == b'\t').count()
}

pub(crate) fn file_name(path: &Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |name| name.to_string_lossy().into_owned(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_header_collapses_case_and_separators() {
        assert_eq!(
            normalize_header("PeptideProphet Q-value"),
            "peptideprophet q value"
        );
        assert_eq!(normalize_header("  Q_Value "), "q value");
        assert_eq!(normalize_header("\u{feff}Spectrum"), "spectrum");
        assert_eq!(normalize_header("M/Z"), "m z");
    }

    #[test]
    fn resolve_columns_reports_variants_and_found_headers() {
        const SPECS: &[ColumnSpec] = &[
            ColumnSpec::required("charge", &["charge", "assumed charge"]),
            ColumnSpec::optional("q_value", &["q value"]),
        ];
        let headers = vec!["Spectrum".to_string(), "Peptide".to_string()];
        let mut problems = Vec::new();
        let resolved = resolve_columns("psm.tsv", &headers, SPECS, &mut problems);
        assert!(resolved.is_empty());
        assert_eq!(
            problems,
            vec![
                "psm.tsv: no required column for 'charge' (accepted headers: 'charge', \
                 'assumed charge'); found: 'Spectrum', 'Peptide'"
                    .to_string()
            ]
        );
    }

    #[test]
    fn resolve_columns_first_duplicate_wins() {
        const SPECS: &[ColumnSpec] = &[ColumnSpec::required("q", &["q value"])];
        let headers = vec!["Q-Value".to_string(), "q_value".to_string()];
        let mut problems = Vec::new();
        let resolved = resolve_columns("x.tsv", &headers, SPECS, &mut problems);
        assert_eq!(resolved.get("q"), Some(&0));
    }

    #[test]
    fn scan_counts_rows_and_accepts_crlf_and_missing_final_newline() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.tsv");
        std::fs::write(&path, "a\tb\r\n1\t2\r\n3\t4").unwrap();
        assert_eq!(scan_report_structure(&path), (2, vec![]));
    }

    #[test]
    fn read_columns_trims_values_and_strips_crlf() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.tsv");
        std::fs::write(&path, "a\tb\r\n 1 \tx\r\n2\ty\r\n").unwrap();
        let resolved = ResolvedColumns::from([("a", 0), ("b", 1)]);
        let columns = read_columns(&path, &resolved).unwrap();
        assert_eq!(columns.rows(), 2);
        assert_eq!(columns.get("a").unwrap(), ["1", "2"]);
        assert_eq!(columns.value("b", 1), "y");
        assert_eq!(columns.value("absent", 0), "");
    }
}
