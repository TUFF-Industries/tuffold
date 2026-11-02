//! Typed diagnostics shared by importers.
//!
//! Issue codes are a stable vocabulary: callers (the CLI, the desktop app,
//! later the MCP server) branch on [`IssueCode`], never on message text.
//! Messages name the affected file, the logical field, the accepted variants,
//! and the consequence for the loaded data.

use std::fmt;

use serde::Serialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum IssueCode {
    // Errors: the input cannot be imported.
    /// The path does not exist or is not a directory.
    InvalidPath,
    /// A required report is absent.
    MissingRequiredFile,
    /// The bundle has no experiment subfolder.
    NoExperimentFolders,
    /// A required report exists but has no header row.
    EmptyFile,
    /// A report or the selected FASTA is not valid UTF-8.
    InvalidEncoding,
    /// A data row's field count does not match the header, or a row is blank.
    MalformedTsv,
    /// A report lacks a column for a required logical field.
    MissingRequiredColumns,
    /// The bundle or one of its files cannot be read.
    UnreadableFile,
    /// The peptide report lists a sequence more than once.
    DuplicatePeptides,

    // Warnings: import proceeds with a stated limitation.
    /// No protein FASTA; sequence coverage is unavailable.
    MissingFasta,
    /// Several FASTA candidates; a deterministic one is selected.
    AmbiguousFasta,
    /// A report has a header but zero data rows.
    EmptyReport,
    /// `psm.tsv` is already filtered, so validating it is not independent.
    FilteredPsms,
    /// Subfolders without PSMs, or PSM sample ids without a subfolder.
    UnmatchedSamples,
    /// Referenced protein ids that appear in no protein group.
    UnmatchedProteinIds,
    /// Peptides with no supporting PSM row.
    UnmatchedPeptides,
    /// PSM sequences absent from the peptide report.
    UnlistedPsmPeptides,
    /// Protein-group members absent from the selected FASTA.
    UnmatchedFastaProteins,
}

impl IssueCode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::InvalidPath => "invalid-path",
            Self::MissingRequiredFile => "missing-required-file",
            Self::NoExperimentFolders => "no-experiment-folders",
            Self::EmptyFile => "empty-file",
            Self::InvalidEncoding => "invalid-encoding",
            Self::MalformedTsv => "malformed-tsv",
            Self::MissingRequiredColumns => "missing-required-columns",
            Self::UnreadableFile => "unreadable-file",
            Self::DuplicatePeptides => "duplicate-peptides",
            Self::MissingFasta => "missing-fasta",
            Self::AmbiguousFasta => "ambiguous-fasta",
            Self::EmptyReport => "empty-report",
            Self::FilteredPsms => "filtered-psms",
            Self::UnmatchedSamples => "unmatched-samples",
            Self::UnmatchedProteinIds => "unmatched-protein-ids",
            Self::UnmatchedPeptides => "unmatched-peptides",
            Self::UnlistedPsmPeptides => "unlisted-psm-peptides",
            Self::UnmatchedFastaProteins => "unmatched-fasta-proteins",
        }
    }

    pub fn severity(self) -> Severity {
        match self {
            Self::InvalidPath
            | Self::MissingRequiredFile
            | Self::NoExperimentFolders
            | Self::EmptyFile
            | Self::InvalidEncoding
            | Self::MalformedTsv
            | Self::MissingRequiredColumns
            | Self::UnreadableFile
            | Self::DuplicatePeptides => Severity::Error,
            Self::MissingFasta
            | Self::AmbiguousFasta
            | Self::EmptyReport
            | Self::FilteredPsms
            | Self::UnmatchedSamples
            | Self::UnmatchedProteinIds
            | Self::UnmatchedPeptides
            | Self::UnlistedPsmPeptides
            | Self::UnmatchedFastaProteins => Severity::Warning,
        }
    }
}

impl fmt::Display for IssueCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Severity {
    Error,
    Warning,
}

/// One diagnostic finding: a stable code plus a human-readable message.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct InspectionIssue {
    pub code: IssueCode,
    pub message: String,
}

impl InspectionIssue {
    pub fn new(code: IssueCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

/// Join names for a message, eliding everything past `limit`.
pub(crate) fn format_names<S: AsRef<str>>(names: &[S], limit: usize) -> String {
    let shown: Vec<&str> = names.iter().take(limit).map(AsRef::as_ref).collect();
    let mut rendered = shown.join(", ");
    if names.len() > limit {
        rendered.push_str(&format!(", ... and {} more", names.len() - limit));
    }
    rendered
}

/// Quote names the way the messages render headers and variants: `'a', 'b'`.
pub(crate) fn quote_list<S: AsRef<str>>(names: &[S]) -> String {
    names
        .iter()
        .map(|name| format!("'{}'", name.as_ref()))
        .collect::<Vec<_>>()
        .join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serialized_code_matches_as_str() {
        for code in [IssueCode::FilteredPsms, IssueCode::UnmatchedFastaProteins] {
            let json = serde_json::to_value(code).unwrap();
            assert_eq!(json, serde_json::Value::String(code.as_str().into()));
        }
    }

    #[test]
    fn format_names_elides_past_limit() {
        let names: Vec<String> = (1..=12).map(|i| format!("n{i}")).collect();
        assert_eq!(
            format_names(&names, 10),
            "n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, ... and 2 more"
        );
        assert_eq!(format_names(&["a", "b"], 10), "a, b");
    }
}
