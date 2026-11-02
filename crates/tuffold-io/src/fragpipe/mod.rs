//! FragPipe/Philosopher results bundles.
//!
//! A bundle is a directory of TSV reports:
//!
//! ```text
//! <bundle>/
//!     psm.tsv                one row per PSM, all experiments combined
//!     modified_peptide.tsv   peptide sequences including modifications
//!     ion.tsv                PSMs collapsed to peptide ions
//!     combined_peptide.tsv   stripped peptide sequences, all experiments
//!     combined_protein.tsv   protein groups, all experiments
//!     <experiment>/          one subfolder per experiment (sample)
//!     *.fas | *.fasta | *.fa optional protein database
//! ```
//!
//! [`inspect_bundle`] owns the layout rules and returns typed diagnostics
//! without ever failing on a bad bundle; [`parse_bundle`] gates on it and maps
//! columns onto the model. Column names drift between Philosopher versions, so
//! headers are normalized (lower-cased, separators collapsed) and matched
//! against the variant tables below.

mod inspection;
mod msfragger;
mod parser;

pub use inspection::{BundleInspection, inspect_bundle};
pub use msfragger::{SearchFile, find_search_files, import_search_files};
pub use parser::{ParseError, import_bundle, import_inspected};

use tuffold_core::BioSample;

use crate::tsv::ColumnSpec;

pub const PSM_FILE: &str = "psm.tsv";
pub const MODIFIED_PEPTIDE_FILE: &str = "modified_peptide.tsv";
pub const ION_FILE: &str = "ion.tsv";
pub const PEPTIDE_FILE: &str = "combined_peptide.tsv";
pub const PROTEIN_GROUP_FILE: &str = "combined_protein.tsv";
pub const BUNDLE_FILES: &[&str] = &[
    PSM_FILE,
    MODIFIED_PEPTIDE_FILE,
    ION_FILE,
    PEPTIDE_FILE,
    PROTEIN_GROUP_FILE,
];

pub(crate) const PSM_COLUMNS: &[ColumnSpec] = &[
    ColumnSpec::required("psm_id", &["spectrum"]),
    ColumnSpec::required("sample_id", &["spectrum file", "sample"]),
    ColumnSpec::required(
        "sequence",
        &["peptide", "peptide sequence", "stripped peptide"],
    ),
    ColumnSpec::required("charge", &["charge", "assumed charge"]),
    ColumnSpec::required("score", &["hyperscore", "score", "xcorr"]),
    ColumnSpec::required("protein_id", &["protein id"]),
    ColumnSpec::optional("mapped_proteins", &["mapped proteins"]),
    ColumnSpec::optional("assigned_modifications", &["assigned modifications"]),
    ColumnSpec::optional("retention", &["retention", "retention time"]),
    ColumnSpec::optional("observed_mz", &["observed m z"]),
    ColumnSpec::optional("calculated_mz", &["calculated m z"]),
    ColumnSpec::optional("protein_start", &["protein start"]),
    ColumnSpec::optional("protein_end", &["protein end"]),
    ColumnSpec::optional("prev_aa", &["prev aa"]),
    ColumnSpec::optional("next_aa", &["next aa"]),
    ColumnSpec::optional("protein_description", &["protein description"]),
    ColumnSpec::optional("nextscore", &["nextscore"]),
    ColumnSpec::optional("expectation", &["expectation", "expect"]),
];

pub(crate) const PEPTIDE_COLUMNS: &[ColumnSpec] = &[
    ColumnSpec::required(
        "sequence",
        &["peptide sequence", "peptide", "stripped peptide"],
    ),
    ColumnSpec::required("protein_id", &["protein id"]),
    ColumnSpec::optional("mapped_proteins", &["mapped proteins"]),
];

pub(crate) const PROTEIN_GROUP_COLUMNS: &[ColumnSpec] = &[
    ColumnSpec::required("protein_id", &["protein id"]),
    ColumnSpec::optional(
        "indistinguishable_proteins",
        &["indistinguishable proteins", "indistinguishable protein"],
    ),
];

/// The reports whose columns are mapped, with their column tables.
pub(crate) const MAPPED_REPORTS: &[(&str, &[ColumnSpec])] = &[
    (PSM_FILE, PSM_COLUMNS),
    (PEPTIDE_FILE, PEPTIDE_COLUMNS),
    (PROTEIN_GROUP_FILE, PROTEIN_GROUP_COLUMNS),
];

/// Combined reports prefix identification-file names with `interact-`, and
/// some versions keep the raw-file extension; both are stripped when linking a
/// PSM to its sample subfolder.
const IDENTIFICATION_FILE_PREFIXES: &[&str] = &["interact-"];
const IDENTIFICATION_FILE_EXTENSIONS: &[&str] = &[
    ".pep.xml",
    ".pep.xml.gz",
    ".raw",
    ".mzml",
    ".mzxml",
    ".mgf",
    ".ms2",
    ".ms1",
    ".d",
];

/// Philosopher joins protein lists with spaces, commas, or semicolons.
pub(crate) fn split_protein_list(value: &str) -> impl Iterator<Item = &str> {
    value
        .split(|ch: char| ch == ';' || ch == ',' || ch.is_whitespace())
        .filter(|item| !item.is_empty())
}

/// Derive a sample id from a `Spectrum File` value: basename, minus the
/// `interact-` prefix and one known identification/raw-file extension.
pub(crate) fn sample_id_from_spectrum_file(value: &str) -> String {
    let normalized = value.trim().replace('\\', "/");
    let mut cleaned = normalized.rsplit('/').next().unwrap_or_default();
    if let Some(prefix) = IDENTIFICATION_FILE_PREFIXES
        .iter()
        .find(|prefix| starts_with_ignore_ascii_case(cleaned, prefix))
    {
        cleaned = &cleaned[prefix.len()..];
    }
    if let Some(extension) = IDENTIFICATION_FILE_EXTENSIONS
        .iter()
        .find(|extension| ends_with_ignore_ascii_case(cleaned, extension))
    {
        cleaned = &cleaned[..cleaned.len() - extension.len()];
    }
    cleaned.to_string()
}

fn starts_with_ignore_ascii_case(value: &str, prefix: &str) -> bool {
    value
        .get(..prefix.len())
        .is_some_and(|head| head.eq_ignore_ascii_case(prefix))
}

fn ends_with_ignore_ascii_case(value: &str, suffix: &str) -> bool {
    value
        .len()
        .checked_sub(suffix.len())
        .and_then(|start| value.get(start..))
        .is_some_and(|tail| tail.eq_ignore_ascii_case(suffix))
}

pub(crate) const UNKNOWN_CONDITION: &str = "unknown";

/// An experiment's biosample: its category is the name's condition prefix.
pub(crate) fn biosample(name: &str) -> BioSample {
    let condition = condition_from_sample_name(name);
    BioSample {
        name: name.to_string(),
        category: (condition != UNKNOWN_CONDITION).then_some(condition),
    }
}

/// Experiment names conventionally start with the condition: `control-1` ->
/// `control`. Names without a `-`/`_` separator get `unknown`.
pub fn condition_from_sample_name(name: &str) -> String {
    match name.split_once(['-', '_']) {
        Some((prefix, _)) => prefix.to_string(),
        None => UNKNOWN_CONDITION.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_id_strips_paths_prefix_and_one_extension() {
        assert_eq!(sample_id_from_spectrum_file("control-1"), "control-1");
        assert_eq!(
            sample_id_from_spectrum_file("interact-control-1.raw"),
            "control-1"
        );
        assert_eq!(
            sample_id_from_spectrum_file("C:\\runs\\Interact-A.pep.xml"),
            "A"
        );
        assert_eq!(sample_id_from_spectrum_file("/data/x/B.pep.xml.gz"), "B");
        assert_eq!(sample_id_from_spectrum_file("run.mzML"), "run");
        assert_eq!(sample_id_from_spectrum_file("run.mzML.raw"), "run.mzML");
    }

    #[test]
    fn protein_lists_split_on_any_separator() {
        let ids: Vec<&str> = split_protein_list(" P1, P2;P3  P4 ").collect();
        assert_eq!(ids, ["P1", "P2", "P3", "P4"]);
        assert_eq!(split_protein_list("").count(), 0);
    }

    #[test]
    fn condition_is_prefix_before_first_separator() {
        assert_eq!(condition_from_sample_name("control-1"), "control");
        assert_eq!(condition_from_sample_name("treat_a-2"), "treat");
        assert_eq!(condition_from_sample_name("single"), "unknown");
    }
}
