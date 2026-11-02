//! Read-only bundle inspection.
//!
//! [`inspect_bundle`] never fails on a bad bundle: every problem comes back as
//! a typed issue, errors meaning "cannot import" and warnings meaning "imports
//! with a stated limitation". Messages state consequences (which rows stay
//! loaded, what becomes unavailable) so users are never left guessing.
//!
//! Orderings are deterministic: reports in [`BUNDLE_FILES`] order, files and
//! folders sorted, FASTA candidates by case-insensitive file name. Hidden
//! (dot-prefixed) entries such as `.DS_Store`, `.meta/`, or macOS `._*`
//! resource forks are not part of a bundle and are ignored.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use serde::Serialize;

use super::{
    BUNDLE_FILES, MAPPED_REPORTS, PEPTIDE_FILE, PROTEIN_GROUP_FILE, PSM_FILE,
    sample_id_from_spectrum_file, split_protein_list,
};
use crate::diagnostics::{InspectionIssue, IssueCode, format_names};
use crate::fasta::{
    FASTA_EXTENSIONS, FastaRecord, build_identifier_index, find_fasta_candidates, parse_fasta,
};
use crate::tsv::{
    Columns, ResolvedColumns, read_columns, read_header, resolve_columns, scan_report_structure,
};

const NAME_LIMIT: usize = 10;

/// Read-only snapshot of a bundle's layout plus every finding.
#[derive(Debug, Clone, Serialize)]
pub struct BundleInspection {
    pub path: PathBuf,
    pub required_files: Vec<String>,
    pub missing_files: Vec<String>,
    /// Every non-hidden top-level file, sorted.
    pub discovered_files: Vec<String>,
    /// Every non-hidden top-level directory, sorted; one sample each.
    pub experiment_folders: Vec<String>,
    /// FASTA candidates in selection order; the first is `selected_fasta`.
    pub fasta_candidates: Vec<String>,
    pub selected_fasta: Option<PathBuf>,
    pub fasta_record_count: Option<usize>,
    /// Data-row counts (header excluded) of structurally valid reports.
    pub row_counts: BTreeMap<String, usize>,
    pub warnings: Vec<InspectionIssue>,
    pub errors: Vec<InspectionIssue>,
    /// Column resolution of each mapped report, reused by the parser.
    #[serde(skip)]
    pub(crate) resolved: BTreeMap<&'static str, ResolvedColumns>,
    /// Header row of each mapped report, for verbatim score names.
    #[serde(skip)]
    pub(crate) headers: BTreeMap<&'static str, Vec<String>>,
}

impl BundleInspection {
    fn empty(path: PathBuf, errors: Vec<InspectionIssue>) -> Self {
        Self {
            path,
            required_files: BUNDLE_FILES.iter().map(|name| name.to_string()).collect(),
            missing_files: Vec::new(),
            discovered_files: Vec::new(),
            experiment_folders: Vec::new(),
            fasta_candidates: Vec::new(),
            selected_fasta: None,
            fasta_record_count: None,
            row_counts: BTreeMap::new(),
            warnings: Vec::new(),
            errors,
            resolved: BTreeMap::new(),
            headers: BTreeMap::new(),
        }
    }

    /// The verbatim header a logical field resolved to in a mapped report.
    pub(crate) fn header(&self, file: &str, logical: &str) -> Option<&str> {
        let index = *self.resolved.get(file)?.get(logical)?;
        self.headers.get(file)?.get(index).map(String::as_str)
    }

    /// `true` when the bundle can be parsed.
    pub fn ok(&self) -> bool {
        self.errors.is_empty()
    }

    /// Discovered files that are not required reports.
    pub fn optional_files(&self) -> Vec<&str> {
        self.discovered_files
            .iter()
            .filter(|name| !BUNDLE_FILES.contains(&name.as_str()))
            .map(String::as_str)
            .collect()
    }

    /// Render the errors as one actionable message naming the bundle path.
    pub fn format_errors(&self) -> String {
        let mut rendered = format!("Invalid FragPipe bundle at '{}':", self.path.display());
        for issue in &self.errors {
            let _ = write!(rendered, "\n  [{}] {}", issue.code, issue.message);
        }
        rendered
    }
}

fn is_hidden(name: &str) -> bool {
    name.starts_with('.')
}

/// Inspect a FragPipe bundle and return typed diagnostics.
pub fn inspect_bundle(path: impl AsRef<Path>) -> BundleInspection {
    let bundle = path.as_ref().to_path_buf();
    if !bundle.exists() {
        let message = format!("bundle path does not exist: '{}'", bundle.display());
        return BundleInspection::empty(
            bundle,
            vec![InspectionIssue::new(IssueCode::InvalidPath, message)],
        );
    }
    if !bundle.is_dir() {
        let message = format!(
            "bundle path must be a directory, but is a file: '{}'",
            bundle.display()
        );
        return BundleInspection::empty(
            bundle,
            vec![InspectionIssue::new(IssueCode::InvalidPath, message)],
        );
    }
    let (files, folders) = match list_entries(&bundle) {
        Ok(entries) => entries,
        Err(error) => {
            let message = format!("bundle directory cannot be read: {error}");
            return BundleInspection::empty(
                bundle,
                vec![InspectionIssue::new(IssueCode::UnreadableFile, message)],
            );
        }
    };

    let mut inspection = BundleInspection {
        discovered_files: files,
        experiment_folders: folders,
        ..BundleInspection::empty(bundle, Vec::new())
    };
    check_layout(&mut inspection);
    let sound_reports = check_report_structure(&mut inspection);
    check_columns(&mut inspection, &sound_reports);
    let fasta_records = check_fasta(&mut inspection);
    check_cross_reports(&mut inspection, fasta_records.as_deref());
    inspection
}

fn list_entries(bundle: &Path) -> std::io::Result<(Vec<String>, Vec<String>)> {
    let mut files = Vec::new();
    let mut folders = Vec::new();
    for entry in std::fs::read_dir(bundle)? {
        let path = entry?.path();
        let name = crate::tsv::file_name(&path);
        if is_hidden(&name) {
            continue;
        }
        if path.is_dir() {
            folders.push(name);
        } else if path.is_file() {
            files.push(name);
        }
    }
    files.sort();
    folders.sort();
    Ok((files, folders))
}

fn check_layout(inspection: &mut BundleInspection) {
    for &name in BUNDLE_FILES {
        if !inspection.discovered_files.iter().any(|file| file == name) {
            inspection.missing_files.push(name.to_string());
            inspection.errors.push(InspectionIssue::new(
                IssueCode::MissingRequiredFile,
                format!("missing required report '{name}'"),
            ));
        }
    }
    if inspection.experiment_folders.is_empty() {
        inspection.errors.push(InspectionIssue::new(
            IssueCode::NoExperimentFolders,
            "no experiment subfolder found (exactly one per experiment is required)",
        ));
    }
}

/// Scan every present report; returns the structurally sound ones.
fn check_report_structure(inspection: &mut BundleInspection) -> BTreeSet<&'static str> {
    let mut sound = BTreeSet::new();
    for &name in BUNDLE_FILES {
        if inspection
            .missing_files
            .iter()
            .any(|missing| missing == name)
        {
            continue;
        }
        let (rows, issues) = scan_report_structure(&inspection.path.join(name));
        if !issues.is_empty() {
            inspection.errors.extend(issues);
            continue;
        }
        sound.insert(name);
        inspection.row_counts.insert(name.to_string(), rows);
        if rows == 0 {
            inspection.warnings.push(InspectionIssue::new(
                IssueCode::EmptyReport,
                format!("{name}: report has a header but no data rows"),
            ));
        }
    }
    sound
}

fn check_columns(inspection: &mut BundleInspection, sound_reports: &BTreeSet<&'static str>) {
    for &(file_name, specs) in MAPPED_REPORTS {
        if !sound_reports.contains(file_name) {
            continue;
        }
        let headers = match read_header(&inspection.path.join(file_name)) {
            Ok(headers) => headers,
            Err(error) => {
                inspection.errors.push(InspectionIssue::new(
                    IssueCode::UnreadableFile,
                    format!("{file_name}: header cannot be read: {error}"),
                ));
                continue;
            }
        };
        let mut problems = Vec::new();
        let resolved = resolve_columns(file_name, &headers, specs, &mut problems);
        inspection.errors.extend(
            problems
                .into_iter()
                .map(|problem| InspectionIssue::new(IssueCode::MissingRequiredColumns, problem)),
        );
        inspection.resolved.insert(file_name, resolved);
        inspection.headers.insert(file_name, headers);
    }

    if inspection.resolved.contains_key(PSM_FILE) {
        inspection.warnings.push(InspectionIssue::new(
            IssueCode::FilteredPsms,
            format!(
                "{PSM_FILE} is FragPipe's filtered list (cut at its FDR, decoys removed), \
                 so Tuffold's FDR on it is not an independent check; keep MSFragger's \
                 .pepXML files in the folder and Tuffold validates from them instead"
            ),
        ));
    }
}

fn check_fasta(inspection: &mut BundleInspection) -> Option<Vec<FastaRecord>> {
    let candidates: Vec<PathBuf> = match find_fasta_candidates(&inspection.path) {
        Ok(found) => found
            .into_iter()
            .filter(|path| !is_hidden(&crate::tsv::file_name(path)))
            .collect(),
        Err(error) => {
            inspection.errors.push(InspectionIssue::new(
                IssueCode::UnreadableFile,
                format!("FASTA discovery failed: {error}"),
            ));
            return None;
        }
    };
    inspection.fasta_candidates = candidates
        .iter()
        .map(|path| crate::tsv::file_name(path))
        .collect();
    let Some(selected) = candidates.into_iter().next() else {
        let extensions: Vec<String> = FASTA_EXTENSIONS
            .iter()
            .map(|ext| format!(".{ext}"))
            .collect();
        inspection.warnings.push(InspectionIssue::new(
            IssueCode::MissingFasta,
            format!(
                "no protein FASTA found at the top level of '{}' (accepted extensions: {}); \
                 sequence coverage will be unavailable",
                inspection.path.display(),
                extensions.join(", ")
            ),
        ));
        return None;
    };
    let selected_name = crate::tsv::file_name(&selected);
    if inspection.fasta_candidates.len() > 1 {
        inspection.warnings.push(InspectionIssue::new(
            IssueCode::AmbiguousFasta,
            format!(
                "multiple protein FASTA candidates found ({}); selected '{selected_name}' \
                 (alphabetically first file name, case-insensitive)",
                inspection.fasta_candidates.join(", ")
            ),
        ));
    }
    inspection.selected_fasta = Some(selected.clone());
    match parse_fasta(&selected) {
        Ok(records) => {
            inspection.fasta_record_count = Some(records.len());
            Some(records)
        }
        Err(error) if error.kind() == std::io::ErrorKind::InvalidData => {
            inspection.errors.push(InspectionIssue::new(
                IssueCode::InvalidEncoding,
                format!(
                    "{selected_name}: file is not valid UTF-8: {error}; \
                     re-export the protein database with UTF-8 encoding"
                ),
            ));
            None
        }
        Err(error) => {
            inspection.errors.push(InspectionIssue::new(
                IssueCode::UnreadableFile,
                format!("{selected_name}: file cannot be read: {error}"),
            ));
            None
        }
    }
}

/// Read the requested resolved columns of a mapped report, if it resolved.
fn report_columns(
    inspection: &mut BundleInspection,
    file_name: &'static str,
    logicals: &[&'static str],
) -> Option<Columns> {
    let resolved = inspection.resolved.get(file_name)?;
    let wanted: ResolvedColumns = logicals
        .iter()
        .filter_map(|logical| resolved.get(logical).map(|&index| (*logical, index)))
        .collect();
    if wanted.is_empty() {
        return None;
    }
    match read_columns(&inspection.path.join(file_name), &wanted) {
        Ok(columns) => Some(columns),
        Err(error) => {
            inspection.errors.push(InspectionIssue::new(
                IssueCode::UnreadableFile,
                format!("{file_name}: columns could not be read: {error}"),
            ));
            None
        }
    }
}

fn collect_proteins<'a>(columns: &'a Columns, logicals: &[&str], into: &mut BTreeSet<&'a str>) {
    for logical in logicals {
        for value in columns.get(logical).unwrap_or_default() {
            into.extend(split_protein_list(value));
        }
    }
}

/// Cross-report consistency checks; every finding is a warning.
fn check_cross_reports(inspection: &mut BundleInspection, fasta_records: Option<&[FastaRecord]>) {
    let psm = report_columns(
        inspection,
        PSM_FILE,
        &["sample_id", "protein_id", "mapped_proteins", "sequence"],
    );
    let peptide = report_columns(
        inspection,
        PEPTIDE_FILE,
        &["sequence", "protein_id", "mapped_proteins"],
    );
    let protein = report_columns(
        inspection,
        PROTEIN_GROUP_FILE,
        &["protein_id", "indistinguishable_proteins"],
    );
    let mut warnings = Vec::new();

    if let Some(sample_ids) = psm.as_ref().and_then(|columns| columns.get("sample_id")) {
        warnings.extend(check_samples(&inspection.experiment_folders, sample_ids));
    }

    let mut group_members = BTreeSet::new();
    if let Some(columns) = &protein {
        collect_proteins(
            columns,
            &["protein_id", "indistinguishable_proteins"],
            &mut group_members,
        );
    }
    let mut referenced = BTreeSet::new();
    for columns in [&psm, &peptide].into_iter().flatten() {
        collect_proteins(columns, &["protein_id", "mapped_proteins"], &mut referenced);
    }
    if !group_members.is_empty() {
        let unmatched: Vec<&str> = referenced.difference(&group_members).copied().collect();
        if !unmatched.is_empty() {
            warnings.push(InspectionIssue::new(
                IssueCode::UnmatchedProteinIds,
                format!(
                    "protein IDs referenced by {PSM_FILE} or {PEPTIDE_FILE} rows ('Protein ID' \
                     or 'Mapped Proteins' columns) but present in no {PROTEIN_GROUP_FILE} group \
                     ({}): {}; the referencing rows remain loaded and usable",
                    unmatched.len(),
                    format_names(&unmatched, NAME_LIMIT)
                ),
            ));
        }
    }

    let psm_sequences = psm.as_ref().and_then(|columns| columns.get("sequence"));
    let peptide_sequences = peptide.as_ref().and_then(|columns| columns.get("sequence"));
    if let Some(peptide_sequences) = peptide_sequences
        && let Some(issue) = check_duplicate_peptides(peptide_sequences)
    {
        inspection.errors.push(issue);
    }
    if let (Some(psm_sequences), Some(peptide_sequences)) = (psm_sequences, peptide_sequences) {
        let supported: BTreeSet<&str> = non_blank(psm_sequences).collect();
        let listed: BTreeSet<&str> = non_blank(peptide_sequences).collect();
        let unlisted: Vec<&str> = supported.difference(&listed).copied().collect();
        if !unlisted.is_empty() {
            warnings.push(InspectionIssue::new(
                IssueCode::UnlistedPsmPeptides,
                format!(
                    "peptide sequences of {PSM_FILE} rows absent from {PEPTIDE_FILE} ({}): {}; \
                     the PSM rows remain loaded, but their peptide_id links to no peptide record",
                    unlisted.len(),
                    format_names(&unlisted, NAME_LIMIT)
                ),
            ));
        }
        let unsupported: Vec<&str> = listed.difference(&supported).copied().collect();
        if !unsupported.is_empty() {
            warnings.push(InspectionIssue::new(
                IssueCode::UnmatchedPeptides,
                format!(
                    "peptides in {PEPTIDE_FILE} with no supporting PSM row in {PSM_FILE} ({}): \
                     {}; they remain loaded with no linked PSMs (empty psm_ids)",
                    unsupported.len(),
                    format_names(&unsupported, NAME_LIMIT)
                ),
            ));
        }
    }

    if let (Some(records), Some(selected)) = (fasta_records, &inspection.selected_fasta)
        && !group_members.is_empty()
    {
        let known = build_identifier_index(records);
        let unmatched: Vec<&str> = group_members
            .iter()
            .filter(|member| !known.contains_key(**member))
            .copied()
            .collect();
        if !unmatched.is_empty() {
            warnings.push(InspectionIssue::new(
                IssueCode::UnmatchedFastaProteins,
                format!(
                    "protein-group members with no record in the selected FASTA '{}' ({}): \
                         {}; they stay grouped but get no sequence or coverage",
                    crate::tsv::file_name(selected),
                    unmatched.len(),
                    format_names(&unmatched, NAME_LIMIT)
                ),
            ));
        }
    }
    inspection.warnings.extend(warnings);
}

/// Peptide ids derive from the sequence, so a repeated sequence would give two
/// records one id; that is an error, reported with every offending line.
fn check_duplicate_peptides(sequences: &[String]) -> Option<InspectionIssue> {
    let mut lines_by_sequence: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (row, sequence) in sequences.iter().enumerate() {
        if !sequence.is_empty() {
            lines_by_sequence.entry(sequence).or_default().push(row + 2);
        }
    }
    let duplicates: Vec<String> = lines_by_sequence
        .into_iter()
        .filter(|(_, lines)| lines.len() > 1)
        .map(|(sequence, lines)| format!("{sequence} (lines {})", join_lines(&lines)))
        .collect();
    (!duplicates.is_empty()).then(|| {
        InspectionIssue::new(
            IssueCode::DuplicatePeptides,
            format!(
                "{PEPTIDE_FILE} lists the same peptide sequence on several rows ({}): {}; \
                 each sequence must appear once",
                duplicates.len(),
                format_names(&duplicates, NAME_LIMIT)
            ),
        )
    })
}

fn join_lines(lines: &[usize]) -> String {
    match lines.split_last() {
        Some((last, rest)) if !rest.is_empty() => {
            let rest: Vec<String> = rest.iter().map(usize::to_string).collect();
            format!("{} and {last}", rest.join(", "))
        }
        _ => lines.iter().map(usize::to_string).collect(),
    }
}

fn non_blank(values: &[String]) -> impl Iterator<Item = &str> {
    values
        .iter()
        .map(String::as_str)
        .filter(|value| !value.trim().is_empty())
}

fn check_samples(folders: &[String], sample_column: &[String]) -> Vec<InspectionIssue> {
    let folder_set: BTreeSet<&str> = folders.iter().map(String::as_str).collect();
    let psm_samples: BTreeSet<String> = non_blank(sample_column)
        .map(sample_id_from_spectrum_file)
        .collect();
    let mut issues = Vec::new();
    let orphans: Vec<&str> = folder_set
        .iter()
        .filter(|folder| !psm_samples.contains(**folder))
        .copied()
        .collect();
    if !orphans.is_empty() {
        issues.push(InspectionIssue::new(
            IssueCode::UnmatchedSamples,
            format!(
                "experiment subfolders with no PSM rows in {PSM_FILE}: {}; each listed subfolder \
                 is still parsed as a sample with zero PSMs",
                format_names(&orphans, NAME_LIMIT)
            ),
        ));
    }
    let unknown: Vec<&str> = psm_samples
        .iter()
        .map(String::as_str)
        .filter(|sample| !folder_set.contains(sample))
        .collect();
    if !unknown.is_empty() {
        issues.push(InspectionIssue::new(
            IssueCode::UnmatchedSamples,
            format!(
                "PSM sample ids in {PSM_FILE} with no matching experiment subfolder: {}; the PSM \
                 rows remain loaded with those sample ids",
                format_names(&unknown, NAME_LIMIT)
            ),
        ));
    }
    issues
}
