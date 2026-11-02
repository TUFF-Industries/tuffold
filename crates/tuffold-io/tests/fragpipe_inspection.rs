//! Bundle inspection tests.
//!
//! Expectations are derived by hand from fixtures/fragpipe_synthetic and the
//! inline bundles, never copied from inspector output. The fixture holds at the
//! top level: psm.tsv (5 rows), modified_peptide.tsv / ion.tsv /
//! combined_peptide.tsv (4 rows each), combined_protein.tsv (2 rows), and
//! protein.fas with 5 records (P00001, P00002, P00003, the rev_ decoy of
//! P00001, and the unidentified P00004), plus subfolders control-1/ and
//! treatment-1/.

mod common;

use std::collections::BTreeMap;
use std::fs;

use common::*;
use tuffold_io::fragpipe::import_inspected;
use tuffold_io::{BundleInspection, InspectionIssue, IssueCode, import_bundle, inspect_bundle};

fn codes(issues: &[InspectionIssue]) -> Vec<IssueCode> {
    issues.iter().map(|issue| issue.code).collect()
}

fn with_code(issues: &[InspectionIssue], code: IssueCode) -> Vec<&str> {
    issues
        .iter()
        .filter(|issue| issue.code == code)
        .map(|issue| issue.message.as_str())
        .collect()
}

fn assert_clean(inspection: &BundleInspection) {
    assert!(
        inspection.errors.is_empty(),
        "unexpected errors: {:?}",
        inspection.errors
    );
}

#[test]
fn valid_synthetic_bundle() {
    let bundle = fixture_bundle();
    let inspection = inspect_bundle(&bundle);

    assert!(inspection.ok());
    assert_eq!(codes(&inspection.warnings), [IssueCode::FilteredPsms]);
    assert_eq!(
        inspection.required_files,
        [
            "psm.tsv",
            "modified_peptide.tsv",
            "ion.tsv",
            "combined_peptide.tsv",
            "combined_protein.tsv"
        ]
    );
    assert!(inspection.missing_files.is_empty());
    let expected_rows = BTreeMap::from([
        ("psm.tsv".to_string(), 5),
        ("modified_peptide.tsv".to_string(), 4),
        ("ion.tsv".to_string(), 4),
        ("combined_peptide.tsv".to_string(), 4),
        ("combined_protein.tsv".to_string(), 2),
    ]);
    assert_eq!(inspection.row_counts, expected_rows);
    assert_eq!(inspection.experiment_folders, ["control-1", "treatment-1"]);
    assert_eq!(inspection.fasta_candidates, ["protein.fas"]);
    assert_eq!(inspection.selected_fasta, Some(bundle.join("protein.fas")));
    assert_eq!(inspection.fasta_record_count, Some(5));
    assert_eq!(inspection.optional_files(), ["protein.fas"]);
}

#[test]
fn missing_required_files() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::remove_file(bundle.join("combined_protein.tsv")).unwrap();
    fs::remove_file(bundle.join("ion.tsv")).unwrap();

    let inspection = inspect_bundle(&bundle);

    let missing = with_code(&inspection.errors, IssueCode::MissingRequiredFile);
    assert_eq!(missing.len(), 2);
    assert!(missing.iter().any(|m| m.contains("combined_protein.tsv")));
    assert!(missing.iter().any(|m| m.contains("ion.tsv")));
    let counted: Vec<&str> = inspection.row_counts.keys().map(String::as_str).collect();
    assert_eq!(
        counted,
        ["combined_peptide.tsv", "modified_peptide.tsv", "psm.tsv"]
    );
    assert_eq!(codes(&inspection.warnings), [IssueCode::FilteredPsms]);

    let formatted = inspection.format_errors();
    assert!(formatted.contains(&bundle.display().to_string()));
    assert!(formatted.contains("combined_protein.tsv"));
}

#[test]
fn no_experiment_subfolders() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::remove_dir_all(bundle.join("control-1")).unwrap();
    fs::remove_dir_all(bundle.join("treatment-1")).unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::NoExperimentFolders]);
    assert!(
        inspection.errors[0]
            .message
            .contains("no experiment subfolder")
    );
}

#[test]
fn missing_fasta_warns_and_parsing_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::remove_file(bundle.join("protein.fas")).unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    assert_eq!(
        codes(&inspection.warnings),
        [IssueCode::FilteredPsms, IssueCode::MissingFasta]
    );
    assert_eq!(inspection.selected_fasta, None);
    assert_eq!(inspection.fasta_record_count, None);
    let batch = import_inspected(&inspection).unwrap();
    assert!(
        batch
            .proteins
            .iter()
            .all(|protein| protein.sequence.is_none())
    );
}

#[test]
fn multiple_fasta_candidates_select_deterministically() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::write(bundle.join("aa_database.fa"), ">P00001\nMS\n").unwrap();
    fs::write(bundle.join("zz_backup.FAS"), ">P00001\nMS\n").unwrap();

    let inspection = inspect_bundle(&bundle);

    // Case-insensitive order: aa_database.fa < protein.fas < zz_backup.FAS.
    assert_eq!(
        inspection.fasta_candidates,
        ["aa_database.fa", "protein.fas", "zz_backup.FAS"]
    );
    assert_eq!(
        inspection.selected_fasta,
        Some(bundle.join("aa_database.fa"))
    );
    // The dummy database holds only P00001, so P00002/P00003 lose sequences.
    assert_eq!(
        codes(&inspection.warnings),
        [
            IssueCode::FilteredPsms,
            IssueCode::AmbiguousFasta,
            IssueCode::UnmatchedFastaProteins
        ]
    );
    assert!(inspection.warnings[1].message.contains("aa_database.fa"));

    // The selected database holds only P00001 (sequence "MS"), so only that
    // protein gets a sequence.
    let batch = import_inspected(&inspection).unwrap();
    let sequences: Vec<(&str, Option<&str>)> = batch
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p.sequence.as_deref()))
        .collect();
    assert_eq!(
        sequences,
        [("P00001", Some("MS")), ("P00002", None), ("P00003", None)]
    );
}

#[test]
fn empty_required_report_warns_and_parsing_succeeds() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::write(
        bundle.join("combined_protein.tsv"),
        "Protein ID\tIndistinguishable Proteins\n",
    )
    .unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    assert_eq!(
        codes(&inspection.warnings),
        [IssueCode::EmptyReport, IssueCode::FilteredPsms]
    );
    assert_eq!(inspection.row_counts["combined_protein.tsv"], 0);
    // combined_protein.tsv is validated, not imported: PSMs still load.
    assert_eq!(import_inspected(&inspection).unwrap().psms.len(), 5);
}

#[test]
fn zero_byte_report_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::write(bundle.join("psm.tsv"), b"").unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::EmptyFile]);
    assert!(inspection.errors[0].message.contains("psm.tsv"));
    assert!(!inspection.row_counts.contains_key("psm.tsv"));
    assert!(
        import_bundle(&bundle)
            .unwrap_err()
            .to_string()
            .contains("psm.tsv")
    );
}

#[test]
fn invalid_encoding_names_file_and_line() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::write(
        bundle.join("combined_peptide.tsv"),
        b"Peptide Sequence\tProtein ID\nPEPTIDEK\tP00001\n\xff\xfe\tP00002\n",
    )
    .unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::InvalidEncoding]);
    let message = &inspection.errors[0].message;
    assert!(message.contains("combined_peptide.tsv"));
    assert!(message.contains("line 3"));
}

#[test]
fn ragged_row_names_line_and_field_counts() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    let path = bundle.join("psm.tsv");
    let text = fs::read_to_string(&path).unwrap();
    let mut lines: Vec<String> = text.lines().map(str::to_string).collect();
    // Drop the trailing (empty) Mapped Proteins field of the first data row.
    let trimmed = lines[1].rsplit_once('\t').unwrap().0.to_string();
    lines[1] = trimmed;
    fs::write(&path, lines.join("\n") + "\n").unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::MalformedTsv]);
    let message = &inspection.errors[0].message;
    assert!(message.contains("psm.tsv"));
    assert!(message.contains("line 2"));
    assert!(message.contains("37 fields") && message.contains("38"));
}

#[test]
fn blank_data_line_is_malformed() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    let path = bundle.join("combined_protein.tsv");
    let text = fs::read_to_string(&path).unwrap();
    fs::write(&path, text.replacen('\n', "\n\n", 1)).unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::MalformedTsv]);
    assert!(inspection.errors[0].message.contains("line 2 is blank"));
}

#[test]
fn psm_tables_are_flagged_filtered_and_lose_post_processor_columns() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    assert_eq!(codes(&inspection.warnings), [IssueCode::FilteredPsms]);
    assert!(inspection.warnings[0].message.contains(".pepXML"));
    let batch = import_inspected(&inspection).unwrap();
    assert_eq!(batch.psms.len(), 5);
    assert!(batch.psms.iter().all(|psm| {
        !psm.scores
            .iter()
            .any(|s| s.name.starts_with("PeptideProphet"))
    }));
}

#[test]
fn missing_required_columns_name_field_variants_and_headers() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    drop_column(&bundle.join("psm.tsv"), "Charge");
    drop_column(&bundle.join("psm.tsv"), "Hyperscore");
    drop_column(&bundle.join("combined_protein.tsv"), "Protein ID");

    let inspection = inspect_bundle(&bundle);

    let column_errors = with_code(&inspection.errors, IssueCode::MissingRequiredColumns);
    assert_eq!(column_errors.len(), 3);
    let charge = column_errors
        .iter()
        .find(|m| m.contains("'charge'"))
        .unwrap();
    assert!(charge.contains("psm.tsv"));
    assert!(charge.contains("'charge', 'assumed charge'"));
    assert!(charge.contains("'Spectrum'"));
    let score = column_errors
        .iter()
        .find(|m| m.contains("'score'"))
        .unwrap();
    assert!(score.contains("'hyperscore', 'score', 'xcorr'"));
    assert!(inspection.format_errors().contains("protein id"));
}

#[test]
fn orphan_subfolder_stays_a_zero_psm_sample() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::create_dir(bundle.join("ghost-1")).unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    let warnings = with_code(&inspection.warnings, IssueCode::UnmatchedSamples);
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("ghost-1"));
    assert!(warnings[0].contains("zero PSMs"));
    let batch = import_inspected(&inspection).unwrap();
    assert!(
        batch
            .biosamples
            .iter()
            .any(|sample| sample.name == "ghost-1")
    );
    assert!(batch.ms_runs.iter().all(|run| run.biosample != "ghost-1"));
}

#[test]
fn psm_samples_without_subfolder_warn_both_ways() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = minimal_bundle(dir.path());
    write_tsv(
        &bundle.join("psm.tsv"),
        &[
            "Spectrum",
            "Spectrum File",
            "Peptide",
            "Charge",
            "Hyperscore",
            "PeptideProphet Probability",
            "Protein ID",
        ],
        &[&[
            "control-1.1.1.2",
            "interact-nowhere-1.raw",
            "PEPTIDEK",
            "2",
            "9.9",
            "0.99",
            "P1",
        ]],
    );
    write_tsv(
        &bundle.join("combined_peptide.tsv"),
        &["Peptide Sequence", "Protein ID"],
        &[&["PEPTIDEK", "P1"]],
    );
    write_tsv(
        &bundle.join("combined_protein.tsv"),
        &["Protein ID"],
        &[&["P1"]],
    );

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    // control-1 has no PSM rows, and the derived id nowhere-1 has no subfolder.
    let warnings = with_code(&inspection.warnings, IssueCode::UnmatchedSamples);
    assert_eq!(warnings.len(), 2);
    assert!(warnings.iter().any(|m| m.contains("control-1")));
    let unknown = warnings.iter().find(|m| m.contains("nowhere-1")).unwrap();
    assert!(unknown.contains("remain loaded"));
}

#[test]
fn unmatched_protein_ids_and_unsupported_peptides() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    append_row(
        &bundle.join("combined_peptide.tsv"),
        &[
            ("Peptide Sequence", "UNSUPPORTEDPEP"),
            ("Protein ID", "P99999"),
        ],
    );

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    let proteins = with_code(&inspection.warnings, IssueCode::UnmatchedProteinIds);
    assert_eq!(proteins.len(), 1);
    assert!(proteins[0].contains("P99999"));
    assert!(proteins[0].contains("Mapped Proteins"));
    assert!(proteins[0].contains("remain loaded"));
    let peptides = with_code(&inspection.warnings, IssueCode::UnmatchedPeptides);
    assert_eq!(peptides.len(), 1);
    assert!(peptides[0].contains("UNSUPPORTEDPEP"));
    assert!(peptides[0].contains("no linked PSMs"));
}

#[test]
fn unmatched_mapped_protein_is_reported_without_its_primary() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    append_row(
        &bundle.join("combined_peptide.tsv"),
        &[
            ("Peptide Sequence", "MAPPEDPEPK"),
            ("Protein ID", "P00001"),
            ("Mapped Proteins", "tr|Q9MAPPED|Q9MAPPED_STREN"),
        ],
    );

    let inspection = inspect_bundle(&bundle);

    let proteins = with_code(&inspection.warnings, IssueCode::UnmatchedProteinIds);
    assert_eq!(proteins.len(), 1);
    assert!(proteins[0].contains("tr|Q9MAPPED|Q9MAPPED_STREN"));
    assert!(!proteins[0].contains("P00001"));
}

#[test]
fn group_members_missing_from_fasta() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::write(
        bundle.join("protein.fas"),
        ">sp|P00001|P1_HUMAN Protein one\nAGSTPV\n",
    )
    .unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    let warnings = with_code(&inspection.warnings, IssueCode::UnmatchedFastaProteins);
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("P00002") && warnings[0].contains("P00003"));
}

#[test]
fn invalid_fasta_encoding_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::write(bundle.join("protein.fas"), b">P00001\n\xff\xfeMS\n").unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::InvalidEncoding]);
    assert!(inspection.errors[0].message.contains("protein.fas"));
}

#[test]
fn invalid_paths() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("not-a-bundle.tsv");
    fs::write(&file, "Spectrum\n").unwrap();

    let inspection = inspect_bundle(&file);
    assert_eq!(codes(&inspection.errors), [IssueCode::InvalidPath]);
    assert!(inspection.errors[0].message.contains("directory"));

    let missing = inspect_bundle(dir.path().join("does-not-exist"));
    assert_eq!(codes(&missing.errors), [IssueCode::InvalidPath]);
    assert!(missing.errors[0].message.contains("does not exist"));
    assert!(missing.discovered_files.is_empty());
    assert!(missing.row_counts.is_empty());
}

#[test]
fn import_inspected_equals_import_bundle() {
    let inspection = inspect_bundle(fixture_bundle());
    assert_eq!(
        import_inspected(&inspection).unwrap(),
        import_bundle(fixture_bundle()).unwrap()
    );
}

#[test]
fn format_errors_lists_every_error() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    fs::remove_file(bundle.join("ion.tsv")).unwrap();
    fs::remove_dir_all(bundle.join("control-1")).unwrap();
    fs::remove_dir_all(bundle.join("treatment-1")).unwrap();

    let formatted = inspect_bundle(&bundle).format_errors();

    assert!(formatted.starts_with(&format!(
        "Invalid FragPipe bundle at '{}':",
        bundle.display()
    )));
    assert!(formatted.contains("[missing-required-file] missing required report 'ion.tsv'"));
    assert!(formatted.contains("[no-experiment-folders]"));
}

#[test]
fn inspection_serializes_for_the_ui() {
    let json = serde_json::to_value(inspect_bundle(fixture_bundle())).unwrap();
    assert_eq!(json["row_counts"]["psm.tsv"], 5);
    assert_eq!(json["warnings"].as_array().unwrap().len(), 1);
    assert_eq!(json["warnings"][0]["code"], "filtered-psms");
    assert!(json.get("resolved").is_none());
}

#[test]
fn psm_peptides_missing_from_the_peptide_report_warn() {
    // Common when PSM- and peptide-level FDR filters differ: the PSM stays
    // loaded but its peptide_id links to no peptide record.
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    let path = bundle.join("combined_peptide.tsv");
    let text = fs::read_to_string(&path).unwrap();
    let kept: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with("FOURTHPEP\t"))
        .collect();
    fs::write(&path, kept.join("\n") + "\n").unwrap();

    let inspection = inspect_bundle(&bundle);

    assert_clean(&inspection);
    let warnings = with_code(&inspection.warnings, IssueCode::UnlistedPsmPeptides);
    assert_eq!(warnings.len(), 1);
    assert!(warnings[0].contains("FOURTHPEP"));
    assert!(warnings[0].contains("remain loaded"));
}

#[test]
fn duplicate_peptide_rows_block_import() {
    // Peptide ids derive from the sequence, so a repeated row would create two
    // records with one id.
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bundle");
    append_row(
        &bundle.join("combined_peptide.tsv"),
        &[("Peptide Sequence", "PEPTIDEK"), ("Protein ID", "P00001")],
    );

    let inspection = inspect_bundle(&bundle);

    assert_eq!(codes(&inspection.errors), [IssueCode::DuplicatePeptides]);
    let message = &inspection.errors[0].message;
    assert!(
        message.contains("PEPTIDEK") && message.contains("lines 2 and 6"),
        "{message}"
    );
}
