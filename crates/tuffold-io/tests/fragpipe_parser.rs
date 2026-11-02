//! FragPipe importer tests.
//!
//! The oracle expectations are derived BY HAND from
//! fixtures/fragpipe_synthetic/psm.tsv and protein.fas by reading the files
//! and tracing columns. Never regenerate them by running the importer;
//! re-derive them by hand when the fixture or mapping conventions change.

mod common;

use std::fs;
use std::path::{Path, PathBuf};

use common::*;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, MsRun, Psm, Score, Spectrum, SpectrumKey,
};
use tuffold_io::{ParseError, import_bundle};

fn psm_at<'a>(batch: &'a ImportBatch, run: &str, scan: u32) -> &'a Psm {
    batch
        .psms
        .iter()
        .find(|psm| psm.ms_run == run && psm.spectrum == SpectrumKey::Scan(scan))
        .unwrap()
}

fn score(name: &str, value: f64) -> Score {
    Score {
        name: name.into(),
        accession: None,
        value,
    }
}

#[test]
fn matches_hand_derived_oracle() {
    // psm.tsv, 5 rows. Traced columns: Spectrum | Spectrum File | Peptide |
    // Charge | Retention | Observed M/Z | Calculated M/Z | Expectation |
    // Hyperscore | Nextscore | PeptideProphet Probability | PeptideProphet
    // Q-value | Protein Start | Protein End | Prev AA | Next AA | Protein ID |
    // Mapped Proteins (Assigned Modifications is empty on every row):
    //   control-1.1001.1001.2   control-1   PEPTIDEK   2 600.1234 473.2528 472.7492 0.0001 42.1 38.2 0.999 0.001 15 22 K A P00001 -
    //   control-1.1002.1002.3   control-1   SEQUENCEK  3 800.5678 334.5241 334.1891 0.0003 38.5 33.1 0.996 0.004 40 48 K A P00001 P00002
    //   control-1.1003.1003.2   control-1   FOURTHPEP  2 700.2345 546.2921 545.7893 0.0009 29.4 25.6 0.988 0.012 77 85 K A P00002 -
    //   treatment-1.2001.2001.2 treatment-1 PEPTIDEK   2 610.2345 473.2528 472.7492 0.0001 45.2 41.3 0.998 0.002 15 22 K A P00001 -
    //   treatment-1.2002.2002.2 treatment-1 ANOTHERPEP 2 900.3456 580.3278 579.8254 0.0012 31.7 28.9 0.991 0.009 30 39 K A P00003 -
    // => runs control-1 and treatment-1 (the Spectrum prefix), each in the
    //    biosample of the same name (Spectrum File); 5 spectra keyed by start
    //    scan; 4 distinct peptide forms; proteins in first-seen order P00001,
    //    P00002 (row 2, mapped), P00003.
    let batch = import_bundle(fixture_bundle()).unwrap();

    assert_eq!(batch.software.name, "FragPipe");
    let biosample = |name: &str, category: &str| BioSample {
        name: name.into(),
        category: Some(category.into()),
    };
    assert_eq!(
        batch.biosamples,
        [
            biosample("control-1", "control"),
            biosample("treatment-1", "treatment")
        ]
    );
    let run = |name: &str| MsRun {
        name: name.into(),
        biosample: name.into(),
        spectra_file: None,
    };
    assert_eq!(batch.ms_runs, [run("control-1"), run("treatment-1")]);
    assert_eq!(batch.spectra.len(), 5);
    assert_eq!(
        batch.spectra[0],
        Spectrum {
            ms_run: "control-1".into(),
            key: SpectrumKey::Scan(1001),
            source_id: Some("control-1.1001.1001.2".into()),
            precursor_mz: Some(473.2528),
            precursor_charge: Some(2),
            rt_seconds: Some(600.1234),
        }
    );

    // Full check of row 2: the engine's score columns under their verbatim
    // headers, in import order (primary score, Nextscore, Expectation);
    // PeptideProphet's probability and q-value are not imported.
    let row2 = psm_at(&batch, "control-1", 1002);
    assert_eq!(row2.rank, 1);
    assert_eq!(row2.peptide.proforma(), "SEQUENCEK");
    assert_eq!(row2.charge, 3);
    assert_eq!(row2.calculated_mz, Some(334.1891));
    assert_eq!(row2.experimental_mz, Some(334.5241));
    assert_eq!(row2.pass_threshold, None);
    assert_eq!(
        row2.scores,
        [
            score("Hyperscore", 38.5),
            score("Nextscore", 33.1),
            score("Expectation", 0.0003),
        ]
    );
    assert_eq!(
        row2.evidence,
        [
            Evidence {
                accession: "P00001".into(),
                start: Some(40),
                end: Some(48),
                pre: Some('K'),
                post: Some('A'),
            },
            Evidence {
                accession: "P00002".into(),
                start: None,
                end: None,
                pre: None,
                post: None,
            },
        ]
    );

    let mut peptides: Vec<String> = batch
        .psms
        .iter()
        .map(|psm| psm.peptide.proforma())
        .collect();
    peptides.sort();
    peptides.dedup();
    assert_eq!(
        peptides,
        ["ANOTHERPEP", "FOURTHPEP", "PEPTIDEK", "SEQUENCEK"]
    );

    // Proteins: descriptions from the primary row's Protein Description;
    // sequences from protein.fas, where the decoy rev_sp|P00001 is listed
    // first and must lose "P00001" to the forward entry. P00003 is stored
    // lower-case in the FASTA and is folded. P00004 is referenced nowhere.
    let proteins: Vec<(&str, Option<&str>, bool)> = batch
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p.description.as_deref(), p.is_decoy))
        .collect();
    assert_eq!(
        proteins,
        [
            ("P00001", Some("Protein one OS=Homo sapiens OX=9606"), false),
            ("P00002", Some("Protein two OS=Homo sapiens OX=9606"), false),
            (
                "P00003",
                Some("Protein three OS=Homo sapiens OX=9606"),
                false
            ),
        ]
    );
    let p1 = format!(
        "MAGTASTPGVQAGAPEPTIDEKGVSTAGSPGVQTAGVQASEQUENCEK{}AG",
        "AGSTPV".repeat(50)
    );
    let p3 = format!(
        "{}AGSTPANOTHERPEP{}A",
        "AGSTPV".repeat(4),
        "AGSTPV".repeat(40)
    );
    assert_eq!(batch.proteins[0].sequence.as_deref(), Some(p1.as_str()));
    assert_eq!(batch.proteins[2].sequence.as_deref(), Some(p3.as_str()));
    assert!(batch.psms.iter().all(|psm| !batch.psm_is_decoy(psm)));
}

#[test]
fn assigned_modifications_become_peptide_forms() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        &[
            "Spectrum",
            "Spectrum File",
            "Peptide",
            "Charge",
            "Hyperscore",
            "Protein ID",
            "Assigned Modifications",
        ],
        &[&[
            "r1.7.7.2",
            "control-1",
            "MPEPTIDEK",
            "2",
            "9.9",
            "P1",
            "1M(15.9949), N-term(42.0106)",
        ]],
    );

    let batch = import_bundle(&bundle).unwrap();

    let form = &batch.psms[0].peptide;
    assert_eq!(form.proforma(), "[+42.0106]-M[+15.9949]PEPTIDEK");
    assert_eq!(form.modifications[1].site, ModSite::Residue(1));
}

#[test]
fn malformed_modifications_name_their_line() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        &[
            "Spectrum",
            "Spectrum File",
            "Peptide",
            "Charge",
            "Hyperscore",
            "Protein ID",
            "Assigned Modifications",
        ],
        &[&[
            "r1.7.7.2",
            "control-1",
            "MPEPTIDEK",
            "2",
            "9.9",
            "P1",
            "M(oops)",
        ]],
    );
    let error = import_bundle(&bundle).unwrap_err();
    assert!(
        matches!(
            error,
            ParseError::InvalidValue {
                line: 2,
                column: "assigned_modifications",
                ..
            }
        ),
        "{error}"
    );
}

#[test]
fn fractions_of_one_experiment_are_separate_runs_in_one_biosample() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        PSM_HEADER,
        &[
            &[
                "liver_F1.10.10.2",
                "interact-control-1.pep.xml",
                "PEPTIDEK",
                "2",
                "9.9",
                "0.01",
                "P1",
            ],
            &[
                "liver_F2.10.10.2",
                "interact-control-1.pep.xml",
                "PEPTIDEK",
                "2",
                "9.9",
                "0.01",
                "P1",
            ],
        ],
    );

    let batch = import_bundle(&bundle).unwrap();

    let runs: Vec<(&str, &str)> = batch
        .ms_runs
        .iter()
        .map(|run| (run.name.as_str(), run.biosample.as_str()))
        .collect();
    assert_eq!(runs, [("liver_F1", "control-1"), ("liver_F2", "control-1")]);
    assert_eq!(batch.biosamples.len(), 1);
}

#[test]
fn missing_required_columns_are_reported_together() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "broken");
    drop_column(&bundle.join("psm.tsv"), "Charge");
    drop_column(&bundle.join("psm.tsv"), "Hyperscore");
    drop_column(&bundle.join("combined_protein.tsv"), "Protein ID");

    let message = import_bundle(&bundle).unwrap_err().to_string();
    for needle in [
        "psm.tsv",
        "combined_protein.tsv",
        "'charge'",
        "'hyperscore'",
        "'protein id'",
    ] {
        assert!(message.contains(needle), "{needle} missing from: {message}");
    }
}

#[test]
fn variant_headers_are_normalized() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        &[
            "Spectrum",
            "Sample",
            "Peptide",
            "Assumed Charge",
            "Score",
            "PeptideProphet Probability",
            "Q_Value",
            "Protein ID",
            "Mapped Proteins",
        ],
        &[&[
            "control-1.5.5.2",
            "interact-control-1.raw",
            "PEPTIDEK",
            "2",
            "9.9",
            "0.99",
            "0.01",
            "P1",
            "P2",
        ]],
    );

    let batch = import_bundle(&bundle).unwrap();

    let psm = &batch.psms[0];
    assert_eq!(psm.ms_run, "control-1");
    assert_eq!(psm.spectrum, SpectrumKey::Scan(5));
    assert_eq!(psm.charge, 2);
    // The probability and Q_Value columns are PeptideProphet's: not imported.
    assert_eq!(psm.scores, [score("Score", 9.9)]);
    let accessions: Vec<&str> = psm.evidence.iter().map(|e| e.accession.as_str()).collect();
    assert_eq!(accessions, ["P1", "P2"]);
}

#[test]
fn non_numeric_values_name_file_line_and_column() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        PSM_HEADER,
        &[
            &["s.1.1.2", "control-1", "PEPTIDEK", "2", "9.9", "0.01", "P1"],
            &[
                "s.2.2.2",
                "control-1",
                "PEPTIDEK",
                "two",
                "9.9",
                "0.01",
                "P1",
            ],
        ],
    );
    let error = import_bundle(&bundle).unwrap_err();
    assert!(matches!(
        error,
        ParseError::InvalidValue {
            line: 3,
            column: "charge",
            ..
        }
    ));
    assert!(error.to_string().contains("psm.tsv: line 3"));
}

#[test]
fn blank_key_cells_are_rejected_with_their_line() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        PSM_HEADER,
        &[
            &["s.1.1.2", "control-1", "PEPTIDEK", "2", "9.9", "0.01", "P1"],
            &["", "control-1", "PEPTIDEK", "2", "9.9", "0.01", "P1"],
        ],
    );
    let error = import_bundle(&bundle).unwrap_err();
    assert!(
        matches!(
            error,
            ParseError::InvalidValue {
                file: "psm.tsv",
                line: 3,
                column: "psm_id",
                ..
            }
        ),
        "{error}"
    );
}

#[test]
fn non_finite_numbers_are_rejected_not_recorded() {
    // NaN would otherwise count as a recorded score while meaning unknown.
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        PSM_HEADER,
        &[&["s.1.1.2", "control-1", "PEPTIDEK", "2", "NaN", "0.01", "P1"]],
    );
    let error = import_bundle(&bundle).unwrap_err();
    assert!(
        matches!(
            error,
            ParseError::InvalidValue {
                column: "score",
                line: 2,
                ..
            }
        ),
        "{error}"
    );
}

#[test]
fn post_processor_q_values_are_not_imported() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        PSM_HEADER,
        &[&["s.1.1.2", "control-1", "PEPTIDEK", "2", "9.9", "0.01", "P1"]],
    );
    let batch = import_bundle(&bundle).unwrap();
    assert_eq!(batch.psms[0].scores, [score("Hyperscore", 9.9)]);
}

#[test]
fn missing_files_and_folders_are_reported_together() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = dir.path().join("incomplete");
    fs::create_dir(&bundle).unwrap();
    fs::write(bundle.join("psm.tsv"), "Spectrum\tSpectrum File\n").unwrap();

    let message = import_bundle(&bundle).unwrap_err().to_string();
    assert!(message.contains("combined_protein.tsv"));
    assert!(message.contains("ion.tsv"));
    assert!(message.contains("no experiment subfolder"));
}

#[test]
fn requires_a_directory() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("not-a-bundle.tsv");
    fs::write(&file, "Spectrum\n").unwrap();
    assert!(
        import_bundle(&file)
            .unwrap_err()
            .to_string()
            .contains("directory")
    );
}

#[test]
fn utf8_byte_order_marks_are_tolerated() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "bom");
    for name in ["psm.tsv", "combined_peptide.tsv", "combined_protein.tsv"] {
        let path = bundle.join(name);
        let mut bytes = b"\xef\xbb\xbf".to_vec();
        bytes.extend(fs::read(&path).unwrap());
        fs::write(&path, bytes).unwrap();
    }
    let batch = import_bundle(&bundle).unwrap();
    assert_eq!(batch.psms.len(), 5);
    assert_eq!(
        batch.spectra[0].source_id.as_deref(),
        Some("control-1.1001.1001.2")
    );
}

#[test]
fn bundle_without_fasta_has_no_sequences() {
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "no-fasta");
    fs::remove_file(bundle.join("protein.fas")).unwrap();
    let batch = import_bundle(&bundle).unwrap();
    assert!(
        batch
            .proteins
            .iter()
            .all(|protein| protein.sequence.is_none())
    );
    assert_eq!(batch.proteins.len(), 3);
}

#[test]
fn hidden_entries_are_not_part_of_a_bundle() {
    // .DS_Store, macOS ._* resource forks, and dot-directories (e.g.
    // Philosopher's .meta/) are ignored instead of becoming files, FASTA
    // candidates, or biosamples.
    let dir = tempfile::tempdir().unwrap();
    let bundle = copy_bundle(dir.path(), "hidden");
    fs::create_dir(bundle.join(".meta")).unwrap();
    fs::write(bundle.join(".DS_Store"), b"\x00\x01").unwrap();
    fs::write(bundle.join("._protein.fas"), b"\x00\x05\x16\x07\xff").unwrap();

    let batch = import_bundle(&bundle).unwrap();
    let samples: Vec<&str> = batch.biosamples.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(samples, ["control-1", "treatment-1"]);
    assert!(batch.proteins[0].sequence.is_some());
}

const PSM_HEADER: &[&str] = &[
    "Spectrum",
    "Spectrum File",
    "Peptide",
    "Charge",
    "Hyperscore",
    "Q-Value",
    "Protein ID",
];

/// A minimal valid bundle whose psm.tsv rows are supplied by the test.
fn bundle_with_psm_rows(dir: &Path, header: &[&str], rows: &[&[&str]]) -> PathBuf {
    let bundle = minimal_bundle(dir);
    write_tsv(&bundle.join("psm.tsv"), header, rows);
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
    bundle
}

#[test]
fn modification_positions_outside_the_peptide_are_rejected() {
    // PEPTIDEK has 8 residues; position 9 would otherwise vanish from the
    // peptide key and merge this PSM into the unmodified form.
    let dir = tempfile::tempdir().unwrap();
    let bundle = bundle_with_psm_rows(
        dir.path(),
        &[
            "Spectrum",
            "Spectrum File",
            "Peptide",
            "Charge",
            "Hyperscore",
            "Protein ID",
            "Assigned Modifications",
        ],
        &[&[
            "r1.7.7.2",
            "control-1",
            "PEPTIDEK",
            "2",
            "9.9",
            "P1",
            "9K(42.0106)",
        ]],
    );
    let error = import_bundle(&bundle).unwrap_err();
    assert!(
        matches!(
            error,
            ParseError::InvalidValue {
                line: 2,
                column: "assigned_modifications",
                ..
            }
        ),
        "{error}"
    );
}
