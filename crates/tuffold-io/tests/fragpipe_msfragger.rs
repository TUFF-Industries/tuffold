//! FragPipe folders that still hold MSFragger's own pepXML files.
//!
//! The folders are built here; every expectation is read off the files
//! written below.

use std::fs;
use std::path::Path;

use tuffold_core::{BioSample, SourceFormat, SpectrumKey};
use tuffold_io::{find_search_files, import_path};

/// One MSFragger pepXML run. MSFragger calls itself `X! Tandem` in
/// `search_engine` and names itself in the version. Each hit is
/// (scan, peptide, protein, hyperscore).
fn msfragger_pepxml(
    engine: &str,
    version: &str,
    run: &str,
    hits: &[(u32, &str, &str, f64)],
) -> String {
    let mut xml = format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<msms_pipeline_analysis>
<msms_run_summary base_name="/data/{run}" raw_data_type="mzML" raw_data=".mzML">
<search_summary base_name="/data/{run}" precursor_mass_type="monoisotopic" search_engine="{engine}" search_engine_version="{version}" fragment_mass_type="monoisotopic" search_id="1">
</search_summary>
"#
    );
    for (scan, peptide, protein, hyperscore) in hits {
        xml.push_str(&format!(
            r#"<spectrum_query spectrum="{run}.{scan}.{scan}.2" start_scan="{scan}" end_scan="{scan}" assumed_charge="2" index="{scan}">
<search_result>
<search_hit hit_rank="1" peptide="{peptide}" peptide_prev_aa="K" peptide_next_aa="A" protein="{protein}" num_tot_proteins="1">
<search_score name="hyperscore" value="{hyperscore}"/>
<search_score name="nextscore" value="10.0"/>
<search_score name="expect" value="0.5"/>
</search_hit>
</search_result>
</spectrum_query>
"#
        ));
    }
    xml.push_str("</msms_run_summary>\n</msms_pipeline_analysis>\n");
    xml
}

const TARGET: &str = "sp|P00001|P1_HUMAN";
const DECOY: &str = "rev_sp|P00001|P1_HUMAN";

/// control-1/A.pepXML: scan 1 PEPTIDEK (target, 30.5), scan 2 KEDITPEP
/// (decoy, 12.0); treatment-1/B.pepXML: scan 7 PEPTIDEK (target, 28.0).
/// Next to them, files that must not be read: PeptideProphet's
/// `interact-A.pep.xml`, MSFragger's `B_rank2.pepXML`, a hidden pepXML
/// (all three not even XML), and a FASTA in control-1.
fn folder(root: &Path) -> std::path::PathBuf {
    let folder = root.join("fragpipe");
    let control = folder.join("control-1");
    let treatment = folder.join("treatment-1");
    fs::create_dir_all(&control).unwrap();
    fs::create_dir_all(&treatment).unwrap();
    let a = msfragger_pepxml(
        "X! Tandem",
        "MSFragger-4.1",
        "A",
        &[(1, "PEPTIDEK", TARGET, 30.5), (2, "KEDITPEP", DECOY, 12.0)],
    );
    let b = msfragger_pepxml(
        "X! Tandem",
        "MSFragger-4.1",
        "B",
        &[(7, "PEPTIDEK", TARGET, 28.0)],
    );
    fs::write(control.join("A.pepXML"), a).unwrap();
    fs::write(treatment.join("B.pepXML"), b).unwrap();
    fs::write(control.join("interact-A.pep.xml"), "not xml").unwrap();
    fs::write(treatment.join("B_rank2.pepXML"), "not xml").unwrap();
    fs::write(folder.join(".A.pepXML"), "not xml").unwrap();
    fs::write(
        control.join("protein.fas"),
        format!(">{TARGET} Protein one\nMKPEPTIDEKA\n"),
    )
    .unwrap();
    fs::write(folder.join("combined_protein.tsv"), "Protein ID\n").unwrap();
    folder
}

#[test]
fn only_the_engines_own_files_are_found() {
    let dir = tempfile::tempdir().unwrap();
    let folder = folder(dir.path());
    let found = find_search_files(&folder).unwrap();
    let names: Vec<(String, Option<&str>)> = found
        .iter()
        .map(|file| {
            (
                file.path
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .into_owned(),
                file.experiment.as_deref(),
            )
        })
        .collect();
    assert_eq!(
        names,
        [
            ("A.pepXML".to_string(), Some("control-1")),
            ("B.pepXML".to_string(), Some("treatment-1"))
        ]
    );
}

#[test]
fn msfragger_results_are_imported_instead_of_psm_tsv() {
    let dir = tempfile::tempdir().unwrap();
    let imported = import_path(&folder(dir.path())).unwrap();

    assert_eq!(imported.format, SourceFormat::FragPipeBundle);
    assert_eq!(imported.batches.len(), 1);
    assert_eq!(imported.notes.len(), 1);
    assert!(
        imported.notes[0].starts_with("read MSFragger's own results (2 pepXML files"),
        "{}",
        imported.notes[0]
    );
    let batch = &imported.batches[0];
    assert_eq!(batch.software.name, "MSFragger");
    assert_eq!(batch.software.version.as_deref(), Some("4.1"));
    // Experiments are biosamples, their category the name's prefix.
    assert_eq!(
        batch.biosamples,
        [
            BioSample {
                name: "control-1".into(),
                category: Some("control".into())
            },
            BioSample {
                name: "treatment-1".into(),
                category: Some("treatment".into())
            }
        ]
    );
    let runs: Vec<(&str, &str, Option<&str>)> = batch
        .ms_runs
        .iter()
        .map(|run| {
            (
                run.name.as_str(),
                run.biosample.as_str(),
                run.spectra_file.as_deref(),
            )
        })
        .collect();
    assert_eq!(
        runs,
        [
            ("A", "control-1", Some("A.mzML")),
            ("B", "treatment-1", Some("B.mzML"))
        ]
    );
    // Every match, the decoy included, with MSFragger's scores only.
    let psms: Vec<(&str, &SpectrumKey, &str)> = batch
        .psms
        .iter()
        .map(|psm| {
            (
                psm.ms_run.as_str(),
                &psm.spectrum,
                psm.peptide.sequence.as_str(),
            )
        })
        .collect();
    assert_eq!(
        psms,
        [
            ("A", &SpectrumKey::Scan(1), "PEPTIDEK"),
            ("A", &SpectrumKey::Scan(2), "KEDITPEP"),
            ("B", &SpectrumKey::Scan(7), "PEPTIDEK")
        ]
    );
    let names: Vec<&str> = batch.psms[0]
        .scores
        .iter()
        .map(|s| s.name.as_str())
        .collect();
    assert_eq!(names, ["hyperscore", "nextscore", "expect"]);
    assert_eq!(batch.psms[2].scores[0].value, 28.0);
    assert!(batch.psm_is_decoy(&batch.psms[1]));
    assert!(!batch.psm_is_decoy(&batch.psms[0]));
    // The target once (both runs name it), its sequence from control-1's
    // FASTA; the decoy is in no FASTA.
    let proteins: Vec<(&str, Option<&str>, bool)> = batch
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p.sequence.as_deref(), p.is_decoy))
        .collect();
    assert_eq!(
        proteins,
        [(TARGET, Some("MKPEPTIDEKA"), false), (DECOY, None, true)]
    );
}

#[test]
fn a_run_in_two_files_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let folder = folder(dir.path());
    fs::copy(
        folder.join("control-1/A.pepXML"),
        folder.join("treatment-1/A-copy.pepXML"),
    )
    .unwrap();
    let error = import_path(&folder).unwrap_err().to_string();
    assert!(error.contains("run A is also in"), "{error}");
}

#[test]
fn searches_by_different_engines_are_refused() {
    let dir = tempfile::tempdir().unwrap();
    let folder = folder(dir.path());
    let comet = msfragger_pepxml(
        "Comet",
        "2015.02 rev. 0",
        "C",
        &[(3, "PEPTIDEK", TARGET, 2.5)],
    );
    fs::write(folder.join("treatment-1/C.pepXML"), comet).unwrap();
    let error = import_path(&folder).unwrap_err().to_string();
    assert!(error.contains("is a Comet search"), "{error}");
    assert!(error.contains("MSFragger searches"), "{error}");
}
