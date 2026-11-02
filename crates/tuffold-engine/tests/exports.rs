//! Results exports read back: mzIdentML through Tuffold's own importer, and
//! the spectral library's entry count against the store.

use std::collections::HashMap;

use tuffold_engine::{create_demo_project, write_msp, write_mzid};
use tuffold_store::Project;

#[test]
fn mzidentml_round_trips_through_the_importer() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let project = Project::open(&path).unwrap();
    let file = dir.path().join("demo.mzid");
    let report = write_mzid(&project, &mut std::fs::File::create(&file).unwrap()).unwrap();
    let data = project.exchange().unwrap();
    assert_eq!(report.spectra, data.spectra.len());
    assert!(report.passing_spectra > 0 && report.passing_groups > 0);

    let batch = tuffold_io::import_mzid(&file).unwrap();
    // One rank-1 item per validated spectrum, target and decoy.
    assert_eq!(batch.psms.len(), data.spectra.len());
    assert!(batch.psms.iter().all(|psm| psm.rank == 1));
    let decoys = batch
        .psms
        .iter()
        .filter(|psm| batch.psm_is_decoy(psm))
        .count();
    assert_eq!(decoys, data.spectra.iter().filter(|s| s.is_decoy).count());
    // Passing targets as the report counts them; passing decoys are marked
    // as passing in the file too.
    let passing = |decoy: bool| {
        batch
            .psms
            .iter()
            .filter(|psm| psm.pass_threshold == Some(true) && batch.psm_is_decoy(psm) == decoy)
            .count()
    };
    assert_eq!(passing(false), report.passing_spectra);
    assert_eq!(
        passing(true),
        data.spectra
            .iter()
            .filter(|s| s.passes && s.is_decoy)
            .count()
    );
    // Each spectrum keeps its peptide (modifications included) and q-value.
    let forms: HashMap<i64, String> = data
        .peptides
        .iter()
        .map(|p| (p.id, p.form.proforma()))
        .collect();
    let mut expected: Vec<(String, String, String)> = data
        .spectra
        .iter()
        .map(|s| {
            let run = &data.runs.iter().find(|r| r.id == s.run_id).unwrap().name;
            (run.clone(), s.key.clone(), forms[&s.peptide_id].clone())
        })
        .collect();
    let mut found: Vec<(String, String, String)> = batch
        .psms
        .iter()
        .map(|psm| {
            (
                psm.ms_run.clone(),
                psm.spectrum.canonical(),
                psm.peptide.proforma(),
            )
        })
        .collect();
    expected.sort();
    found.sort();
    assert_eq!(found, expected);
    let q = batch.psms[0]
        .scores
        .iter()
        .find(|s| s.accession.as_deref() == Some("MS:1002354"))
        .unwrap();
    assert!((0.0..=1.0).contains(&q.value));
}

/// What the mzIdentML 1.2.0 XSD rejected in the first exports (checked
/// with `xmllint --schema`; CI has no schema validator): every data file
/// needs a `FileFormat` before its `SpectrumIDFormat`, and `location` is an
/// `xs:anyURI`, so no spaces.
#[test]
fn mzidentml_inputs_follow_the_schema() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let project = Project::open(&path).unwrap();
    let mut xml = Vec::new();
    write_mzid(&project, &mut xml).unwrap();
    let xml = String::from_utf8(xml).unwrap();
    let inputs = &xml[xml.find("<Inputs>").unwrap()..xml.find("</Inputs>").unwrap()];
    let blocks: Vec<&str> = inputs.split("<SpectraData ").skip(1).collect();
    assert_eq!(blocks.len(), 4, "one per demo run");
    for block in &blocks {
        let format = block.find("<FileFormat>").expect("a FileFormat");
        let id_format = block
            .find("<SpectrumIDFormat>")
            .expect("a SpectrumIDFormat");
        assert!(format < id_format, "FileFormat comes first: {block}");
    }
    for location in inputs.split("location=\"").skip(1) {
        let value = &location[..location.find('"').unwrap()];
        assert!(
            !value.is_empty() && !value.contains(' '),
            "{value:?} is no URI"
        );
    }
}

#[test]
fn the_library_has_one_entry_per_passing_peptide_ion_with_peaks() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("demo.tuff");
    create_demo_project(&path).unwrap();
    let project = Project::open(&path).unwrap();
    let mut text = Vec::new();
    let report = write_msp(&project, &mut text).unwrap();
    let library = project.library().unwrap();
    let text = String::from_utf8(text).unwrap();
    assert_eq!(report.entries, text.matches("\nNum peaks: ").count());
    assert_eq!(
        report.entries + report.unnamed_modifications,
        library.entries.len()
    );
    assert!(report.entries > 100, "{}", report.entries);
}
