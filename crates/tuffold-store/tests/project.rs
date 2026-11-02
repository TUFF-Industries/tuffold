//! `.tuff` project file behavior. Batches are built inline so every expected
//! count can be read straight off the test.

use std::path::Path;

use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, SourceFormat, Spectrum, SpectrumKey,
};
use tuffold_store::{ImportSource, Project, StoreError};

fn spectrum(run: &str, scan: u32, mz: Option<f64>) -> Spectrum {
    Spectrum {
        ms_run: run.into(),
        key: SpectrumKey::Scan(scan),
        source_id: None,
        precursor_mz: mz,
        precursor_charge: Some(2),
        rt_seconds: None,
    }
}

fn protein(accession: &str, sequence: Option<&str>) -> Protein {
    Protein {
        accession: accession.into(),
        description: None,
        sequence: sequence.map(str::to_string),
        is_decoy: tuffold_core::decoy::is_decoy_accession(accession),
        decoy_explicit: false,
    }
}

fn psm(run: &str, scan: u32, peptide: PeptideForm, accessions: &[&str], score: (&str, f64)) -> Psm {
    Psm {
        ms_run: run.into(),
        spectrum: SpectrumKey::Scan(scan),
        rank: 1,
        peptide,
        charge: 2,
        calculated_mz: None,
        experimental_mz: None,
        pass_threshold: None,
        scores: vec![Score {
            name: score.0.into(),
            accession: None,
            value: score.1,
        }],
        evidence: accessions
            .iter()
            .map(|accession| Evidence {
                accession: (*accession).into(),
                start: None,
                end: None,
                pre: None,
                post: None,
            })
            .collect(),
    }
}

fn oxidized(sequence: &str, position: u32) -> PeptideForm {
    PeptideForm::new(
        sequence,
        vec![Modification {
            site: ModSite::Residue(position),
            mass_delta: 15.994_915,
            name: Some("Oxidation".into()),
            accession: Some("UNIMOD:35".into()),
        }],
    )
    .unwrap()
}

/// Engine A: run r1, scans 1-3. Scan 3 is a decoy-only hit.
fn engine_a() -> ImportBatch {
    ImportBatch {
        software: Software {
            name: "EngineA".into(),
            version: Some("1.0".into()),
        },
        biosamples: vec![BioSample {
            name: "liver".into(),
            category: Some("control".into()),
        }],
        ms_runs: vec![MsRun {
            name: "r1".into(),
            biosample: "liver".into(),
            spectra_file: None,
        }],
        spectra: vec![
            spectrum("r1", 1, None),
            spectrum("r1", 2, Some(500.25)),
            spectrum("r1", 3, None),
        ],
        proteins: vec![protein("P1", Some("MPEPTIDEK")), protein("rev_P1", None)],
        psms: vec![
            psm(
                "r1",
                1,
                PeptideForm::unmodified("PEPTIDEK"),
                &["P1"],
                ("score", 10.0),
            ),
            psm("r1", 2, oxidized("MPEPTIDEK", 1), &["P1"], ("score", 9.0)),
            psm(
                "r1",
                3,
                PeptideForm::unmodified("KEDITPEP"),
                &["rev_P1"],
                ("score", 2.0),
            ),
        ],
    }
}

/// Engine B searched the same run: scan 1 again (same peptide), scan 2 with a
/// precursor m/z engine A lacked for scan 1, and a new scan 4.
fn engine_b() -> ImportBatch {
    ImportBatch {
        software: Software {
            name: "EngineB".into(),
            version: None,
        },
        biosamples: vec![],
        ms_runs: vec![MsRun {
            name: "r1".into(),
            biosample: "other-name-ignored".into(),
            spectra_file: Some("r1.mzML".into()),
        }],
        spectra: vec![
            spectrum("r1", 1, Some(472.75)),
            spectrum("r1", 2, Some(999.0)),
            spectrum("r1", 4, None),
        ],
        proteins: vec![protein("P1", Some("SHOULDNOTREPLACE")), protein("P2", None)],
        psms: vec![
            psm(
                "r1",
                1,
                PeptideForm::unmodified("PEPTIDEK"),
                &["P1"],
                ("xcorr", 3.1),
            ),
            psm(
                "r1",
                2,
                oxidized("MPEPTIDEK", 1),
                &["P1", "P2"],
                ("xcorr", 2.4),
            ),
            psm(
                "r1",
                4,
                PeptideForm::unmodified("NEWPEPK"),
                &["P2"],
                ("xcorr", 1.9),
            ),
        ],
    }
}

fn source(name: &str, format: SourceFormat) -> ImportSource {
    ImportSource {
        path: format!("/data/{name}"),
        format,
    }
}

fn query<T: rusqlite::types::FromSql>(path: &Path, sql: &str) -> T {
    rusqlite::Connection::open(path)
        .unwrap()
        .query_row(sql, [], |row| row.get(0))
        .unwrap()
}

#[test]
fn create_import_and_summarize() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();

    let report = project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    assert_eq!(
        (
            report.psms,
            report.new_spectra,
            report.new_peptides,
            report.new_proteins,
            report.new_ms_runs
        ),
        (3, 3, 3, 2, 1)
    );
    let summary = project.summary().unwrap();
    assert_eq!(summary.psms, 3);
    assert_eq!(summary.decoy_psms, 1);
    assert_eq!(summary.decoy_proteins, 1);
    assert_eq!(summary.proteins_with_sequence, 1);
    assert_eq!(summary.runs[0].biosample, "liver");
}

#[test]
fn second_engine_shares_spectra_peptides_and_proteins() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    let report = project
        .import(&source("b.pep.xml", SourceFormat::PepXml), &engine_b())
        .unwrap();

    // Only scan 4, NEWPEPK, and P2 are new.
    assert_eq!(
        (
            report.psms,
            report.new_spectra,
            report.new_peptides,
            report.new_proteins,
            report.new_ms_runs
        ),
        (3, 1, 1, 1, 0)
    );
    let summary = project.summary().unwrap();
    assert_eq!(
        (
            summary.spectra,
            summary.psms,
            summary.peptides,
            summary.proteins
        ),
        (4, 6, 4, 3)
    );
    assert_eq!(summary.runs.len(), 1);
    // The run keeps its first biosample; the later name is not applied.
    assert_eq!(summary.runs[0].biosample, "liver");
    assert_eq!(summary.runs[0].psms, 6);
    let records = project.imports().unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(
        (records[1].software.as_str(), records[1].psms),
        ("EngineB", 3)
    );
}

#[test]
fn later_imports_fill_unknowns_but_never_overwrite() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();
    project
        .import(&source("b.pep.xml", SourceFormat::PepXml), &engine_b())
        .unwrap();
    drop(project);

    // Scan 1 had no m/z from A, so B's value fills it; scan 2 keeps A's value.
    let scan1: f64 = query(
        &path,
        "SELECT precursor_mz FROM spectra WHERE spectrum_key = 'scan=1'",
    );
    let scan2: f64 = query(
        &path,
        "SELECT precursor_mz FROM spectra WHERE spectrum_key = 'scan=2'",
    );
    assert_eq!((scan1, scan2), (472.75, 500.25));
    let sequence: String = query(
        &path,
        "SELECT sequence FROM proteins WHERE accession = 'P1'",
    );
    assert_eq!(sequence, "MPEPTIDEK");
    let spectra_file: String = query(&path, "SELECT spectra_file FROM ms_runs WHERE name = 'r1'");
    assert_eq!(spectra_file, "r1.mzML");
}

#[test]
fn modifications_scores_and_evidence_are_stored() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();
    project
        .import(&source("b.pep.xml", SourceFormat::PepXml), &engine_b())
        .unwrap();
    drop(project);

    let proforma: String = query(
        &path,
        "SELECT proforma FROM peptides WHERE sequence = 'MPEPTIDEK'",
    );
    assert_eq!(proforma, "M[+15.9949]PEPTIDEK");
    let accession: String = query(&path, "SELECT accession FROM modifications");
    assert_eq!(accession, "UNIMOD:35");
    let xcorr_rows: i64 = query(
        &path,
        "SELECT COUNT(*) FROM psm_scores WHERE name = 'xcorr'",
    );
    assert_eq!(xcorr_rows, 3);
    // Engine B's scan-2 PSM maps to both P1 and P2.
    let shared: i64 = query(
        &path,
        "SELECT COUNT(*) FROM psm_proteins pp JOIN psms p ON p.id = pp.psm_id
         JOIN imports i ON i.id = p.import_id
         JOIN spectra s ON s.id = p.spectrum_id
         WHERE i.software = 'EngineB' AND s.spectrum_key = 'scan=2'",
    );
    assert_eq!(shared, 2);
}

#[test]
fn reimporting_the_same_file_is_refused_and_changes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    let a = source("a.mzid", SourceFormat::MzIdentMl);
    project.import(&a, &engine_a()).unwrap();

    let error = project.import(&a, &engine_a()).unwrap_err();

    assert!(
        matches!(error, StoreError::AlreadyImported { .. }),
        "{error}"
    );
    assert_eq!(project.summary().unwrap().psms, 3);
}

#[test]
fn one_source_with_several_engines_imports_each_as_its_own_engine() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    let scaffold = source("old.sf3", SourceFormat::Sf3);

    let reports = project
        .import_all(&scaffold, &[engine_a(), engine_b()])
        .unwrap();

    assert_eq!(reports.len(), 2);
    let paths: String = query(
        &path,
        "SELECT group_concat(source_path || '=' || software, ', ') FROM
             (SELECT source_path, software FROM imports ORDER BY id)",
    );
    assert_eq!(
        paths,
        "/data/old.sf3#EngineA=EngineA, /data/old.sf3#EngineB=EngineB"
    );
    // The same file again is refused, whichever way it is imported.
    for again in [
        project.import_all(&scaffold, &[engine_a(), engine_b()]),
        project
            .import(&scaffold, &engine_a())
            .map(|report| vec![report]),
    ] {
        assert!(matches!(again, Err(StoreError::AlreadyImported { .. })));
    }
    let imports: i64 = query(&path, "SELECT COUNT(*) FROM imports");
    assert_eq!(imports, 2);
}

#[test]
fn reopening_preserves_everything() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    Project::create(&path)
        .unwrap()
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    let summary = Project::open(&path).unwrap().summary().unwrap();

    assert_eq!((summary.imports, summary.psms, summary.peptides), (1, 3, 3));
}

#[test]
fn create_refuses_to_overwrite() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    std::fs::write(&path, b"precious").unwrap();
    assert!(matches!(
        Project::create(&path),
        Err(StoreError::AlreadyExists(_))
    ));
    assert_eq!(std::fs::read(&path).unwrap(), b"precious");
}

#[test]
fn open_rejects_files_that_are_not_projects() {
    let dir = tempfile::tempdir().unwrap();
    let text = dir.path().join("notes.tuff");
    std::fs::write(&text, b"not sqlite at all").unwrap();
    assert!(matches!(
        Project::open(&text),
        Err(StoreError::NotAProject(_))
    ));

    let other_sqlite = dir.path().join("other.tuff");
    rusqlite::Connection::open(&other_sqlite)
        .unwrap()
        .execute_batch("CREATE TABLE t (x)")
        .unwrap();
    assert!(matches!(
        Project::open(&other_sqlite),
        Err(StoreError::NotAProject(_))
    ));

    assert!(matches!(
        Project::open(&dir.path().join("missing.tuff")),
        Err(StoreError::NotAProject(_))
    ));
}

#[test]
fn open_refuses_a_newer_schema() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("future.tuff");
    drop(Project::create(&path).unwrap());
    rusqlite::Connection::open(&path)
        .unwrap()
        .pragma_update(None, "user_version", 99)
        .unwrap();
    assert!(matches!(
        Project::open(&path),
        Err(StoreError::UnsupportedVersion(99))
    ));
}

#[test]
fn psm_for_an_undeclared_spectrum_is_kept() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    let batch = ImportBatch {
        spectra: vec![],
        ms_runs: vec![],
        ..engine_a()
    };

    let report = project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &batch)
        .unwrap();

    // The undeclared run r1 is created (and counted) with its 3 spectra.
    assert_eq!(
        (report.psms, report.new_spectra, report.new_ms_runs),
        (3, 3, 1)
    );
    assert_eq!(project.summary().unwrap().runs[0].biosample, "r1");
}

fn peak_spectrum(scan: u32, rt: Option<f64>) -> tuffold_core::PeakSpectrum {
    tuffold_core::PeakSpectrum {
        key: SpectrumKey::Scan(scan),
        native_id: format!("scan={scan}"),
        ms_level: 2,
        rt_seconds: rt,
        precursor_mz: Some(123.4),
        precursor_charge: Some(3),
        mz: vec![100.5, 200.25],
        intensity: vec![10.0, 20.0],
        profile: false,
        parent: None,
        isolation: None,
        precursor_mobility: None,
        mobility: Vec::new(),
    }
}

#[test]
fn spectra_attach_only_to_identified_spectra_of_their_run() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    // Engine A identified scans 1-3 of r1. The file also has scan 99.
    let mut attach = project.attach_spectra("r1").unwrap();
    for scan in [1, 2, 99] {
        attach.add(&peak_spectrum(scan, Some(61.5))).unwrap();
    }
    let report = attach.finish().unwrap();

    assert_eq!(
        (
            report.read,
            report.attached,
            report.identified_without_peaks
        ),
        (3, 2, 1)
    );
    assert_eq!(
        project.peaks("r1", &SpectrumKey::Scan(1)).unwrap(),
        Some(tuffold_store::Peaks {
            mz: vec![100.5, 200.25],
            intensity: vec![10.0, 20.0]
        })
    );
    assert_eq!(project.peaks("r1", &SpectrumKey::Scan(3)).unwrap(), None);
    drop(project);
    // Scan 2 keeps engine A's precursor m/z; its unknown RT is filled.
    let (mz, rt): (f64, f64) = rusqlite::Connection::open(&path)
        .unwrap()
        .query_row(
            "SELECT precursor_mz, rt_seconds FROM spectra WHERE spectrum_key = 'scan=2'",
            [],
            |row| Ok((row.get(0)?, row.get(1)?)),
        )
        .unwrap();
    assert_eq!((mz, rt), (500.25, 61.5));
}

/// An MS3 scan acquired from `parent` (SPS-MS3).
fn ms3_spectrum(scan: u32, parent: u32) -> tuffold_core::PeakSpectrum {
    tuffold_core::PeakSpectrum {
        ms_level: 3,
        mz: vec![126.1277, 127.1248],
        intensity: vec![5.0, 7.0],
        parent: Some(format!("scan={parent}")),
        ..peak_spectrum(scan, Some(62.0))
    }
}

#[test]
fn ms3_scans_are_stored_with_the_spectrum_they_name() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    // Engine A identified scans 1-3 of r1. Scan 101 is scan 1's MS3, 102
    // names an unidentified scan, 103 is a second MS3 for scan 1.
    let mut attach = project.attach_spectra("r1").unwrap();
    attach.add(&peak_spectrum(1, Some(61.5))).unwrap();
    for (scan, parent) in [(101, 1), (102, 99), (103, 1)] {
        attach.add(&ms3_spectrum(scan, parent)).unwrap();
    }
    let report = attach.finish().unwrap();
    assert_eq!(
        (
            report.read,
            report.attached,
            report.ms3_scans,
            report.ms3_attached
        ),
        (1, 1, 3, 1)
    );

    // Attached again with scan 2's MS3: it is stored though scan 2 has no
    // peaks, and scan 1 keeps its first.
    let mut again = project.attach_spectra("r1").unwrap();
    again.add(&ms3_spectrum(104, 2)).unwrap();
    again.add(&ms3_spectrum(105, 1)).unwrap();
    let report = again.finish().unwrap();
    assert_eq!((report.ms3_scans, report.ms3_attached), (2, 1));
    drop(project);

    let stored: Vec<(String, bool)> = rusqlite::Connection::open(&path)
        .unwrap()
        .prepare("SELECT spectrum_key, ms3_peaks IS NOT NULL FROM spectra ORDER BY spectrum_key")
        .unwrap()
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    // MS3 scans are never spectra of their own.
    assert_eq!(
        stored,
        [
            ("scan=1".to_string(), true),
            ("scan=2".to_string(), true),
            ("scan=3".to_string(), false)
        ]
    );
}

#[test]
fn attaching_to_an_unknown_run_lists_the_known_runs() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    let error = project.attach_spectra("HeLa_02").err().unwrap();

    assert!(matches!(error, StoreError::UnknownRun { .. }));
    assert!(error.to_string().contains("known runs: r1"), "{error}");
}

#[test]
fn an_unfinished_attach_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();

    {
        let mut attach = project.attach_spectra("r1").unwrap();
        attach.add(&peak_spectrum(1, None)).unwrap();
        // Dropped without finish(): rolled back.
    }

    assert_eq!(project.peaks("r1", &SpectrumKey::Scan(1)).unwrap(), None);
}

fn single_psm_batch(accession: &str, is_decoy: bool, decoy_explicit: bool) -> ImportBatch {
    ImportBatch {
        proteins: vec![Protein {
            accession: accession.into(),
            description: None,
            sequence: None,
            is_decoy,
            decoy_explicit,
        }],
        psms: vec![psm(
            "r1",
            1,
            PeptideForm::unmodified("PEPTIDEK"),
            &[accession],
            ("score", 1.0),
        )],
        ..ImportBatch::default()
    }
}

#[test]
fn explicit_decoy_flags_beat_the_accession_convention() {
    // "rev_P9" looks like a decoy, but mzIdentML says isDecoy="false".
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    project
        .import(
            &source("a.mzid", SourceFormat::MzIdentMl),
            &single_psm_batch("rev_P9", false, true),
        )
        .unwrap();
    // A later pepXML infers decoy from the same accession: it must not win.
    project
        .import(
            &source("b.pep.xml", SourceFormat::PepXml),
            &single_psm_batch("rev_P9", true, false),
        )
        .unwrap();

    let summary = project.summary().unwrap();
    assert_eq!((summary.decoy_proteins, summary.decoy_psms), (0, 0));
}

#[test]
fn explicit_decoy_flags_apply_to_earlier_psms() {
    // pepXML first: "SHD1" carries no decoy label, so it is a target...
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    project
        .import(
            &source("b.pep.xml", SourceFormat::PepXml),
            &single_psm_batch("SHD1", false, false),
        )
        .unwrap();
    assert_eq!(project.summary().unwrap().decoy_psms, 0);
    // ...until an mzIdentML flags it isDecoy="true": both PSMs become decoys.
    project
        .import(
            &source("a.mzid", SourceFormat::MzIdentMl),
            &single_psm_batch("SHD1", true, true),
        )
        .unwrap();

    let summary = project.summary().unwrap();
    assert_eq!((summary.decoy_proteins, summary.decoy_psms), (1, 2));
}

#[test]
fn repeated_loci_are_kept_and_repeated_scores_are_counted() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    let mut project = Project::create(&path).unwrap();
    let mut batch = single_psm_batch("P1", false, false);
    let locus = |start| Evidence {
        accession: "P1".into(),
        start: Some(start),
        end: Some(start + 7),
        pre: None,
        post: None,
    };
    // PEPTIDEK twice in P1 (two loci), the second locus repeated once, and
    // the score named twice.
    batch.psms[0].evidence = vec![locus(10), locus(50), locus(50)];
    batch.psms[0].scores.push(Score {
        name: "score".into(),
        accession: None,
        value: 99.0,
    });

    let report = project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &batch)
        .unwrap();

    assert_eq!((report.duplicate_evidence, report.duplicate_scores), (1, 1));
    drop(project);
    let loci: i64 = query(&path, "SELECT COUNT(*) FROM psm_proteins");
    let kept: f64 = query(&path, "SELECT value FROM psm_scores");
    assert_eq!((loci, kept), (2, 1.0));
}

#[test]
fn attached_peaks_are_never_replaced() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();
    let mut first = project.attach_spectra("r1").unwrap();
    first.add(&peak_spectrum(1, None)).unwrap();
    // The same file repeats scan 1 with different peaks.
    let mut repeat = peak_spectrum(1, None);
    repeat.mz = vec![999.0];
    repeat.intensity = vec![1.0];
    first.add(&repeat).unwrap();
    let report = first.finish().unwrap();
    assert_eq!((report.attached, report.already_attached), (1, 1));

    // A second file attaching scan 1 again changes nothing either.
    let mut second = project.attach_spectra("r1").unwrap();
    second.add(&repeat).unwrap();
    let report = second.finish().unwrap();
    assert_eq!((report.attached, report.already_attached), (0, 1));

    let stored = project.peaks("r1", &SpectrumKey::Scan(1)).unwrap().unwrap();
    assert_eq!(stored.mz, [100.5, 200.25]);
}

#[test]
fn a_precursors_ion_mobility_is_recorded_once_and_filled_in_later() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();
    let with_mobility = |scan: u32, mobility: f64| tuffold_core::PeakSpectrum {
        precursor_mobility: Some(mobility),
        ..peak_spectrum(scan, None)
    };
    let mut first = project.attach_spectra("r1").unwrap();
    first.add(&with_mobility(1, 1.0)).unwrap();
    first.add(&peak_spectrum(2, None)).unwrap();
    first.finish().unwrap();
    // Scan 2's peaks were attached without a mobility; attaching it again
    // fills it in, and scan 1's recorded mobility stays.
    let mut second = project.attach_spectra("r1").unwrap();
    second.add(&with_mobility(1, 1.2)).unwrap();
    second.add(&with_mobility(2, 1.1)).unwrap();
    assert_eq!(second.finish().unwrap().already_attached, 2);
    drop(project);
    let mobility = |key: &str| -> Option<f64> {
        query(
            &path,
            &format!("SELECT ion_mobility FROM spectra WHERE spectrum_key = '{key}'"),
        )
    };
    assert_eq!(mobility("scan=1"), Some(1.0));
    assert_eq!(mobility("scan=2"), Some(1.1));
    assert_eq!(mobility("scan=3"), None);
}

#[test]
fn partial_peaks_fill_gaps_and_give_way_to_a_spectra_file() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();
    let mut truncated = peak_spectrum(1, None);
    truncated.mz = vec![999.0];
    truncated.intensity = vec![1.0];
    let mut from_project = project.attach_partial_peaks("r1").unwrap();
    from_project.add(&truncated).unwrap();
    assert_eq!(from_project.finish().unwrap().attached, 1);
    // More partial peaks don't replace partial ones...
    let mut again = project.attach_partial_peaks("r1").unwrap();
    again.add(&peak_spectrum(1, None)).unwrap();
    assert_eq!(again.finish().unwrap().already_attached, 1);
    let stored = project.peaks("r1", &SpectrumKey::Scan(1)).unwrap().unwrap();
    assert_eq!(stored.mz, [999.0]);

    // ...a spectra file's do, and are then kept like any complete peaks.
    let mut file = project.attach_spectra("r1").unwrap();
    file.add(&peak_spectrum(1, None)).unwrap();
    assert_eq!(file.finish().unwrap().attached, 1);
    let stored = project.peaks("r1", &SpectrumKey::Scan(1)).unwrap().unwrap();
    assert_eq!(stored.mz, [100.5, 200.25]);
    let mut late = project.attach_partial_peaks("r1").unwrap();
    late.add(&truncated).unwrap();
    assert_eq!(late.finish().unwrap().already_attached, 1);
    let stored = project.peaks("r1", &SpectrumKey::Scan(1)).unwrap().unwrap();
    assert_eq!(stored.mz, [100.5, 200.25]);
}

#[test]
fn a_saved_copy_is_a_whole_project_with_an_undo_note() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("p.tuff");
    let project = Project::create(&path).unwrap();
    // SQLite reports the resolved path (macOS temp dirs are symlinked).
    let canonical = |p: &std::path::Path| p.canonicalize().unwrap();
    assert_eq!(project.path().map(canonical), Some(canonical(&path)));
    assert_eq!(project.undo_note().unwrap(), None);

    let copy_path = dir.path().join("p.tuff.undo");
    project.save_copy(&copy_path).unwrap();
    // Never over an existing file.
    assert!(project.save_copy(&copy_path).is_err());

    let mut copy = Project::open(&copy_path).unwrap();
    let note = tuffold_store::UndoNote {
        label: "import 2 files".into(),
        redo: false,
        saved_at: 1_700_000_000,
    };
    copy.set_undo_note(Some(&note)).unwrap();
    assert_eq!(copy.undo_note().unwrap(), Some(note.clone()));
    let redo = tuffold_store::UndoNote { redo: true, ..note };
    copy.set_undo_note(Some(&redo)).unwrap();
    assert_eq!(copy.undo_note().unwrap(), Some(redo));
    copy.set_undo_note(None).unwrap();
    assert_eq!(copy.undo_note().unwrap(), None);
    // The original never carried a note.
    assert_eq!(project.undo_note().unwrap(), None);
}

#[test]
fn copies_share_their_projects_identity() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("p.tuff")).unwrap();
    let id = project.identity().unwrap();
    assert_eq!(project.identity().unwrap(), id, "stable once made");
    project.save_copy(&dir.path().join("copy.tuff")).unwrap();
    let mut copy = Project::open(&dir.path().join("copy.tuff")).unwrap();
    assert_eq!(copy.identity().unwrap(), id);
    let mut other = Project::create(&dir.path().join("other.tuff")).unwrap();
    assert_ne!(other.identity().unwrap(), id);
    // The exclusive lock keeps a second connection from writing meanwhile.
    let second = Project::open(&dir.path().join("p.tuff")).unwrap();
    let blocked = project
        .exclusive(|| Project::open(&dir.path().join("p.tuff")).err())
        .unwrap();
    assert!(
        matches!(blocked, Some(tuffold_store::StoreError::Busy(_))),
        "{blocked:?}"
    );
    drop(second);
}

#[test]
fn removing_a_plex_frees_its_runs_and_keeps_channel_biosamples() {
    use tuffold_core::reporter::Plex;
    use tuffold_store::PlexDefinition;

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("study.tuff");
    let mut project = Project::create(&path).unwrap();
    project
        .import(&source("a.mzid", SourceFormat::MzIdentMl), &engine_a())
        .unwrap();
    let definition = PlexDefinition {
        name: "pool-1".into(),
        plex: Plex::Tmt6,
        runs: vec!["r1".into()],
        channels: vec![("126".into(), "wt-1".into()), ("127".into(), "ko-1".into())],
        impurities: None,
        references: vec!["126".into()],
    };
    project.set_plex(&definition).unwrap();
    assert_eq!(project.plexes().unwrap(), [definition]);

    let err = project.remove_plex("missing").unwrap_err();
    assert!(
        matches!(
            &err,
            StoreError::UnknownPlex { name, known } if name == "missing" && known == &["pool-1"]
        ),
        "{err}"
    );

    project.remove_plex("pool-1").unwrap();
    assert!(project.plexes().unwrap().is_empty());
    drop(project);
    assert_eq!(query::<i64>(&path, "SELECT COUNT(*) FROM plex_runs"), 0);
    assert_eq!(query::<i64>(&path, "SELECT COUNT(*) FROM plex_channels"), 0);
    // Channel biosamples stay as ordinary biosamples.
    let project = Project::open(&path).unwrap();
    let samples: Vec<String> = project
        .protein_list()
        .unwrap()
        .biosamples
        .into_iter()
        .map(|b| b.name)
        .collect();
    assert_eq!(samples, ["ko-1", "liver", "wt-1"]);
}
