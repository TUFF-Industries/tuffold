//! Every tool against the demo project (and the committed fixtures).

use std::sync::OnceLock;

use super::*;

/// The demo project, built once; tests copy it before writing.
fn demo() -> &'static Path {
    static DEMO: OnceLock<(tempfile::TempDir, PathBuf)> = OnceLock::new();
    let (_, path) = DEMO.get_or_init(|| {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("demo.tuff");
        create_demo(&path).unwrap();
        (dir, path)
    });
    path
}

fn copy_of_demo(dir: &Path) -> PathBuf {
    let path = dir.join("copy.tuff");
    std::fs::copy(demo(), &path).unwrap();
    path
}

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures")
}

#[test]
fn the_summary_describes_contents_and_validation() {
    let summary = project_summary(demo()).unwrap();
    assert_eq!(summary["counts"]["biosamples"], 4);
    assert_eq!(summary["biosamples"][0]["name"], "control-A");
    assert_eq!(summary["biosamples"][0]["category"], "control");
    let validation = &summary["validation"];
    assert_eq!(validation["current"], true);
    assert_eq!(validation["engines"][0]["method"], "rescored");
    let groups = validation["passing"]["protein_groups"]["targets"]
        .as_u64()
        .unwrap();
    assert!(groups > 100);
    let missing = std::env::temp_dir().join("tuffold-missing.tuff");
    assert!(project_summary(&missing).is_err());
    // Relative paths would depend on the server's working directory.
    let relative = project_summary(Path::new("demo.tuff")).unwrap_err();
    assert!(relative.contains("absolute"), "{relative}");
}

#[test]
fn groups_page_filter_and_drill_down() {
    let all = list_protein_groups(demo(), None, true, 0, 5).unwrap();
    let total = all["total"].as_u64().unwrap();
    assert!(total > 100);
    assert_eq!(all["page"], json!({"offset": 0, "limit": 5, "returned": 5}));
    // Paging past the end returns nothing, and says where it stopped.
    let past = list_protein_groups(demo(), None, true, 10_000, 5).unwrap();
    assert_eq!(
        past["page"],
        json!({"offset": total, "limit": 5, "returned": 0})
    );
    // A limit over the maximum is capped, and the page says so.
    let capped = list_protein_groups(demo(), None, true, 0, 10 * MAX_ROWS).unwrap();
    assert_eq!(capped["page"]["limit"], MAX_ROWS);

    let first = &all["groups"][0];
    let accession = first["accession"].as_str().unwrap();
    let filtered =
        list_protein_groups(demo(), Some(&accession.to_lowercase()), true, 0, 50).unwrap();
    let groups = filtered["groups"].as_array().unwrap();
    assert!(groups.iter().any(|g| g["accession"] == accession));

    let group = protein_group(demo(), first["group_id"].as_i64().unwrap()).unwrap();
    assert_eq!(group["members"][0]["accession"], accession);
    let coverage = group["coverage"].as_f64().unwrap();
    assert!((coverage - first["coverage"].as_f64().unwrap()).abs() < 1e-9);
    let peptide = group["peptides"][0]["peptide_id"].as_i64().unwrap();
    let spectra = peptide_spectra(demo(), peptide, 2).unwrap();
    assert!(spectra["page"]["returned"].as_u64().unwrap() <= 2);
    assert!(spectra["total"].as_u64().unwrap() >= 1);
    assert!(protein_group(demo(), -1).unwrap_err().contains("-1"));
}

#[test]
fn categories_compare_and_change() {
    let compared = compare_categories(demo(), "control", "treated", 0, 10).unwrap();
    assert_eq!(compared["primary_test"], PrimaryTest::Welch.describe());
    assert_eq!(
        compared["versus_samples"],
        json!(["treated-A", "treated-B"])
    );
    assert_eq!(compared["rows"].as_array().unwrap().len(), 10);
    // The next page starts where this one ended.
    let next = compare_categories(demo(), "control", "treated", 10, 10).unwrap();
    assert_eq!(next["page"]["offset"], 10);
    assert_ne!(next["rows"][0], compared["rows"][0]);
    assert!(compare_categories(demo(), "control", "nope", 0, 10).is_err());
    // By MS1 intensity: the demo's runs carry MS1 scans.
    let by_intensity = compare_intensities(demo(), "control", "treated", 0, 5).unwrap();
    assert_eq!(by_intensity["method"], "ms1_intensity");
    assert_eq!(by_intensity["runs"][0]["ms1_scans"], 1200);
    assert!(by_intensity["significant_at_5pct"].as_u64().unwrap() >= 5);
    assert_eq!(by_intensity["rows"].as_array().unwrap().len(), 5);
    let all = compare_all_categories(demo(), "intensity", None, 0, 3).unwrap();
    assert_eq!(all["categories"][1]["category"], "treated");
    assert!(all["significant_at_5pct"].as_u64().unwrap() >= 5);
    assert!(compare_all_categories(demo(), "vibes", None, 0, 3).is_err());
    let qc = sample_qc(demo(), "counts", None).unwrap();
    let samples = qc["samples"].as_array().unwrap();
    assert_eq!(samples.len(), 4);
    assert_eq!(samples[0]["correlation"].as_array().unwrap().len(), 4);
    assert!(qc["pca"]["variance_explained"][0].as_f64().unwrap() > 0.0);
    assert!(samples[0]["pc1"].is_number());
    assert!(sample_qc(demo(), "vibes", None).is_err());
    let dir = tempfile::tempdir().unwrap();
    let workbook = dir.path().join("demo.xlsx");
    let written = export_workbook(demo(), &workbook).unwrap();
    assert_eq!(written["sheets"][0]["sheet"], "Proteins");
    assert!(
        export_workbook(demo(), &workbook).is_err(),
        "never overwrites"
    );
    assert!(export_workbook(demo(), &dir.path().join("x.csv")).is_err());
    let methods = methods_text(demo()).unwrap();
    assert!(methods["text"].as_str().unwrap().contains("References"));

    let dir = tempfile::tempdir().unwrap();
    let path = copy_of_demo(dir.path());
    let changed = set_categories(&path, &[("treated-B".into(), String::new())]).unwrap();
    assert_eq!(changed["biosamples"][3]["category"], Value::Null);
    // Undoable: back, then forward again.
    let undone = undo_last_change(&path).unwrap();
    assert_eq!(undone["done"], "undid");
    assert_eq!(undone["change"], "change sample groups");
    let summary = project_summary(&path).unwrap();
    assert_eq!(summary["biosamples"][3]["category"], "treated");
    assert_eq!(undo_last_change(&path).unwrap()["done"], "redid");
    // One treated sample left: Fisher on pooled counts decides.
    let compared = compare_categories(&path, "control", "treated", 0, 3).unwrap();
    assert_eq!(compared["primary_test"], PrimaryTest::Fisher.describe());
}

#[test]
fn plex_defines_and_quantifies_reporters() {
    use tuffold_core::reporter::Plex;
    use tuffold_core::{
        BioSample, Evidence, ImportBatch, MsRun, PeakSpectrum, PeptideForm, Protein, Psm, Score,
        Software, SourceFormat, Spectrum, SpectrumKey,
    };
    use tuffold_store::ImportSource;
    use tuffold_validate::ThresholdMode;

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("tmt.tuff");
    let mut project = Project::create(&path).unwrap();
    let run = "pool_1";
    let psm = |scan: u32, seq: &str, acc: &str, score: f64| Psm {
        ms_run: run.into(),
        spectrum: SpectrumKey::Scan(scan),
        rank: 1,
        peptide: PeptideForm::unmodified(seq),
        charge: 2,
        calculated_mz: Some(500.0),
        experimental_mz: Some(500.0),
        pass_threshold: None,
        scores: vec![Score {
            name: "hyperscore".into(),
            accession: None,
            value: score,
        }],
        evidence: vec![Evidence {
            accession: acc.into(),
            start: None,
            end: None,
            pre: None,
            post: None,
        }],
    };
    let batch = ImportBatch {
        software: Software {
            name: "Sage".into(),
            version: None,
        },
        biosamples: vec![BioSample {
            name: run.into(),
            category: None,
        }],
        ms_runs: vec![MsRun {
            name: run.into(),
            biosample: run.into(),
            spectra_file: None,
        }],
        spectra: (1..=3)
            .map(|s| Spectrum {
                ms_run: run.into(),
                key: SpectrumKey::Scan(s),
                source_id: None,
                precursor_mz: Some(500.0),
                precursor_charge: Some(2),
                rt_seconds: Some(f64::from(s) * 10.0),
            })
            .collect(),
        psms: vec![
            psm(1, "PEPTIDEK", "P1", 50.0),
            psm(2, "ANOTHERK", "P2", 40.0),
            psm(3, "DECOYPEPK", "rev_P1", 10.0),
        ],
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: None,
                sequence: None,
                is_decoy: false,
                decoy_explicit: false,
            },
            Protein {
                accession: "P2".into(),
                description: None,
                sequence: None,
                is_decoy: false,
                decoy_explicit: false,
            },
            Protein {
                accession: "rev_P1".into(),
                description: None,
                sequence: None,
                is_decoy: true,
                decoy_explicit: false,
            },
        ],
    };
    project
        .import(
            &ImportSource {
                path: "/data/tmt.pep.xml".into(),
                format: SourceFormat::PepXml,
            },
            &batch,
        )
        .unwrap();
    let mz: Vec<f64> = Plex::Tmt6Original
        .channels()
        .into_iter()
        .take(4)
        .map(|c| c.mz)
        .collect();
    let scan = |s: u32, intensities: [f32; 4]| PeakSpectrum {
        key: SpectrumKey::Scan(s),
        native_id: format!("scan={s}"),
        ms_level: 2,
        rt_seconds: None,
        precursor_mz: Some(500.0),
        precursor_charge: Some(2),
        mz: mz.clone(),
        intensity: intensities.to_vec(),
        profile: false,
        parent: None,
        isolation: None,
        precursor_mobility: None,
        mobility: Vec::new(),
    };
    // P1: [1000, 2000, 4000, 16000]; P2: [2000, 1000, 2000, 2000].
    // Channel 129 has 2x loading relative to 128 (shift_128 - shift_129 = +1.0 log2),
    // and after normalization P1 is +1.0 log2 and P2 is -1.0 log2 with nonzero within-group spread.
    let mut attacher = project.attach_spectra(run).unwrap();
    attacher
        .add(&scan(1, [1000.0, 2000.0, 4000.0, 16000.0]))
        .unwrap();
    attacher
        .add(&scan(2, [2000.0, 1000.0, 2000.0, 2000.0]))
        .unwrap();
    attacher.finish().unwrap();
    validate_project(
        &mut project,
        &Settings {
            rescore: false,
            thresholds: Thresholds {
                mode: ThresholdMode::Fdr,
                psm: 0.9,
                peptide: 0.9,
                protein: 0.9,
                min_peptides: 1,
            },
        },
    )
    .unwrap();
    drop(project);

    let channels = vec![
        ("126".into(), "wt-1".into()),
        ("127".into(), "wt-2".into()),
        ("128".into(), "ko-1".into()),
        ("129".into(), "ko-2".into()),
    ];
    // Defined with tmt6 first: 127 and 128 are at different m/z from tmt6-original,
    // so 0 PSMs have all 4 channels and better_reagents points to tmt6-original (2 PSMs).
    // Blank names and assignments are the engine's to trim and the store's to refuse.
    let blank = [(" 126 ".to_string(), "  ".to_string())];
    assert_eq!(
        set_plex(&path, "  ", "tmt6", &[], &channels, None, None).unwrap_err(),
        "name the plex"
    );
    assert_eq!(
        set_plex(&path, "pool", "tmt6", &[], &blank, None, None).unwrap_err(),
        "assign at least one channel to a biosample"
    );
    let defined = set_plex(&path, " pool ", "tmt6", &[], &channels, None, None).unwrap();
    assert_eq!(defined["plexes"][0]["name"], "pool");
    assert_eq!(defined["plexes"][0]["label"], "tmt6");
    assert_eq!(defined["plexes"][0]["runs"], json!(["pool_1"]));
    let wrong = quantify_reporters(&path, None, None, None, 0, 50).unwrap();
    assert_eq!(wrong["psms_used"], 0);
    assert_eq!(wrong["better_reagents"][0]["label"], "tmt6-original");
    assert_eq!(wrong["better_reagents"][0]["complete_psms"], 2);

    // Redefine with tmt6-original and a pure lot's certificate (all zeros,
    // so the numbers below are unchanged), then categorize the channel
    // biosamples.
    let certificate = path.with_extension("certificate.txt");
    std::fs::write(
        &certificate,
        "126 0 0 0 0\n127 0 0 0 0\n128 0 0 0 0\n129 0 0 0 0\n130 0 0 0 0\n131 0 0 0 0\n",
    )
    .unwrap();
    let file = certificate.to_str().unwrap();
    let defined = set_plex(
        &path,
        "pool",
        "tmt6-original",
        &[],
        &channels,
        Some(file),
        None,
    )
    .unwrap();
    assert_eq!(defined["plexes"][0]["impurity_corrected"], true);
    // Saved again without one, the plex keeps it.
    let kept = set_plex(&path, "pool", "tmt6-original", &[], &channels, None, None).unwrap();
    assert_eq!(kept["plexes"][0]["impurity_corrected"], true);
    // Reference channels are set, then kept when omitted (a lone plex still
    // quantifies its reference as a sample).
    let reference = ["126".to_string()];
    let referenced = set_plex(
        &path,
        "pool",
        "tmt6-original",
        &[],
        &channels,
        None,
        Some(&reference),
    )
    .unwrap();
    assert_eq!(referenced["plexes"][0]["references"], json!(["126"]));
    let kept = set_plex(&path, "pool", "tmt6-original", &[], &channels, None, None).unwrap();
    assert_eq!(kept["plexes"][0]["references"], json!(["126"]));
    set_categories(
        &path,
        &[
            ("wt-1".into(), "wt".into()),
            ("wt-2".into(), "wt".into()),
            ("ko-1".into(), "ko".into()),
            ("ko-2".into(), "ko".into()),
        ],
    )
    .unwrap();

    let listed = quantify_reporters(&path, Some("pool"), None, None, 0, 50).unwrap();
    assert_eq!(listed["psms_used"], 2);
    assert_eq!(listed["better_reagents"], json!([]));
    let shift_128 = listed["channels"][2]["log2_shift"].as_f64().unwrap();
    let shift_129 = listed["channels"][3]["log2_shift"].as_f64().unwrap();
    assert!((shift_128 - shift_129 - 1.0).abs() < 1e-9);

    let compared = quantify_reporters(&path, None, Some("wt"), Some("ko"), 0, 10).unwrap();
    assert_eq!(compared["method"], "reporters");
    let rows = compared["rows"].as_array().unwrap();
    let fold_of = |acc: &str| {
        rows.iter()
            .find(|r| r["accession"] == acc)
            .and_then(|r| r["log2_fold"].as_f64())
            .unwrap()
    };
    assert!((fold_of("P1") - 1.0).abs() < 1e-9);
    assert!((fold_of("P2") + 1.0).abs() < 1e-9);

    let anova = compare_all_categories(&path, "reporters", Some("pool"), 0, 10).unwrap();
    assert_eq!(anova["test"], "one-way ANOVA on log2 reporter intensities");
    assert_eq!(anova["tested"], 2);

    assert!(remove_plex(&path, "missing").is_err());
    let removed = remove_plex(&path, "pool").unwrap();
    assert_eq!(removed["plexes"], json!([]));
    assert!(
        removed["biosamples"]
            .as_array()
            .unwrap()
            .iter()
            .any(|b| b["name"] == "wt-1")
    );
}

#[test]
fn validation_reruns_and_thresholds_apply() {
    let dir = tempfile::tempdir().unwrap();
    let path = copy_of_demo(dir.path());
    let strict = ThresholdRequest {
        fdr: [None, None, Some(0.001)],
        ..ThresholdRequest::default()
    };
    let applied = set_thresholds(&path, &strict).unwrap();
    assert_eq!(applied["thresholds"]["by"], "fdr");
    assert_eq!(applied["thresholds"]["protein"], 0.001);
    // Probability thresholds: values are probabilities, and each level says
    // what FDR the decoys put on what passes.
    let by_probability = ThresholdRequest {
        by: Some(tuffold_validate::ThresholdMode::Probability),
        ..ThresholdRequest::default()
    };
    let applied = set_thresholds(&path, &by_probability).unwrap();
    assert_eq!(applied["thresholds"]["by"], "probability");
    assert_eq!(applied["thresholds"]["protein"], 0.99);
    assert!(applied["passing"]["protein_groups"]["estimated_fdr"].is_number());
    let groups = list_protein_groups(&path, None, true, 0, 1).unwrap();
    assert!(groups["groups"][0]["probability"].as_f64().unwrap() >= 0.99);
    let ranked = validate(&path, &ThresholdRequest::default(), false).unwrap();
    assert_eq!(ranked["rescore"], false);
    assert_eq!(ranked["engines"][0]["method"], "best_score");
}

#[test]
fn imports_report_each_input_and_leave_nothing_when_all_fail() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("new.tuff");
    let imported = import_data(
        &path,
        &[
            fixtures().join("third_party/pyteomics/test.pep.xml"),
            dir.path().join("nothing.txt"),
        ],
    )
    .unwrap();
    let statuses: Vec<&str> = imported["outcomes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|o| o["status"].as_str().unwrap())
        .collect();
    assert_eq!(statuses, ["failed", "imported"]);
    // test.pep.xml: 6 PSMs (see the importer oracle tests).
    assert_eq!(imported["outcomes"][1]["psms"], 6);

    // Nothing importable: an error, and no empty project left behind.
    let stray = dir.path().join("stray.tuff");
    let failed = import_data(&stray, &[dir.path().join("nothing.txt")]).unwrap_err();
    assert!(failed.starts_with("nothing was imported"), "{failed}");
    assert!(!stray.exists());
}

#[test]
fn searches_that_match_nothing_are_errors_and_leave_nothing() {
    // test.mgf: 2 spectra of 6 peaks, too few to search; protein.fas has
    // 4 targets and 1 decoy.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("searched.tuff");
    let failed = search_spectra(
        &path,
        &fixtures().join("fragpipe_synthetic/protein.fas"),
        &[fixtures().join("third_party/pyteomics/test.mgf")],
        SearchSettings::default(),
    )
    .unwrap_err();
    assert!(failed.starts_with("nothing was imported"), "{failed}");
    assert!(failed.contains("\"target_proteins\":4"), "{failed}");
    assert!(failed.contains("\"not_searched\":2"), "{failed}");
    assert!(!path.exists());

    assert!(create_demo(demo()).unwrap_err().contains("already exists"));
}
