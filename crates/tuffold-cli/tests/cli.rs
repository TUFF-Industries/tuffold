//! End-to-end CLI checks against the committed fixtures.
//!
//! Expected counts are hand-derived (see the importer oracle tests):
//! - FragPipe fixture: 2 biosamples, 2 runs, 5 spectra, 5 PSMs, 4 peptides,
//!   3 proteins (all with FASTA sequences), no decoys.
//! - pyteomics test.pep.xml: 1 run, 5 spectra, 6 PSMs (2 decoy), 8 proteins
//!   (3 decoy).
//! - pyteomics test.mzid: 1 run, 14 spectra, 22 PSMs (5 decoy), 19 peptides,
//!   7 proteins (2 decoy, all with sequences).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures")
}

fn tuffold(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_tuffold"))
        .args(args)
        .output()
        .expect("tuffold binary runs")
}

fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).unwrap()
}

fn stderr(output: &Output) -> String {
    String::from_utf8(output.stderr.clone()).unwrap()
}

fn path(relative: &str) -> String {
    fixtures().join(relative).to_str().unwrap().to_string()
}

#[test]
fn import_fragpipe_bundle_into_a_new_project() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("study.tuff");

    let output = tuffold(&[
        "import",
        &path("fragpipe_synthetic"),
        "--into",
        project.to_str().unwrap(),
    ]);

    assert!(output.status.success(), "{}", stderr(&output));
    let text = stdout(&output);
    assert!(
        text.contains("5 PSMs; new: 2 runs, 5 spectra, 4 peptides, 3 proteins"),
        "{text}"
    );
    assert!(text.contains(
        "project: 1 imports, 2 biosamples, 2 runs, 5 spectra, 5 PSMs (0 decoy), 4 peptides, \
         3 proteins (0 decoy, 3 with sequence)"
    ));
    assert!(text.contains("run control-1 (biosample control-1): 3 spectra, 3 PSMs"));
}

#[test]
fn import_all_three_formats_into_one_project() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("study.tuff");
    let project = project.to_str().unwrap();

    let output = tuffold(&[
        "import",
        &path("fragpipe_synthetic"),
        &path("third_party/pyteomics/test.pep.xml"),
        &path("third_party/pyteomics/test.mzid"),
        "-o",
        project,
    ]);

    assert!(output.status.success(), "{}", stderr(&output));
    // No shared runs or proteins between the three sources: totals add up.
    // PSMs 5 + 6 + 22 = 33; decoys 0 + 2 + 5 = 7; proteins 3 + 8 + 7 = 18.
    let text = stdout(&output);
    assert!(text.contains("project: 3 imports"), "{text}");
    assert!(text.contains("33 PSMs (7 decoy)"), "{text}");
    assert!(text.contains("18 proteins (5 decoy"), "{text}");

    let info = tuffold(&["info", project]);
    assert!(info.status.success());
    let info = stdout(&info);
    assert!(info.contains("[pepxml] SEQUEST"), "{info}");
    assert!(info.contains("[mzidentml] Sequest + Mascot"), "{info}");
}

#[test]
fn reimporting_the_same_input_fails_cleanly() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("study.tuff");
    let project = project.to_str().unwrap();
    let mzid = path("third_party/pyteomics/test.mzid");
    assert!(tuffold(&["import", &mzid, "-o", project]).status.success());

    let again = tuffold(&["import", &mzid, "-o", project]);

    assert!(!again.status.success());
    assert!(stderr(&again).contains("was already imported"));
}

#[test]
fn unknown_inputs_are_rejected_before_a_project_is_created() {
    let dir = tempfile::tempdir().unwrap();
    let notes = dir.path().join("notes.txt");
    std::fs::write(&notes, "hi").unwrap();
    let project = dir.path().join("study.tuff");

    let output = tuffold(&[
        "import",
        notes.to_str().unwrap(),
        "-o",
        project.to_str().unwrap(),
    ]);

    assert!(!output.status.success());
    assert!(stderr(&output).contains("unrecognized input"));
    assert!(
        !project.exists(),
        "no project file is created for a rejected input"
    );
}

#[test]
fn inspect_reports_ok_for_the_fixture() {
    let output = tuffold(&["inspect", &path("fragpipe_synthetic")]);
    assert!(output.status.success());
    assert!(stdout(&output).trim_end().ends_with("OK"));
}

#[test]
fn inspect_json_is_machine_readable() {
    let output = tuffold(&["inspect", &path("fragpipe_synthetic"), "--json"]);
    let json: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(
        json["experiment_folders"],
        serde_json::json!(["control-1", "treatment-1"])
    );
}

#[test]
fn invalid_bundles_fail_with_a_nonzero_exit() {
    let dir = tempfile::tempdir().unwrap();
    let output = tuffold(&["inspect", dir.path().to_str().unwrap()]);
    assert!(!output.status.success());
    assert!(stdout(&output).contains("error [missing-required-file]"));
}

#[test]
fn info_rejects_files_that_are_not_projects() {
    let output = tuffold(&["info", &path("third_party/pyteomics/test.mzid")]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("is not a Tuffold project"));
}

#[test]
fn spectra_attach_to_the_identified_scans_of_their_run() {
    // A pepXML whose run (base_name) matches tiny.pwiz.1.1.mzML and that
    // identifies scan 20 (its only MS/MS spectrum) and scan 99 (absent).
    let dir = tempfile::tempdir().unwrap();
    let pepxml = dir.path().join("search.pep.xml");
    std::fs::write(
        &pepxml,
        r#"<?xml version="1.0"?>
<msms_pipeline_analysis><msms_run_summary base_name="/data/tiny.pwiz.1.1" raw_data=".mzML">
<search_summary search_engine="Comet"/>
<spectrum_query spectrum="a" start_scan="20" assumed_charge="2"><search_result><search_hit hit_rank="1" peptide="PEPTIDEK" protein="P1"/></search_result></spectrum_query>
<spectrum_query spectrum="b" start_scan="99" assumed_charge="2"><search_result><search_hit hit_rank="1" peptide="PEPTIDER" protein="P1"/></search_result></spectrum_query>
</msms_run_summary></msms_pipeline_analysis>
"#,
    )
    .unwrap();
    let project = dir.path().join("study.tuff");

    let output = tuffold(&[
        "import",
        &path("third_party/pyteomics/tiny.pwiz.1.1.mzML"),
        pepxml.to_str().unwrap(),
        "-o",
        project.to_str().unwrap(),
    ]);

    assert!(output.status.success(), "{}", stderr(&output));
    let text = stdout(&output);
    assert!(
        text.contains(
            "to run tiny.pwiz.1.1: 1 of 1 MS/MS spectra matched identified spectra; \
             1 identified spectra still without peaks"
        ),
        "{text}"
    );
}

#[test]
fn spectra_for_an_unknown_run_name_the_known_runs() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("study.tuff");
    let output = tuffold(&[
        "import",
        &path("third_party/pyteomics/test.mzid"),
        &path("third_party/pyteomics/test.mgf"),
        "-o",
        project.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(
        stderr(&output).contains("no MS run named 'test'")
            && stderr(&output).contains("Fraction_X"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn validate_explains_engines_it_cannot_rank_and_info_shows_it() {
    // The pyteomics files report no score on every rank-1 PSM of an engine
    // (see tuffold-engine's validation tests), so nothing can be ranked.
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("study.tuff");
    let project = project.to_str().unwrap();
    let pep_xml = path("third_party/pyteomics/test.pep.xml");
    let mzid = path("third_party/pyteomics/test.mzid");
    assert!(
        tuffold(&["import", &pep_xml, &mzid, "--into", project])
            .status
            .success()
    );

    let info = stdout(&tuffold(&["info", project]));
    assert!(info.contains("validation: not yet run"), "{info}");

    let output = tuffold(&["validate", project, "--peptide-fdr", "5%"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let text = stdout(&output);
    assert!(
        text.contains("engine SEQUEST: 5 rank-1 PSMs (1 decoy)"),
        "{text}"
    );
    assert!(text.contains("thresholds: PSM FDR 1%, peptide FDR 5%, protein FDR 1%"));
    assert!(text.contains("note [engine-unscored] Sequest + Mascot PSMs were not validated"));

    let info = stdout(&tuffold(&["info", project]));
    assert!(info.contains("validation (rescoring on):"), "{info}");
}

#[test]
fn validating_without_installed_models_says_so() {
    // TUFFOLD_MODELS names an empty folder: no models there. A real (not
    // synthetic) project whose identified spectrum has peaks (scan 20 of
    // tiny.pwiz.1.1.mzML) would have used them, so the note says they
    // weren't; `--no-deep` (asked for), rescoring off, or a project without
    // peaks (nothing to predict) leave it out.
    let dir = tempfile::tempdir().unwrap();
    let models = dir.path().join("models");
    std::fs::create_dir(&models).unwrap();
    let pepxml = dir.path().join("search.pep.xml");
    std::fs::write(
        &pepxml,
        r#"<?xml version="1.0"?>
<msms_pipeline_analysis><msms_run_summary base_name="/data/tiny.pwiz.1.1" raw_data=".mzML">
<search_summary search_engine="Comet"/>
<spectrum_query spectrum="a" start_scan="20" assumed_charge="2"><search_result><search_hit hit_rank="1" peptide="PEPTIDEK" protein="P1"/></search_result></spectrum_query>
</msms_run_summary></msms_pipeline_analysis>
"#,
    )
    .unwrap();
    let with_peaks = dir.path().join("peaks.tuff");
    let without_peaks = dir.path().join("plain.tuff");
    let mzml = path("third_party/pyteomics/tiny.pwiz.1.1.mzML");
    for (project, inputs) in [
        (&with_peaks, vec![mzml.as_str(), pepxml.to_str().unwrap()]),
        (&without_peaks, vec![pepxml.to_str().unwrap()]),
    ] {
        let mut args = vec!["import", "--into", project.to_str().unwrap()];
        args.extend(inputs);
        let output = tuffold(&args);
        assert!(output.status.success(), "{}", stderr(&output));
    }
    let validate = |project: &Path, extra: &[&str]| {
        let output = Command::new(env!("CARGO_BIN_EXE_tuffold"))
            .arg("validate")
            .arg(project)
            .args(extra)
            .env("TUFFOLD_MODELS", &models)
            .output()
            .expect("tuffold binary runs");
        assert!(output.status.success(), "{}", stderr(&output));
        stdout(&output)
    };

    let text = validate(&with_peaks, &[]);
    let note = format!(
        "note [deep-learning] No deep-learning models are installed (looked in {})",
        models.display()
    );
    assert!(text.contains(&note), "{text}");
    for extra in [["--no-deep"], ["--no-rescore"]] {
        let text = validate(&with_peaks, &extra);
        assert!(!text.contains("[deep-learning]"), "{extra:?}: {text}");
    }
    let text = validate(&without_peaks, &[]);
    assert!(!text.contains("[deep-learning]"), "{text}");
}

#[test]
fn thresholds_need_a_validated_project_and_proper_rates() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("study.tuff");
    let project = project.to_str().unwrap();
    assert!(
        tuffold(&["import", &path("fragpipe_synthetic"), "--into", project])
            .status
            .success()
    );

    let output = tuffold(&["thresholds", project, "--protein-fdr", "5%"]);
    assert!(!output.status.success());
    assert!(
        stderr(&output).contains("has not been validated"),
        "{}",
        stderr(&output)
    );

    let output = tuffold(&["validate", project, "--psm-fdr", "1"]);
    assert!(!output.status.success());
    assert!(
        stderr(&output).contains("write 1% or 0.01"),
        "{}",
        stderr(&output)
    );
}

#[test]
fn demo_creates_a_validated_synthetic_project_once() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("demo.tuff");
    let project = project.to_str().unwrap();
    let output = tuffold(&["demo", "--into", project]);
    assert!(output.status.success(), "{}", stderr(&output));
    let text = stdout(&output);
    assert!(text.contains("from synthetic data"), "{text}");
    assert!(
        text.contains("engine Tuffold demo (synthetic): 4800 rank-1 PSMs"),
        "{text}"
    );
    assert!(
        text.contains("PTM localization:") && text.contains("localized at >= 95% site probability"),
        "{text}"
    );

    let loc_tsv = dir.path().join("localized.tsv");
    let loc_out = tuffold(&[
        "localize",
        project,
        "--site-prob",
        "90%",
        "--export",
        loc_tsv.to_str().unwrap(),
    ]);
    assert!(loc_out.status.success(), "{}", stderr(&loc_out));
    assert!(
        stdout(&loc_out).contains("localized at >= 90% site probability"),
        "{}",
        stdout(&loc_out)
    );
    // 20 ppm, from demo.rs: a correct spectrum keeps each ion with
    // probability `quality` (0.15..0.85, mean 0.5) at 3 ppm sigma, so 20 ppm
    // finds ~0.5 of its ions; a wrong one (peaks of an unseen peptide) ~0.
    // Correct are 0.55 of spectra and the 0.45 wrong are split between targets
    // and decoys, so ~0.55 / 0.775 = 0.71 of the target spectra sampled are
    // correct: 20 ppm finds ~0.35 of the sampled ions. 0.5 Da adds what an
    // empty 1 Th window catches by chance: 40 noise peaks plus at most ~60
    // other fragments over at least 400 Th (100 to max(top, 500)), at most
    // 0.25 of the ions. 2 x 0.35 = 0.7 > 0.35 + 0.25: 20 ppm.
    assert!(
        stdout(&loc_out).contains("PTM sites were localized with fragments matched within 20 ppm"),
        "{}",
        stdout(&loc_out)
    );
    assert!(text.contains("note [localization] PTM sites"), "{text}");
    let tsv = std::fs::read_to_string(&loc_tsv).unwrap();
    assert!(
        tsv.starts_with(
            "run\tbiosample\tspectrum\tcharge\tprecursor_mz\trt_seconds\tpeptide\tlocalized_peptide\tsite_probabilities\tbest_ascore\tq_value\tprobability\n"
        ),
        "{tsv}"
    );

    let again = tuffold(&["demo", "--into", project]);
    assert!(!again.status.success());
    assert!(
        stderr(&again).contains("already exists"),
        "{}",
        stderr(&again)
    );
}

#[test]
fn quant_compares_categories_and_exports_a_table() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("demo.tuff");
    let project = project.to_str().unwrap();
    assert!(tuffold(&["demo", "--into", project]).status.success());

    // The demo's biosamples come categorized; recategorize one and back.
    let output = tuffold(&["samples", project, "--set", "treated-B=knockout"]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(stdout(&output).contains("treated-B\tknockout"));
    assert!(
        tuffold(&["samples", project, "--set", "treated-B=treated"])
            .status
            .success()
    );

    let table = dir.path().join("quant.tsv");
    let output = tuffold(&[
        "quant",
        project,
        "--compare",
        "control",
        "treated",
        "--export",
        table.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{}", stderr(&output));
    let text = stdout(&output);
    assert!(
        text.contains("treated (treated-A, treated-B) vs control"),
        "{text}"
    );
    let tsv = std::fs::read_to_string(&table).unwrap();
    assert!(tsv.starts_with("accession\tcount:control-A"), "{tsv}");

    let output = tuffold(&["quant", project, "--compare", "control", "nope"]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("no biosample is in category 'nope'"));
    let output = tuffold(&["samples", project, "--set", "nonsense"]);
    assert!(stderr(&output).contains("write NAME=CATEGORY"));

    // A plex can be defined, listed, and removed with `--remove` (and undone).
    let defined = tuffold(&[
        "plex",
        project,
        "pool-1",
        "--label",
        "tmt6",
        "--run",
        "control-A",
        "126=c1",
        "127=t1",
    ]);
    assert!(defined.status.success(), "{}", stderr(&defined));
    assert!(stdout(&defined).contains("pool-1 (tmt6): runs control-A"));
    let missing = tuffold(&["plex", project, "nope", "--remove"]);
    assert!(!missing.status.success());
    assert!(stderr(&missing).contains("no plex named 'nope'"));
    let removed = tuffold(&["plex", project, "pool-1", "--remove"]);
    assert!(removed.status.success(), "{}", stderr(&removed));
    assert!(stdout(&removed).contains("no plexes; define one with --label"));
    assert!(tuffold(&["undo", project]).status.success());
    assert!(stdout(&tuffold(&["plex", project])).contains("pool-1 (tmt6)"));
}

#[test]
fn search_checks_inputs_then_searches_and_imports() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("searched.tuff");
    let project = project.to_str().unwrap();
    let fasta = path("fragpipe_synthetic/protein.fas");

    // A non-spectra input is refused before any project exists.
    let refused = tuffold(&[
        "search",
        &path("third_party/pyteomics/test.mzid"),
        "--fasta",
        &fasta,
        "-o",
        project,
    ]);
    assert!(!refused.status.success());
    assert!(
        stderr(&refused).contains("not a spectra file"),
        "{}",
        stderr(&refused)
    );
    assert!(!Path::new(project).exists());

    // Search options are checked before any project exists, too.
    let bad = tuffold(&[
        "search",
        &path("third_party/pyteomics/test.mgf"),
        "--fasta",
        &fasta,
        "-o",
        project,
        "--enzyme",
        "pepsin",
    ]);
    assert!(!bad.status.success());
    assert!(
        stderr(&bad).contains("unknown enzyme 'pepsin'"),
        "{}",
        stderr(&bad)
    );
    assert!(!Path::new(project).exists());
    // A fragment tolerance in ppm or in Da, never both.
    let both = tuffold(&[
        "search",
        &path("third_party/pyteomics/test.mgf"),
        "--fasta",
        &fasta,
        "-o",
        project,
        "--fragment-ppm",
        "20",
        "--fragment-da",
        "0.5",
    ]);
    assert!(!both.status.success());
    assert!(
        stderr(&both).contains("cannot be used with"),
        "{}",
        stderr(&both)
    );
    assert!(!Path::new(project).exists());

    // test.mgf holds 2 MS/MS spectra of 6 peaks each: too few to score
    // (15), so both are skipped and nothing is identified.
    let mgf = path("third_party/pyteomics/test.mgf");
    let output = tuffold(&["search", &mgf, "--fasta", &fasta, "-o", project]);
    assert!(output.status.success(), "{}", stderr(&output));
    let text = stdout(&output);
    // protein.fas: 5 entries, one of them a rev_ decoy (left out, regenerated).
    assert!(text.contains("digested 4 proteins"), "{text}");
    assert!(text.contains("search: trypsin, up to 2 missed"), "{text}");
    assert!(
        text.contains("1 decoy entries in the FASTA were left out"),
        "{text}"
    );
    assert!(
        text.contains("(run test): 0 spectra; no matches, nothing imported"),
        "{text}"
    );
    assert!(
        text.contains("not searched: 2 with fewer than 15 peaks"),
        "{text}"
    );
    assert!(text.contains("next: tuffold validate"), "{text}");

    // Nothing was imported, so the file can be searched again; a missing
    // file fails on its own, with a failure exit.
    let missing = dir.path().join("missing.mgf");
    let missing = missing.to_str().unwrap();
    let again = tuffold(&["search", &mgf, missing, "--fasta", &fasta, "-o", project]);
    assert!(!again.status.success());
    let text = stdout(&again);
    assert!(text.contains("(run test): 0 spectra"), "{text}");
    assert!(
        text.contains("failed") && text.contains("missing.mgf"),
        "{text}"
    );
}

#[test]
fn undo_reverses_the_last_import_then_redoes_it() {
    let dir = tempfile::tempdir().unwrap();
    let project = dir.path().join("p.tuff");
    let project = project.to_str().unwrap();
    let nothing = tuffold(&["undo", project]);
    assert!(!nothing.status.success());

    let pepxml = path("third_party/pyteomics/test.pep.xml");
    assert!(
        tuffold(&["import", &pepxml, "-o", project])
            .status
            .success()
    );
    let undone = tuffold(&["undo", project]);
    assert!(undone.status.success(), "{}", stderr(&undone));
    assert_eq!(
        stdout(&undone).trim(),
        "undid: import test.pep.xml (run tuffold undo again to redo it)"
    );
    // `info` lists one "import N: ..." line per import: none now.
    let info = stdout(&tuffold(&["info", project]));
    assert!(!info.contains("import 1:"), "{info}");
    let redone = tuffold(&["undo", project]);
    assert!(
        stdout(&redone).starts_with("redid: import test.pep.xml"),
        "{}",
        stdout(&redone)
    );
    let info = stdout(&tuffold(&["info", project]));
    assert!(
        info.contains("import 1:") && info.contains("test.pep.xml"),
        "{info}"
    );
}
