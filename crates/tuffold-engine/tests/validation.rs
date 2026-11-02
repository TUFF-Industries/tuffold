//! Validation through a real `.tuff` project: import, validate, change
//! thresholds, import again.

use std::path::Path;

use tuffold_core::{
    BioSample, Evidence, ImportBatch, MsRun, PeptideForm, Protein, Psm, Score, Software,
    SourceFormat, SpectrumKey,
};
use tuffold_engine::{ValidateError, apply_thresholds, import_inputs, validate_project};
use tuffold_store::{ImportSource, Project};
use tuffold_validate::{
    Level, NoteCode, SIBLINGS, SIBLINGS_IMPUTED, Settings, ThresholdMode, Thresholds,
};

const LETTERS: &[u8] = b"ACDEFGHILMNQSTVWY";
/// Proteins whose peptides the correct matches come from.
const PRESENT: u32 = 40;
const PROTEINS: u32 = 400;

fn noise(seed: u64) -> f64 {
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    ((z ^ (z >> 31)) >> 11) as f64 / (1_u64 << 53) as f64
}

fn sequence(protein: u32, k: u32, decoy: bool) -> String {
    let mut body: Vec<u8> = Vec::new();
    let mut rest = protein * 5 + k;
    for _ in 0..6 {
        body.push(LETTERS[(rest % 17) as usize]);
        rest /= 17;
    }
    if decoy {
        body.reverse();
    }
    format!("{}K", String::from_utf8(body).unwrap())
}

fn accession(protein: u32, decoy: bool) -> String {
    format!("{}sp|T{protein:03}|SYN", if decoy { "rev_" } else { "" })
}

/// A search of `run`: even scans are correct (scores 15..35), odd scans are
/// random targets or decoys (scores 5..25).
fn search(run: &str, scans: u32) -> ImportBatch {
    let mut psms = Vec::new();
    let mut proteins = Vec::new();
    for scan in 1..=scans {
        let seed = u64::from(scan) * 7 + u64::from(run.len() as u32);
        let correct = scan % 2 == 0;
        let decoy = !correct && noise(seed) < 0.5;
        let protein = if correct {
            (scan / 2) % PRESENT
        } else {
            (noise(seed + 1) * f64::from(PROTEINS)) as u32
        };
        let k = (scan / 2 / PRESENT) % 5;
        let hyperscore = if correct { 15.0 } else { 5.0 } + 20.0 * noise(seed + 2);
        let expect = 10_f64.powf(-hyperscore / 5.0 + noise(seed + 3));
        proteins.push(Protein {
            accession: accession(protein, decoy),
            description: None,
            sequence: None,
            is_decoy: decoy,
            decoy_explicit: false,
        });
        psms.push(Psm {
            ms_run: run.into(),
            spectrum: SpectrumKey::Scan(scan),
            rank: 1,
            peptide: PeptideForm::unmodified(&sequence(protein, k, decoy)),
            charge: 2,
            calculated_mz: None,
            experimental_mz: None,
            pass_threshold: None,
            scores: vec![
                Score {
                    name: "hyperscore".into(),
                    accession: None,
                    value: hyperscore,
                },
                Score {
                    name: "expect".into(),
                    accession: None,
                    value: expect,
                },
            ],
            evidence: vec![Evidence {
                accession: accession(protein, decoy),
                start: None,
                end: None,
                pre: None,
                post: None,
            }],
        });
    }
    ImportBatch {
        software: Software {
            name: "Synthetic".into(),
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
        spectra: Vec::new(),
        proteins,
        psms,
    }
}

fn project_with(dir: &Path, runs: &[&str]) -> Project {
    let mut project = Project::create(&dir.join("study.tuff")).unwrap();
    for run in runs {
        let source = ImportSource {
            path: format!("/searches/{run}.mzid"),
            format: SourceFormat::MzIdentMl,
        };
        project.import(&source, &search(run, 1200)).unwrap();
    }
    project
}

#[test]
fn a_validated_project_reports_counts_at_its_thresholds() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = project_with(dir.path(), &["r1"]);
    assert_eq!(project.validation_summary().unwrap(), None);

    let summary = validate_project(&mut project, &Settings::default()).unwrap();
    assert!(summary.validated_at.is_some());
    assert_eq!(summary.engines.len(), 1);
    let engine = &summary.engines[0];
    assert_eq!((engine.name.as_str(), engine.psms), ("Synthetic", 1200));
    // -log10(expect) is hyperscore / 5 minus noise, so hyperscore carries
    // all the signal and no trained model can beat it: rescoring must keep
    // hyperscore rather than cost identifications.
    assert_eq!(engine.method, "best_score", "{}", engine.detail);
    assert_eq!(engine.features, ["hyperscore"]);
    assert!(engine.detail.contains("fewer targets"), "{}", engine.detail);
    let codes: Vec<NoteCode> = summary.notes.iter().map(|note| note.code).collect();
    assert_eq!(codes, [NoteCode::RescoringFallback], "{:?}", summary.notes);
    // Correct matches come from 40 proteins with several peptides each, so
    // their peptides lift one another (decoys `rev_…` pair with their
    // targets), and the summary says what was pooled.
    let scoring = summary
        .peptide_scoring
        .as_ref()
        .expect("peptides were scored");
    assert_eq!(scoring.method, "pooled", "{}", scoring.detail);
    assert!(scoring.features.iter().any(|f| f == SIBLINGS));
    assert!(
        scoring.detail.contains("other identified peptides"),
        "{}",
        scoring.detail
    );
    // 600 correct spectra, whose scores overlap the incorrect ones by half
    // their range: a good share clears 1% FDR, with few decoys along.
    assert!(summary.spectra.targets > 250, "{:?}", summary.spectra);
    assert!(summary.spectra.decoys * 100 <= summary.spectra.targets * 2);
    // Correct matches come from 40 proteins (plus the odd random target a
    // 1% peptide FDR lets in; pooled, one decoy peptide passes among 201, as
    // 1% allows) and no decoy group outscores its target, so every group's
    // q-value is (0 + 1) / groups ≤ 1 / 40 = 2.5%, yet above 1% with fewer
    // than 100 groups: nothing can pass 1%, and the summary says how low the
    // estimate can go.
    assert!(summary.peptides.decoys * 100 <= summary.peptides.targets);
    assert_eq!(summary.protein_groups.targets, 0);
    let lowest = summary.protein_groups.lowest_q.unwrap();
    assert!(lowest > 0.01 && lowest <= 1.0 / 40.0, "{lowest}");

    // Results survive reopening.
    drop(project);
    let reopened = Project::open(&dir.path().join("study.tuff")).unwrap();
    assert_eq!(reopened.validation_summary().unwrap().unwrap(), summary);
}

/// Decoys named as TPP numbers them (`reverse-12`) name no target, so their
/// sibling counts are imputed from low-scoring targets'; the methods text
/// says so, and protein groups rebuilt for a new threshold rank as the
/// validation's did (by the stored group score). Wrong targets here land on
/// peptides the sample lacks (the body rotated by one: same length,
/// composition and missed cleavages, so no feature tells them from decoys),
/// so in a present protein they have siblings that decoys naming nothing
/// could not get.
#[test]
fn decoys_naming_no_target_get_imputed_siblings_and_groups_survive_rethresholding() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("study.tuff")).unwrap();
    let mut batch = search("r1", 1200);
    let rename = |accession: &mut String| {
        if let Some(rest) = accession.strip_prefix("rev_sp|T") {
            *accession = format!("reverse-{}", &rest[..3]);
        }
    };
    for protein in &mut batch.proteins {
        rename(&mut protein.accession);
    }
    for psm in &mut batch.psms {
        let SpectrumKey::Scan(scan) = psm.spectrum else {
            unreachable!("scans")
        };
        let decoy = psm.evidence[0].accession.starts_with("rev_");
        if scan % 2 == 1 && !decoy {
            let sequence = &psm.peptide.sequence;
            let (body, end) = sequence.split_at(sequence.len() - 1);
            psm.peptide = PeptideForm::unmodified(&format!("{}{}{end}", &body[1..], &body[..1]));
        }
        for evidence in &mut psm.evidence {
            rename(&mut evidence.accession);
        }
    }
    let source = ImportSource {
        path: "/searches/r1.mzid".into(),
        format: SourceFormat::MzIdentMl,
    };
    project.import(&source, &batch).unwrap();
    let settings = Settings {
        thresholds: Thresholds::default().with(Level::Protein, 0.05),
        ..Settings::default()
    };
    let validated = validate_project(&mut project, &settings).unwrap();
    assert!(
        validated
            .notes
            .iter()
            .any(|note| note.code == NoteCode::SiblingsImputed),
        "{:?}",
        validated.notes
    );
    let scoring = validated.peptide_scoring.as_ref().expect("scored");
    assert_eq!(scoring.method, "pooled", "{}", scoring.detail);
    assert!(scoring.features.iter().any(|f| f == SIBLINGS_IMPUTED));
    assert!(!scoring.features.iter().any(|f| f == SIBLINGS));
    assert!(scoring.detail.contains("imputed"), "{}", scoring.detail);
    let text = tuffold_engine::methods(&project).unwrap().text();
    assert!(text.contains("low-scoring targets"), "{text}");

    let loose = settings.thresholds.with(Level::Peptide, 0.2);
    apply_thresholds(&mut project, &loose).unwrap();
    let back = apply_thresholds(&mut project, &settings.thresholds).unwrap();
    assert_eq!(back.protein_groups, validated.protein_groups);
    assert!(
        validated.protein_groups.targets > 0,
        "{:?}",
        validated.protein_groups
    );
}

#[test]
fn thresholds_change_counts_without_revalidating() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = project_with(dir.path(), &["r1"]);
    let strict = validate_project(&mut project, &Settings::default()).unwrap();

    let loose = Thresholds::default()
        .with(Level::Psm, 0.05)
        .with(Level::Protein, 0.05);
    let summary = apply_thresholds(&mut project, &loose).unwrap();
    assert_eq!(summary.settings.thresholds, loose);
    assert_eq!(
        summary.validated_at, strict.validated_at,
        "no new validation run"
    );
    assert!(summary.spectra.targets > strict.spectra.targets);
    assert_eq!(
        summary.peptides, strict.peptides,
        "peptide threshold unchanged"
    );
    // The 40 proteins correct matches come from, and any random targets a
    // 5% PSM threshold lets in (no decoy group passes, so each is at
    // q = 1 / targets).
    assert!(
        summary.protein_groups.targets >= 40,
        "{:?}",
        summary.protein_groups
    );

    // Requiring more peptides than any group has passes nothing.
    let many = Thresholds {
        min_peptides: 6,
        ..loose
    };
    let summary = apply_thresholds(&mut project, &many).unwrap();
    assert_eq!(
        summary.protein_groups.targets, 0,
        "5 peptides per protein at most"
    );

    // A looser peptide threshold re-infers proteins from more peptides.
    let peptides = Thresholds::default().with(Level::Peptide, 0.2);
    let summary = apply_thresholds(&mut project, &peptides).unwrap();
    assert!(summary.peptides.targets > strict.peptides.targets);

    let bad = Thresholds::default().with(Level::Protein, 0.0);
    let error = apply_thresholds(&mut project, &bad).unwrap_err();
    assert_eq!(
        error.to_string(),
        "protein FDR must be above 0% and below 100%, not 0%"
    );

    // Probability thresholds read PEPs; switching mode re-infers proteins.
    let by_probability =
        apply_thresholds(&mut project, &Thresholds::probability_default()).unwrap();
    assert_eq!(
        by_probability.settings.thresholds.mode,
        ThresholdMode::Probability
    );
    assert_eq!(by_probability.validated_at, strict.validated_at);
    assert!(by_probability.spectra.targets > 0);
    // Every protein group here sits at q = 1 / groups (about 2.5%), and a PEP
    // is never below its q-value, so none can be 99% probable; at 97% the
    // 40 correct proteins' groups pass.
    assert_eq!(by_probability.protein_groups.targets, 0);
    let at_97 = Thresholds::probability_default().with(Level::Protein, 0.97);
    let summary = apply_thresholds(&mut project, &at_97).unwrap();
    assert!(
        summary.protein_groups.targets >= 40,
        "{:?}",
        summary.protein_groups
    );
}

#[test]
fn importing_more_data_makes_results_stale() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = project_with(dir.path(), &["r1"]);
    validate_project(&mut project, &Settings::default()).unwrap();

    let source = ImportSource {
        path: "/searches/r2.mzid".into(),
        format: SourceFormat::MzIdentMl,
    };
    project.import(&source, &search("r2", 1200)).unwrap();
    let stale = project.validation_summary().unwrap().unwrap();
    assert_eq!(stale.validated_at, None);
    assert_eq!(stale.spectra.targets, 0);
    assert!(stale.engines.is_empty());
    assert_eq!(stale.peptide_scoring, None);
    assert!(matches!(
        apply_thresholds(&mut project, &Thresholds::default()),
        Err(ValidateError::NotValidated)
    ));

    // Validating again covers both runs with the saved settings.
    let summary = validate_project(&mut project, &stale.settings).unwrap();
    assert_eq!(summary.engines[0].psms, 2400);
}

#[test]
fn the_pyteomics_fixtures_explain_why_they_cannot_be_validated() {
    let fixtures =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/third_party/pyteomics");
    let dir = tempfile::tempdir().unwrap();
    let mut project = Project::create(&dir.path().join("fixtures.tuff")).unwrap();
    let outcomes = import_inputs(
        &mut project,
        &[fixtures.join("test.pep.xml"), fixtures.join("test.mzid")],
        None,
    );
    assert!(outcomes.iter().all(|outcome| !outcome.is_failure()));

    // Neither file reports one score on all of an engine's rank-1 PSMs:
    // test.pep.xml has 5 rank-1 hits from mixed engines, and its most common
    // score, `expect`, is on 3 of them; test.mzid has 19, with
    // ProteinScape:IntensityCoverage on 11, SequestMetaScore on 9, and
    // MascotScore on 2. Nothing can be ranked, and each engine says why.
    let summary = validate_project(&mut project, &Settings::default()).unwrap();
    let methods: Vec<(&str, &str)> = summary
        .engines
        .iter()
        .map(|engine| (engine.name.as_str(), engine.method.as_str()))
        .collect();
    assert_eq!(
        methods,
        [("SEQUEST", "unscored"), ("Sequest + Mascot", "unscored")]
    );
    assert!(
        summary
            .notes
            .iter()
            .all(|note| note.code == NoteCode::EngineUnscored)
    );
    assert_eq!(summary.spectra.targets, 0);
}
