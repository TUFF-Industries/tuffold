//! Browse queries on a hand-built project with hand-built validation
//! results, so every count below can be read off the setup.
//!
//! Protein sp|P1|ONE = MPEPTIDEKAGSAMPLERK (19 residues); its decoy rev_P1.
//! - control/r1: scan 1 and scan 2 assigned PEPTIDEK (q 0.001 and 0.02);
//! - treated/r2: scan 1 assigned AGSAMPLER (q 0.005); scan 2 the decoy
//!   KEDITPEP (q 0.5).
//!
//! At the default 1% thresholds scan r1/2 fails, so P1 has 1 control and 1
//! treated spectrum. Coverage: PEPTIDEK at 2-9 and AGSAMPLER at 10-18 cover
//! 17 of 19 residues.

use std::path::Path;

use rusqlite::Connection;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, MsRun, PeakSpectrum, PeptideForm, Protein, Psm, Score,
    Software, SourceFormat, SpectrumKey,
};
use tuffold_store::{ImportSource, Project, StoredPeaks};
use tuffold_validate::{
    MemberRole, PeptideResult, PeptideRole, ProteinGroup, ProteinInference, Settings,
    SpectrumMatch, Validation,
};

const SEQUENCE: &str = "MPEPTIDEKAGSAMPLERK";

fn psm(run: &str, scan: u32, peptide: &str, accession: &str) -> Psm {
    Psm {
        ms_run: run.into(),
        spectrum: SpectrumKey::Scan(scan),
        rank: 1,
        peptide: PeptideForm::unmodified(peptide),
        charge: 2,
        calculated_mz: None,
        experimental_mz: None,
        pass_threshold: None,
        scores: vec![Score {
            name: "hyperscore".into(),
            accession: None,
            value: 30.0,
        }],
        evidence: vec![Evidence {
            accession: accession.into(),
            start: None,
            end: None,
            pre: None,
            post: None,
        }],
    }
}

fn batch() -> ImportBatch {
    let run = |name: &str, sample: &str| MsRun {
        name: name.into(),
        biosample: sample.into(),
        spectra_file: None,
    };
    ImportBatch {
        software: Software {
            name: "Engine".into(),
            version: None,
        },
        biosamples: ["control", "treated"]
            .map(|name| BioSample {
                name: name.into(),
                category: None,
            })
            .to_vec(),
        ms_runs: vec![run("r1", "control"), run("r2", "treated")],
        spectra: Vec::new(),
        proteins: vec![
            Protein {
                accession: "sp|P1|ONE".into(),
                description: Some("Protein one".into()),
                sequence: Some(SEQUENCE.into()),
                is_decoy: false,
                decoy_explicit: false,
            },
            Protein {
                accession: "rev_sp|P1|ONE".into(),
                description: None,
                sequence: None,
                is_decoy: true,
                decoy_explicit: false,
            },
        ],
        psms: vec![
            psm("r1", 1, "PEPTIDEK", "sp|P1|ONE"),
            psm("r1", 2, "PEPTIDEK", "sp|P1|ONE"),
            psm("r2", 1, "AGSAMPLER", "sp|P1|ONE"),
            psm("r2", 2, "KEDITPEP", "rev_sp|P1|ONE"),
        ],
    }
}

fn id(conn: &Connection, sql: &str) -> i64 {
    conn.query_row(sql, [], |row| row.get(0)).unwrap()
}

/// Import the batch, attach peaks to r1 scan 1, and save hand-made results.
fn project(dir: &Path) -> Project {
    let path = dir.join("browse.tuff");
    let mut project = Project::create(&path).unwrap();
    let source = ImportSource {
        path: "/s/a.mzid".into(),
        format: SourceFormat::MzIdentMl,
    };
    project.import(&source, &batch()).unwrap();
    let mut attacher = project.attach_spectra("r1").unwrap();
    attacher
        .add(&PeakSpectrum {
            key: SpectrumKey::Scan(1),
            native_id: "scan=1".into(),
            ms_level: 2,
            rt_seconds: Some(60.0),
            precursor_mz: Some(465.73),
            precursor_charge: Some(2),
            mz: vec![148.06, 227.10],
            intensity: vec![10.0, 20.0],
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        })
        .unwrap();
    attacher.finish().unwrap();

    let conn = Connection::open(&path).unwrap();
    let spectrum = |run: &str, scan: u32| {
        id(
            &conn,
            &format!(
                "SELECT s.id FROM spectra s JOIN ms_runs r ON r.id = s.ms_run_id
             WHERE r.name = '{run}' AND s.spectrum_key = 'scan={scan}'"
            ),
        )
    };
    let peptide = |sequence: &str| {
        id(
            &conn,
            &format!("SELECT id FROM peptides WHERE sequence = '{sequence}'"),
        )
    };
    let (pep, agg, dec) = (
        peptide("PEPTIDEK"),
        peptide("AGSAMPLER"),
        peptide("KEDITPEP"),
    );
    let p1 = id(
        &conn,
        "SELECT id FROM proteins WHERE accession = 'sp|P1|ONE'",
    );
    let found = |spectrum_id, peptide_id, is_decoy, score, q_value| SpectrumMatch {
        spectrum_id,
        peptide_id,
        is_decoy,
        score,
        q_value,
        // PEPs equal to the q-values keep the fixture easy to reason about.
        pep: q_value,
        engines: 1,
    };
    let validation = Validation {
        spectra: vec![
            found(spectrum("r1", 1), pep, false, 5.0, 0.001),
            found(spectrum("r1", 2), pep, false, 2.0, 0.02),
            found(spectrum("r2", 1), agg, false, 3.0, 0.005),
            found(spectrum("r2", 2), dec, true, 0.5, 0.5),
        ],
        peptides: vec![
            PeptideResult {
                peptide_id: pep,
                is_decoy: false,
                score: 5.0,
                group_score: 5.0,
                q_value: 0.001,
                pep: 0.0,
                best_spectrum_id: spectrum("r1", 1),
            },
            PeptideResult {
                peptide_id: agg,
                is_decoy: false,
                score: 3.0,
                group_score: 3.0,
                q_value: 0.005,
                pep: 0.0,
                best_spectrum_id: spectrum("r2", 1),
            },
        ],
        proteins: ProteinInference {
            groups: vec![ProteinGroup {
                leader: p1,
                members: vec![(p1, MemberRole::Leader)],
                peptides: vec![(pep, PeptideRole::Unique), (agg, PeptideRole::Unique)],
                is_decoy: false,
                score: 5.0,
                q_value: Some(0.001),
                pep: Some(0.002),
                distinct_sequences: 2,
            }],
            notes: Vec::new(),
        },
        ..Validation::default()
    };
    project
        .save_validation(&Settings::default(), &validation)
        .unwrap();
    project
}

#[test]
fn protein_list_counts_passing_spectra_per_biosample_and_coverage() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let list = project.protein_list().unwrap();
    let samples: Vec<&str> = list.biosamples.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(samples, ["control", "treated"]);
    assert_eq!(list.rows.len(), 1);
    let row = &list.rows[0];
    assert_eq!(
        (row.accession.as_str(), row.description.as_deref()),
        ("sp|P1|ONE", Some("Protein one"))
    );
    assert_eq!((row.peptides, row.spectra, row.equivalents), (2, 2, 0));
    assert_eq!(row.spectra_by_sample, [1, 1]);
    assert!((row.coverage.unwrap() - 17.0 / 19.0).abs() < 1e-12);
    assert_eq!(row.length, Some(19));
    assert_eq!(row.q_value, Some(0.001));
    // Both peptides are unique: every count is exclusive. One passing
    // spectrum each (r1/1 PEPTIDEK in control, r2/1 AGSAMPLER in treated).
    assert_eq!(row.exclusive_spectra_by_sample, [1, 1]);
    assert_eq!(row.exclusive_peptides_by_sample, [1, 1]);
    assert_eq!(row.peptides_by_sample, [1, 1]);
    // MPEPTIDEK | AGSAMPLER | K: two tryptic peptides of 6-30 residues.
    assert_eq!(row.observable_peptides, Some(2));
}

#[test]
fn razor_peptides_count_but_are_not_exclusive() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let conn = Connection::open(dir.path().join("browse.tuff")).unwrap();
    conn.execute(
        "UPDATE protein_group_peptides SET role = 'razor'
         WHERE peptide_id = (SELECT id FROM peptides WHERE sequence = 'AGSAMPLER')",
        [],
    )
    .unwrap();
    let row = &project.protein_list().unwrap().rows[0];
    // AGSAMPLER (treated's only spectrum) is now razor: still counted,
    // no longer exclusive.
    assert_eq!(row.spectra_by_sample, [1, 1]);
    assert_eq!(row.exclusive_spectra_by_sample, [1, 0]);
    assert_eq!(row.exclusive_peptides_by_sample, [1, 0]);
    assert_eq!(row.peptides_by_sample, [1, 1]);
}

#[test]
fn group_detail_lists_members_peptides_and_where_they_fall() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let group = project.protein_list().unwrap().rows[0].group_id;
    let detail = project.group_detail(group).unwrap().unwrap();
    assert_eq!(detail.members.len(), 1);
    assert_eq!(detail.members[0].role, "leader");
    assert_eq!(detail.sequence.as_deref(), Some(SEQUENCE));
    let peptides: Vec<(&str, &str, usize)> = detail
        .peptides
        .iter()
        .map(|p| (p.sequence.as_str(), p.role.as_str(), p.spectra))
        .collect();
    assert_eq!(
        peptides,
        [("PEPTIDEK", "unique", 1), ("AGSAMPLER", "unique", 1)]
    );
    let spans: Vec<(u32, u32)> = detail.coverage.iter().map(|c| (c.start, c.end)).collect();
    assert_eq!(spans, [(2, 9), (10, 18)]);
    assert_eq!(project.group_detail(group + 100).unwrap(), None);
}

#[test]
fn peptide_spectra_flag_the_ones_passing_the_threshold() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let group = project.protein_list().unwrap().rows[0].group_id;
    let detail = project.group_detail(group).unwrap().unwrap();
    let page = project
        .peptide_spectra(detail.peptides[0].peptide_id, 10)
        .unwrap();
    assert_eq!(page.total, 2);
    assert!(!page.is_capped());
    let spectra = page.rows;
    let rows: Vec<(&str, &str, bool, bool)> = spectra
        .iter()
        .map(|s| (s.run.as_str(), s.key.as_str(), s.passes, s.has_peaks))
        .collect();
    assert_eq!(
        rows,
        [("r1", "scan=1", true, true), ("r1", "scan=2", false, false)]
    );
    assert_eq!(spectra[0].biosample, "control");

    // Capped: the best one only, and the total still says there are two.
    let capped = project
        .peptide_spectra(detail.peptides[0].peptide_id, 1)
        .unwrap();
    assert_eq!(capped.rows, spectra[..1]);
    assert_eq!(capped.total, 2);
    assert!(capped.is_capped());
}

#[test]
fn spectrum_detail_brings_peaks_the_assigned_peptide_and_engine_scores() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let group = project.protein_list().unwrap().rows[0].group_id;
    let peptide = project.group_detail(group).unwrap().unwrap().peptides[0].peptide_id;
    let spectra = project.peptide_spectra(peptide, 10).unwrap().rows;

    let detail = project
        .spectrum_detail(spectra[0].spectrum_id)
        .unwrap()
        .unwrap();
    assert_eq!((detail.run.as_str(), detail.key.as_str()), ("r1", "scan=1"));
    let StoredPeaks::Stored(peaks) = &detail.peaks else {
        panic!("peaks were attached: {:?}", detail.peaks);
    };
    assert_eq!(peaks.mz, [148.06, 227.10]);
    let assigned = detail.assigned.as_ref().unwrap();
    assert_eq!(
        (assigned.form.sequence.as_str(), assigned.charge),
        ("PEPTIDEK", Some(2))
    );
    assert_eq!(detail.psms.len(), 1);
    assert_eq!(detail.psms[0].scores, [("hyperscore".to_string(), 30.0)]);

    let unattached = project
        .spectrum_detail(spectra[1].spectrum_id)
        .unwrap()
        .unwrap();
    assert_eq!(unattached.peaks, StoredPeaks::Missing);
    assert_eq!(project.spectrum_detail(-1).unwrap(), None);
}

#[test]
fn biosample_categories_can_be_set_and_cleared() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = project(dir.path());
    project.set_category("treated", Some(" knockout ")).unwrap();
    project.set_category("control", None).unwrap();
    let categories: Vec<Option<String>> = project
        .protein_list()
        .unwrap()
        .biosamples
        .into_iter()
        .map(|b| b.category)
        .collect();
    assert_eq!(categories, [None, Some("knockout".to_string())]);
    let error = project.set_category("missing", Some("x")).unwrap_err();
    assert_eq!(
        error.to_string(),
        "no biosample named 'missing' in this project (biosamples: control, treated)"
    );
}

#[test]
fn probability_thresholds_filter_on_peps_and_rows_carry_probabilities() {
    use tuffold_validate::{Level, ThresholdMode, Thresholds};
    let dir = tempfile::tempdir().unwrap();
    let mut project = project(dir.path());
    let list = project.protein_list().unwrap();
    assert_eq!(list.thresholds.mode, ThresholdMode::Fdr);
    let probability = list.rows[0].probability.unwrap();
    assert!((probability - 0.998).abs() < 1e-12, "1 - PEP 0.002");
    let peptide = project
        .group_detail(list.rows[0].group_id)
        .unwrap()
        .unwrap()
        .peptides[0]
        .peptide_id;
    let passing = |project: &tuffold_store::Project| -> Vec<bool> {
        let page = project.peptide_spectra(peptide, 10).unwrap();
        page.rows.iter().map(|s| s.passes).collect()
    };
    // 1% FDR: scan=2 (q 0.02) fails.
    assert_eq!(passing(&project), [true, false]);

    // 95% PSM probability passes PEPs up to 0.05: scan=2 (PEP 0.02) too.
    project
        .save_thresholds(&Thresholds::probability_default())
        .unwrap();
    assert_eq!(passing(&project), [true, true]);
    let list = project.protein_list().unwrap();
    assert_eq!(list.thresholds.mode, ThresholdMode::Probability);
    assert_eq!(list.rows.len(), 1, "PEP 0.002 passes 99%");

    // 99.9% protein probability allows PEP 0.001: the group's 0.002 fails.
    let strict = Thresholds::probability_default().with(Level::Protein, 0.999);
    project.save_thresholds(&strict).unwrap();
    assert!(project.protein_list().unwrap().rows.is_empty());
    let summary = project.validation_summary().unwrap().unwrap();
    assert_eq!(summary.settings.thresholds, strict);
    assert_eq!(summary.protein_groups.targets, 0);
}

#[test]
fn ms1_features_become_precursor_intensities_where_a_spectrum_passes() {
    use std::collections::HashMap;
    use tuffold_quant::ms1::Feature;
    let dir = tempfile::tempdir().unwrap();
    let mut project = project(dir.path());
    let conn = Connection::open(dir.path().join("browse.tuff")).unwrap();
    let peptide = |sequence: &str| {
        id(
            &conn,
            &format!("SELECT id FROM peptides WHERE sequence = '{sequence}'"),
        )
    };
    let (pep, agg, dec) = (
        peptide("PEPTIDEK"),
        peptide("AGSAMPLER"),
        peptide("KEDITPEP"),
    );
    let feature = |area: f64| Feature {
        area,
        apex_rt: 60.0,
        apex_intensity: area / 10.0,
        start_rt: 50.0,
        end_rt: 70.0,
        points: 8,
        apex_ppm: 0.0,
        apex_ratio: 0.5,
    };

    // r1 scan 1 is PEPTIDEK at 2+. By hand: residues P 97.05276 ×2,
    // E 129.04259 ×2, T 101.04768, I 113.08406, D 115.02694, K 128.09496
    // = 909.44434, + water 18.01056 = 927.45490 Da; (927.45490 + 2 ×
    // 1.00728) / 2 = 464.7347.
    let mut attacher = project.attach_spectra("r1").unwrap();
    let precursors = attacher.precursors(&SpectrumKey::Scan(1)).to_vec();
    assert_eq!(precursors.len(), 1);
    assert_eq!((precursors[0].peptide_id, precursors[0].charge), (pep, 2));
    assert!((precursors[0].mz - 464.7347).abs() < 1e-3);
    attacher
        .record_ms1(
            10,
            &HashMap::from([
                ((pep, 2), Some(feature(1000.0))),
                // A 3+ trace no passing spectrum is at: never counted.
                ((pep, 3), Some(feature(5000.0))),
            ]),
        )
        .unwrap();
    let report = attacher.finish().unwrap();
    assert_eq!(
        (
            report.ms1_scans,
            report.precursors_traced,
            report.precursors_measured
        ),
        (10, 2, 2)
    );

    // r2: AGSAMPLER measured; the decoy traced, no peak; PEPTIDEK has a peak
    // but no passing r2 spectrum names it, so it doesn't count there.
    let mut attacher = project.attach_spectra("r2").unwrap();
    attacher
        .record_ms1(
            12,
            &HashMap::from([
                ((agg, 2), Some(feature(500.0))),
                ((dec, 2), None),
                ((pep, 2), Some(feature(777.0))),
            ]),
        )
        .unwrap();
    attacher.finish().unwrap();
    // Tracing r2 again replaces its three features, and says so.
    let mut attacher = project.attach_spectra("r2").unwrap();
    attacher
        .record_ms1(
            12,
            &HashMap::from([
                ((agg, 2), Some(feature(500.0))),
                ((dec, 2), None),
                ((pep, 2), Some(feature(777.0))),
            ]),
        )
        .unwrap();
    assert_eq!(attacher.finish().unwrap().precursors_replaced, 3);
    // Attaching a file without MS1 scans changes nothing.
    let mut attacher = project.attach_spectra("r2").unwrap();
    attacher.record_ms1(0, &HashMap::new()).unwrap();
    attacher.finish().unwrap();

    let list = project.intensity_list().unwrap();
    assert!(list.has_ms1());
    let names: Vec<&str> = list.biosamples.iter().map(|b| b.name.as_str()).collect();
    assert_eq!(names, ["control", "treated"]);
    assert_eq!(list.rows.len(), 1);
    let found: Vec<(&str, u8, Vec<Option<f64>>)> = list.rows[0]
        .peptides
        .iter()
        .map(|p| (p.proforma.as_str(), p.charge, p.intensities.clone()))
        .collect();
    assert_eq!(
        found,
        [
            ("AGSAMPLER", 2, vec![None, Some(500.0)]),
            ("PEPTIDEK", 2, vec![Some(1000.0), None]),
        ]
    );
    let runs: Vec<(Option<usize>, usize, usize, usize)> = list
        .runs
        .iter()
        .map(|r| (r.ms1_scans, r.traced, r.measured, r.untraced))
        .collect();
    assert_eq!(runs, [(Some(10), 2, 2, 0), (Some(12), 3, 2, 0)]);
}

#[test]
fn the_report_lists_every_passing_peptide_and_spectrum() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let report = project.report().unwrap();
    // Both target peptides pass 1% (q 0.001 and 0.005); the decoy never
    // appears. Spectra: r1/1 (q 0.001) and r2/1 (q 0.005); r1/2 (q 0.02)
    // fails the 1% PSM FDR.
    let peptides: Vec<(&str, usize, Vec<usize>)> = report
        .peptides
        .iter()
        .map(|p| (p.proforma.as_str(), p.spectra, p.spectra_by_sample.clone()))
        .collect();
    assert_eq!(
        peptides,
        [("PEPTIDEK", 1, vec![1, 0]), ("AGSAMPLER", 1, vec![0, 1])]
    );
    assert_eq!(
        report.peptides[0].groups,
        [("sp|P1|ONE".to_string(), "unique".to_string())]
    );
    let spectra: Vec<(&str, &str, &str)> = report
        .spectra
        .iter()
        .map(|s| (s.run.as_str(), s.key.as_str(), s.proforma.as_str()))
        .collect();
    assert_eq!(
        spectra,
        [("r1", "scan=1", "PEPTIDEK"), ("r2", "scan=1", "AGSAMPLER")]
    );
    // PEPs equal q-values in the fixture.
    assert!((report.spectra[0].probability - 0.999).abs() < 1e-12);
    assert_eq!(report.spectra[0].rt_seconds, Some(60.0));
}

#[test]
fn rows_list_their_modifications_and_stars() {
    let dir = tempfile::tempdir().unwrap();
    let mut project = project(dir.path());
    let conn = Connection::open(dir.path().join("browse.tuff")).unwrap();
    // PEPTIDEK (unique to the group) gets a named and an unnamed delta.
    conn.execute_batch(
        "INSERT INTO modifications (peptide_id, ordinal, site_kind, position, mass_delta, name)
         SELECT id, 90, 'residue', 1, 15.9949, 'Oxidation' FROM peptides WHERE sequence = 'PEPTIDEK';
         INSERT INTO modifications (peptide_id, ordinal, site_kind, position, mass_delta, name)
         SELECT id, 91, 'residue', 2, 0.98402, NULL FROM peptides WHERE sequence = 'PEPTIDEK';",
    )
    .unwrap();
    let row = &project.protein_list().unwrap().rows[0];
    // Sorted: '+' before 'O'; the unnamed delta to three decimals.
    assert_eq!(row.modifications, ["+0.984", "Oxidation"]);
    assert!(!row.starred);

    project.set_starred("sp|P1|ONE", true).unwrap();
    project.set_starred("sp|P1|ONE", true).unwrap();
    assert!(project.protein_list().unwrap().rows[0].starred);
    project.set_starred("sp|P1|ONE", false).unwrap();
    assert!(!project.protein_list().unwrap().rows[0].starred);
}

#[test]
fn groups_sharing_peptides_are_similar_and_clustered() {
    let dir = tempfile::tempdir().unwrap();
    let project = project(dir.path());
    let conn = Connection::open(dir.path().join("browse.tuff")).unwrap();
    // A second passing target group, led by a new protein, holding
    // AGSAMPLER as shared; and a third, failing (q 0.5), holding PEPTIDEK.
    conn.execute_batch(
        "INSERT INTO proteins (accession, description, sequence, is_decoy)
             VALUES ('sp|P2|TWO', 'Protein two', NULL, 0), ('sp|P3|THREE', NULL, NULL, 0);
         INSERT INTO protein_groups (id, leader_protein_id, is_decoy, score, q_value, peptides, pep)
             SELECT 900, id, 0, 2.0, 0.001, 1, 0.002 FROM proteins WHERE accession = 'sp|P2|TWO';
         INSERT INTO protein_groups (id, leader_protein_id, is_decoy, score, q_value, peptides, pep)
             SELECT 901, id, 0, 0.1, 0.5, 1, 0.6 FROM proteins WHERE accession = 'sp|P3|THREE';
         INSERT INTO protein_group_members (group_id, protein_id, role)
             SELECT 900, id, 'leader' FROM proteins WHERE accession = 'sp|P2|TWO';
         INSERT INTO protein_group_members (group_id, protein_id, role)
             SELECT 901, id, 'leader' FROM proteins WHERE accession = 'sp|P3|THREE';
         INSERT INTO protein_group_peptides (group_id, peptide_id, role)
             SELECT 900, id, 'shared' FROM peptides WHERE sequence = 'AGSAMPLER';
         INSERT INTO protein_group_peptides (group_id, peptide_id, role)
             SELECT 901, id, 'shared' FROM peptides WHERE sequence = 'PEPTIDEK';",
    )
    .unwrap();
    let agsampler: i64 = conn
        .query_row(
            "SELECT id FROM peptides WHERE sequence = 'AGSAMPLER'",
            [],
            |r| r.get(0),
        )
        .unwrap();
    let first = project
        .protein_list()
        .unwrap()
        .rows
        .iter()
        .find(|row| row.accession == "sp|P1|ONE")
        .unwrap()
        .group_id;
    let detail = project.group_detail(first).unwrap().unwrap();
    // Only the passing group counts: P2 shares AGSAMPLER; P3 fails at 1%.
    assert_eq!(
        detail.similar,
        [tuffold_store::SimilarGroup {
            group_id: 900,
            accession: "sp|P2|TWO".into(),
            description: Some("Protein two".into()),
            shared: vec![agsampler],
            peptides: 1,
        }]
    );
    assert_eq!(detail.cluster_size, 2);
    let lone = project.group_detail(900).unwrap().unwrap();
    assert_eq!(lone.similar[0].group_id, first);
    assert_eq!(lone.cluster_size, 2);
}
