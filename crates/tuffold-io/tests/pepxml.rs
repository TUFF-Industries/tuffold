//! pepXML importer tests.
//!
//! The oracle is derived BY HAND from
//! fixtures/third_party/pyteomics/test.pep.xml (Apache-2.0, from pyteomics),
//! a deliberately mixed real-world sample: SEQUEST, Mascot, X!Tandem, and
//! Inspect hits under one SEQUEST search summary that declares AVERAGE mass
//! types. Undeclared modification deltas therefore use average residue
//! masses for every hit, as the file states (real X!Tandem or Inspect output
//! would declare monoisotopic masses and get monoisotopic deltas).
//! Masses below are worked out from the residue table in
//! tuffold-core/src/chemistry.rs; never regenerate them by running the code.

use std::io::Write as _;
use std::path::{Path, PathBuf};

use tuffold_core::{Evidence, ImportBatch, Psm, Score, SpectrumKey};
use tuffold_io::{XmlImportError, import_pepxml};

const PROTON: f64 = 1.007_276_466_621;
const RUN: &str = "pps_sl20060731_18mix_25ul_r1_1154456409";

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/third_party/pyteomics/test.pep.xml")
}

fn psm_with<'a>(batch: &'a ImportBatch, sequence: &str) -> &'a Psm {
    batch
        .psms
        .iter()
        .find(|psm| psm.peptide.sequence == sequence)
        .unwrap()
}

fn score(name: &str, value: f64) -> Score {
    Score {
        name: name.into(),
        accession: None,
        value,
    }
}

fn close(actual: Option<f64>, expected: f64) -> bool {
    actual.is_some_and(|value| (value - expected).abs() < 1e-9)
}

#[test]
fn imports_the_pyteomics_sample_by_hand_derived_oracle() {
    let batch = import_pepxml(&fixture()).unwrap();

    // One msms_run_summary: base_name's file name, raw_data ".mzXML".
    assert_eq!(batch.software.name, "SEQUEST");
    assert_eq!(batch.ms_runs.len(), 1);
    assert_eq!(batch.ms_runs[0].name, RUN);
    assert_eq!(
        batch.ms_runs[0].spectra_file.as_deref(),
        Some(&*format!("{RUN}.mzXML"))
    );
    assert_eq!(batch.biosamples[0].name, RUN);

    // Five spectrum_query elements keyed by start_scan, file order.
    let keys: Vec<&SpectrumKey> = batch.spectra.iter().map(|s| &s.key).collect();
    assert_eq!(
        keys,
        [
            &SpectrumKey::Scan(100),
            &SpectrumKey::Scan(40),
            &SpectrumKey::Scan(1366),
            &SpectrumKey::Scan(6862),
            &SpectrumKey::Scan(23)
        ]
    );
    // Scan 100: precursor_neutral_mass 860.3920 at 1+ -> 861.399276...
    assert!(close(batch.spectra[0].precursor_mz, 860.392 + PROTON));
    assert_eq!(batch.spectra[2].rt_seconds, Some(38.426123));
    assert_eq!(batch.spectra[2].source_id.as_deref(), Some("MASCOT"));
    // &quot; entities are unescaped in the Inspect title.
    assert_eq!(
        batch.spectra[4].source_id.as_deref(),
        Some("\"Cmpd 24, +MSn(1085.6886), 1.2 min.23.23.3\"")
    );

    // Six search_hit elements (the Inspect query has ranks 1 and 2).
    assert_eq!(batch.psms.len(), 6);

    let sequest = psm_with(&batch, "SLNGEWR");
    assert_eq!((sequest.rank, sequest.charge), (1, 1));
    assert!(close(sequest.calculated_mz, 860.892 + PROTON));
    assert!(close(sequest.experimental_mz, 860.392 + PROTON));
    assert_eq!(
        sequest.scores,
        [
            score("xcorr", 1.553),
            score("deltacn", 0.081),
            score("deltacnstar", 0.0),
            score("spscore", 894.0),
            score("sprank", 1.0),
            // PeptideProphet's 0.96 is a post-processor's verdict: not imported.
        ]
    );
    assert_eq!(
        sequest.evidence,
        [Evidence {
            accession: "sp|P00722|BGAL_ECOLI".into(),
            start: None,
            end: None,
            pre: Some('R'),
            post: Some('F'),
        }]
    );

    // X!Tandem: M10 at 147.0354 - average M 131.1926 = +15.8428 (average
    // masses declared). The N-term 1.0079 and C-term 17.0031 differ from the
    // average H 1.00794 / OH 17.00734 by under 0.005 and are dropped.
    let tandem = psm_with(&batch, "EVPLNTIIFMGR");
    assert_eq!(tandem.peptide.proforma(), "EVPLNTIIFM[+15.8428]GR");
    let accessions: Vec<&str> = tandem
        .evidence
        .iter()
        .map(|e| e.accession.as_str())
        .collect();
    assert_eq!(
        accessions,
        ["sp|P01008|ANT3_HUMAN", "tr|Q8TCE1|Q8TCE1_HUMAN"]
    );
    assert_eq!(tandem.scores.last(), Some(&score("expect", 0.012)));

    // Inspect rank 1: no declared match for S 166.99803 or C 160.03019 (the
    // SEQUEST C entry is 160.1388), so deltas use average residue masses:
    // S 166.99803 - 87.0782 = +79.9198, C 160.03019 - 103.1439 = +56.8863.
    let inspect = psm_with(&batch, "DQQFDSSSSMALEDCGEETNCQSDFK");
    assert_eq!(
        inspect.peptide.proforma(),
        "DQQFDS[+79.9198]S[+79.9198]SS[+79.9198]MALEDC[+56.8863]GEETNC[+56.8863]QSDFK"
    );
    // Rank 2: Y 243.02933 - 163.1760 = +79.8533; K 170.10596 - 128.1741 =
    // +41.9319; T 181.01368 - 101.1047 = +79.9090; C +56.8863.
    let rank2 = psm_with(&batch, "CENCNYPKEGTHTNQHETLHTSR");
    assert_eq!(rank2.rank, 2);
    assert_eq!(
        rank2.peptide.proforma(),
        "C[+56.8863]ENC[+56.8863]NY[+79.8533]PK[+41.9319]EGT[+79.9090]HT[+79.9090]NQHET[+79.9090]LHT[+79.9090]SR"
    );

    // Proteins in first-seen order; the three Inspect proteins are decoys by
    // their ":reversed" descriptions, so both Inspect PSMs are decoys.
    let proteins: Vec<(&str, bool)> = batch
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p.is_decoy))
        .collect();
    assert_eq!(
        proteins,
        [
            ("sp|P00722|BGAL_ECOLI", false),
            ("gi|3212198|gb|AAC22319.1|", false),
            ("IPI00200898", false),
            ("sp|P01008|ANT3_HUMAN", false),
            ("tr|Q8TCE1|Q8TCE1_HUMAN", false),
            ("BL_ORD_ID:125453", true),
            ("BL_ORD_ID:144314", true),
            ("BL_ORD_ID:154629", true),
        ]
    );
    let decoys = batch
        .psms
        .iter()
        .filter(|psm| batch.psm_is_decoy(psm))
        .count();
    assert_eq!(decoys, 2);
}

fn write(dir: &Path, name: &str, body: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, wrap(body)).unwrap();
    path
}

fn wrap(body: &str) -> String {
    format!(
        "<?xml version=\"1.0\"?>\n<msms_pipeline_analysis xmlns=\"http://regis-web.systemsbiology.net/pepXML\">\n\
         <msms_run_summary base_name=\"/data/runA\" raw_data=\".mzML\">\n\
         <search_summary search_engine=\"Comet\" search_engine_version=\"2024.01\">\n\
         <aminoacid_modification aminoacid=\"C\" massdiff=\"57.021464\" mass=\"160.030649\" variable=\"N\"/>\n\
         </search_summary>\n{body}\n</msms_run_summary>\n</msms_pipeline_analysis>\n"
    )
}

#[test]
fn declared_and_explicit_modification_deltas_win() {
    let dir = tempfile::tempdir().unwrap();
    let path = write(
        dir.path(),
        "a.pep.xml",
        r#"<spectrum_query spectrum="runA.00007.00007.2" start_scan="7" assumed_charge="2" precursor_neutral_mass="1000.5">
<search_result><search_hit hit_rank="1" peptide="MCK" protein="DECOY_P9" peptide_prev_aa="-" peptide_next_aa="A">
<modification_info><mod_aminoacid_mass position="2" mass="160.030649"/><mod_aminoacid_mass position="1" mass="147.0" variable="15.994915"/></modification_info>
<search_score name="xcorr" value="2.5"/><search_score name="flag" value="n/a"/>
</search_hit></search_result></spectrum_query>"#,
    );

    let batch = import_pepxml(&path).unwrap();

    assert_eq!(batch.software.name, "Comet");
    assert_eq!(batch.software.version.as_deref(), Some("2024.01"));
    assert_eq!(batch.ms_runs[0].spectra_file.as_deref(), Some("runA.mzML"));
    let psm = &batch.psms[0];
    // M1 uses its explicit variable delta; C2 matches the declared table.
    assert_eq!(psm.peptide.proforma(), "M[+15.9949]C[+57.0215]K");
    // "-" flanks mean protein terminus: no residue. Text scores are skipped.
    assert_eq!(
        (psm.evidence[0].pre, psm.evidence[0].post),
        (None, Some('A'))
    );
    assert_eq!(psm.scores, [score("xcorr", 2.5)]);
    assert!(batch.psm_is_decoy(psm));
}

#[test]
fn gzip_files_are_read_transparently() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("a.pep.xml.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(&path).unwrap(),
        flate2::Compression::default(),
    );
    encoder
        .write_all(
            wrap(
                r#"<spectrum_query spectrum="s" index="3" assumed_charge="2">
<search_result><search_hit hit_rank="1" peptide="PEPTIDE" protein="P1"/></search_result></spectrum_query>"#,
            )
            .as_bytes(),
        )
        .unwrap();
    encoder.finish().unwrap();

    let batch = import_pepxml(&path).unwrap();

    assert_eq!(batch.psms.len(), 1);
    // No start_scan: the title keys it. The query index (a 1-based ordinal,
    // not a spectrum index) never does.
    assert_eq!(batch.spectra[0].key, SpectrumKey::Native("s".into()));
}

#[test]
fn missing_required_attributes_name_element_and_line() {
    let dir = tempfile::tempdir().unwrap();
    let path = write(
        dir.path(),
        "a.pep.xml",
        "<spectrum_query spectrum=\"s\" start_scan=\"1\">\n</spectrum_query>",
    );

    let error = import_pepxml(&path).unwrap_err();

    let XmlImportError::Invalid {
        line,
        element,
        problem,
        ..
    } = &error
    else {
        panic!("unexpected error: {error}");
    };
    // Line 1 declaration, 2 root, 3 run summary, 4-6 search summary, 7 query.
    assert_eq!((*line, element.as_str()), (7, "spectrum_query"));
    assert!(problem.contains("assumed_charge"), "{problem}");
}

#[test]
fn modification_positions_outside_the_peptide_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = write(
        dir.path(),
        "a.pep.xml",
        r#"<spectrum_query spectrum="s" start_scan="1" assumed_charge="2"><search_result>
<search_hit hit_rank="1" peptide="PEP" protein="P1"><modification_info><mod_aminoacid_mass position="9" mass="100.0"/></modification_info></search_hit>
</search_result></spectrum_query>"#,
    );
    let message = import_pepxml(&path).unwrap_err().to_string();
    assert!(
        message.contains("position 9 outside peptide PEP"),
        "{message}"
    );
}

#[test]
fn malformed_xml_reports_its_line() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("broken.pep.xml");
    std::fs::write(&path, "<?xml version=\"1.0\"?>\n<a>\n<b></c>\n</a>\n").unwrap();
    let error = import_pepxml(&path).unwrap_err();
    assert!(
        matches!(error, XmlImportError::Malformed { line: 3, .. }),
        "{error}"
    );
}

#[test]
fn same_scan_queries_share_a_spectrum_and_mods_follow_their_search_id() {
    // Two search summaries: search 1 declares C+57, search 2 declares C+58
    // (deliberately different, to see which one a hit uses). Scan 7 is
    // queried at charge 2 (search 2) and charge 3 (search 1).
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("multi.pep.xml");
    std::fs::write(
        &path,
        r#"<?xml version="1.0"?>
<msms_pipeline_analysis><msms_run_summary base_name="/data/runB.mzML" raw_data=".mzML">
<search_summary search_engine="Comet" search_id="1" precursor_mass_type="monoisotopic" fragment_mass_type="monoisotopic">
<aminoacid_modification aminoacid="C" massdiff="57.021464" mass="160.030649"/></search_summary>
<search_summary search_engine="X! Tandem" search_id="2" precursor_mass_type="monoisotopic" fragment_mass_type="monoisotopic">
<aminoacid_modification aminoacid="C" massdiff="58.0" mass="160.030649"/></search_summary>
<spectrum_query spectrum="runB.7.7.2" start_scan="7" assumed_charge="2" precursor_neutral_mass="998.0">
<search_result search_id="2"><search_hit hit_rank="1" peptide="CK" protein="P1"><modification_info><mod_aminoacid_mass position="1" mass="160.030649"/></modification_info></search_hit></search_result></spectrum_query>
<spectrum_query spectrum="runB.7.7.3" start_scan="7" assumed_charge="3" precursor_neutral_mass="1497.0" retention_time_sec="12.5">
<search_result search_id="1"><search_hit hit_rank="1" peptide="CK" protein="P1"><modification_info><mod_aminoacid_mass position="1" mass="160.030649"/></modification_info></search_hit></search_result></spectrum_query>
</msms_run_summary></msms_pipeline_analysis>
"#,
    )
    .unwrap();

    let batch = import_pepxml(&path).unwrap();

    // ".mzML" in base_name is not part of the run name.
    assert_eq!(batch.ms_runs[0].name, "runB");
    // One spectrum: charges 2 and 3 conflict, so the charge is unknown; the
    // second query's RT fills in. m/z agrees: (998 + 2p)/2 = (1497 + 3p)/3.
    assert_eq!(batch.spectra.len(), 1);
    let spectrum = &batch.spectra[0];
    assert_eq!(
        (spectrum.precursor_charge, spectrum.rt_seconds),
        (None, Some(12.5))
    );
    assert!((spectrum.precursor_mz.unwrap() - (499.0 + PROTON)).abs() < 1e-9);
    // Each hit uses its own search summary's declared C delta.
    let keys: Vec<String> = batch
        .psms
        .iter()
        .map(|psm| psm.peptide.proforma())
        .collect();
    assert_eq!(keys, ["C[+58.0000]K", "C[+57.0215]K"]);
}
