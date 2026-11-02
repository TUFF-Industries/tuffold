//! mzIdentML importer tests.
//!
//! The oracle is derived BY HAND from fixtures/third_party/pyteomics/test.mzid
//! (Apache-2.0, from pyteomics; the PSI "MPC_use_case_handcrafted" example):
//! two SpectrumIdentificationLists (Sequest, Mascot) over one SpectraData,
//! 7 DBSequences, 19 Peptides, 22 SpectrumIdentificationItems, each with one
//! PeptideEvidenceRef. Never regenerate these values by running the code.

use std::path::{Path, PathBuf};

use tuffold_core::{Evidence, ImportBatch, Psm, Score, SpectrumKey};
use tuffold_io::{XmlImportError, import_mzid};

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/third_party/pyteomics/test.mzid")
}

fn psm_for<'a>(batch: &'a ImportBatch, sequence: &str, calculated_mz: f64) -> &'a Psm {
    batch
        .psms
        .iter()
        .find(|psm| psm.peptide.sequence == sequence && psm.calculated_mz == Some(calculated_mz))
        .unwrap()
}

#[test]
fn imports_the_psi_example_by_hand_derived_oracle() {
    let batch = import_mzid(&fixture()).unwrap();

    // SEQUEST_SW and Mascot_SW, named by their SoftwareName cvParams.
    assert_eq!(batch.software.name, "Sequest + Mascot");
    assert_eq!(batch.software.version, None);

    // SpectraData location ".../Separation_1D_LC/Fraction_X".
    assert_eq!(batch.ms_runs.len(), 1);
    assert_eq!(batch.ms_runs[0].name, "Fraction_X");

    // Sequest results on databasekey 1,2,3,10,11,12,13,15,20; Mascot adds
    // 4,6,35,36,40 and repeats 2,3,11,12 -> 14 distinct spectra.
    assert_eq!(batch.spectra.len(), 14);
    assert_eq!(
        batch.spectra[0].key,
        SpectrumKey::Native("databasekey=1".into())
    );
    assert_eq!(batch.psms.len(), 22);
    let on_spectrum_2 = batch
        .psms
        .iter()
        .filter(|psm| psm.spectrum == SpectrumKey::Native("databasekey=2".into()))
        .count();
    assert_eq!(
        on_spectrum_2, 2,
        "Sequest and Mascot hits share databasekey=2"
    );

    // SEQ_spec1_pep1, every field.
    let first = psm_for(&batch, "AGTQIENIDEDFR", 1507.695);
    assert_eq!(
        (first.rank, first.charge, first.pass_threshold),
        (1, 1, Some(true))
    );
    assert_eq!(first.experimental_mz, Some(1507.696));
    assert_eq!(
        first.scores,
        [
            Score {
                name: "ProteinScape:IntensityCoverage".into(),
                accession: Some("MS:1001505".into()),
                value: 0.391_954_560_380_971_8,
            },
            Score {
                name: "ProteinScape:SequestMetaScore".into(),
                accession: Some("MS:1001506".into()),
                value: 7.594_885_189_034_25,
            },
        ]
    );
    assert_eq!(
        first.evidence,
        [Evidence {
            accession: "IPI00013808.1".into(),
            start: Some(67),
            end: Some(79),
            pre: Some('K'),
            post: Some('D'),
        }]
    );

    // prot1_pep4: Oxidation (UNIMOD:35) at location 1 with no mass stated,
    // so the delta comes from the Unimod table.
    let oxidized = psm_for(&batch, "MLDAEDIVNTARPDEK", 1832.862115);
    assert_eq!(oxidized.peptide.proforma(), "M[+15.9949]LDAEDIVNTARPDEK");
    assert_eq!(
        oxidized.peptide.modifications[0].accession.as_deref(),
        Some("UNIMOD:35")
    );
    assert_eq!(
        oxidized.peptide.modifications[0].name.as_deref(),
        Some("Oxidation")
    );

    // Mas_spec6 has two rank-1 items (a tie); both are kept.
    let spec6: Vec<u32> = batch
        .psms
        .iter()
        .filter(|psm| psm.spectrum == SpectrumKey::Native("databasekey=6".into()))
        .map(|psm| psm.rank)
        .collect();
    assert_eq!(spec6, [1, 1]);

    // 19 distinct peptide forms.
    let mut forms: Vec<String> = batch
        .psms
        .iter()
        .map(|psm| psm.peptide.proforma())
        .collect();
    forms.sort();
    forms.dedup();
    assert_eq!(forms.len(), 19);

    // isDecoy="true" on every evidence of prot4_SHD and prot7_SHD. Their
    // "SHD" accessions carry no decoy label, so the explicit flag decides.
    let decoys: Vec<&str> = batch
        .proteins
        .iter()
        .filter(|protein| protein.is_decoy)
        .map(|protein| protein.accession.as_str())
        .collect();
    assert_eq!(decoys, ["SHD00382470.3", "SHD00644576.1"]);
    assert_eq!(batch.proteins.len(), 7);
    // SEQ_spec15, SEQ_spec20 (prot4); Mas_spec6_pep2, spec36b1, spec40 (prot7).
    let decoy_psms = batch
        .psms
        .iter()
        .filter(|psm| batch.psm_is_decoy(psm))
        .count();
    assert_eq!(decoy_psms, 5);

    // Descriptions drop the leading '>' some exporters keep; sequences load.
    let keratin = batch
        .proteins
        .iter()
        .find(|p| p.accession == "IPI00554648.1")
        .unwrap();
    assert_eq!(
        keratin.description.as_deref(),
        Some("IPI:IPI00554648.1|SWISS-PROT:P05787 Tax_Id=9606 Keratin, type II cytoskeletal 8")
    );
    assert!(
        keratin
            .sequence
            .as_deref()
            .unwrap()
            .starts_with("SIRVTQKSYKVSTSGPRAFSS")
    );
}

fn mzid(dir: &Path, sequence_collection: &str, results: &str) -> PathBuf {
    let path = dir.join("search.mzid");
    std::fs::write(
        &path,
        format!(
            r#"<?xml version="1.0" encoding="UTF-8"?>
<MzIdentML xmlns="http://psidev.info/psi/pi/mzIdentML/1.2" version="1.2.0">
<AnalysisSoftwareList><AnalysisSoftware id="sw" version="v2024"><SoftwareName><cvParam accession="MS:1002048" name="MS-GF+" cvRef="PSI-MS"/></SoftwareName></AnalysisSoftware></AnalysisSoftwareList>
<SequenceCollection>{sequence_collection}</SequenceCollection>
<AnalysisCollection><SpectrumIdentification id="si" spectrumIdentificationProtocol_ref="p" spectrumIdentificationList_ref="l"/></AnalysisCollection>
<AnalysisProtocolCollection><SpectrumIdentificationProtocol id="p" analysisSoftware_ref="sw"/></AnalysisProtocolCollection>
<DataCollection><Inputs><SpectraData id="sd" location="/data/HeLa_01.mzML"/></Inputs>
<AnalysisData><SpectrumIdentificationList id="l">{results}</SpectrumIdentificationList></AnalysisData></DataCollection>
</MzIdentML>
"#
        ),
    )
    .unwrap();
    path
}

const PROTEINS: &str = r#"<DBSequence id="d1" accession="sp|P1|ONE_HUMAN"/><DBSequence id="d2" accession="XXX_sp|P1|ONE_HUMAN"/>"#;

#[test]
fn single_engine_versions_units_and_decoy_conventions() {
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        &format!(
            r#"{PROTEINS}<Peptide id="pep"><PeptideSequence>PEPTIDEK</PeptideSequence>
<Modification location="0" monoisotopicMassDelta="42.010565"><cvParam accession="UNIMOD:1" name="Acetyl" cvRef="UNIMOD"/></Modification></Peptide>
<PeptideEvidence id="e1" peptide_ref="pep" dBSequence_ref="d1" pre="-" post="A"/>
<PeptideEvidence id="e2" peptide_ref="pep" dBSequence_ref="d2"/>"#
        ),
        r#"<SpectrumIdentificationResult id="r" spectrumID="controllerType=0 controllerNumber=1 scan=4021" spectraData_ref="sd">
<SpectrumIdentificationItem id="i" peptide_ref="pep" rank="1" chargeState="2" experimentalMassToCharge="500.25" passThreshold="false">
<PeptideEvidenceRef peptideEvidence_ref="e1"/><PeptideEvidenceRef peptideEvidence_ref="e2"/>
<cvParam accession="MS:1002053" name="MS-GF:EValue" cvRef="PSI-MS" value="1.2e-5"/>
<cvParam accession="MS:1001363" name="peptide unique to one protein" cvRef="PSI-MS"/>
<userParam name="IsotopeError" value="0"/>
</SpectrumIdentificationItem>
<cvParam accession="MS:1000016" name="scan start time" cvRef="PSI-MS" value="12.5" unitAccession="UO:0000031" unitName="minute"/>
</SpectrumIdentificationResult>"#,
    );

    let batch = import_mzid(&path).unwrap();

    assert_eq!(batch.software.name, "MS-GF+");
    assert_eq!(batch.software.version.as_deref(), Some("v2024"));
    assert_eq!(batch.ms_runs[0].name, "HeLa_01");
    let spectrum = &batch.spectra[0];
    assert_eq!(spectrum.key, SpectrumKey::Scan(4021));
    assert_eq!(spectrum.rt_seconds, Some(750.0));
    assert_eq!(
        (spectrum.precursor_mz, spectrum.precursor_charge),
        (Some(500.25), Some(2))
    );
    let psm = &batch.psms[0];
    assert_eq!(psm.pass_threshold, Some(false));
    assert_eq!(psm.peptide.proforma(), "[+42.0106]-PEPTIDEK");
    // Valueless flags are not scores; userParams are, without accession.
    let names: Vec<&str> = psm.scores.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, ["MS-GF:EValue", "IsotopeError"]);
    assert_eq!(
        (psm.evidence[0].pre, psm.evidence[0].post),
        (None, Some('A'))
    );
    // No isDecoy attributes: the XXX_ accession label decides.
    let decoys: Vec<bool> = batch.proteins.iter().map(|p| p.is_decoy).collect();
    assert_eq!(decoys, [false, true]);
    assert!(!batch.psm_is_decoy(psm), "shared with a target protein");
}

#[test]
fn modifications_without_any_mass_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        r#"<DBSequence id="d1" accession="P1"/><Peptide id="pep"><PeptideSequence>PEPTIDEK</PeptideSequence>
<Modification location="3"><cvParam accession="MS:1001460" name="unknown modification" cvRef="PSI-MS"/></Modification></Peptide>
<PeptideEvidence id="e1" peptide_ref="pep" dBSequence_ref="d1"/>"#,
        r#"<SpectrumIdentificationResult id="r" spectrumID="index=1" spectraData_ref="sd">
<SpectrumIdentificationItem id="i" peptide_ref="pep" rank="1" chargeState="2"><PeptideEvidenceRef peptideEvidence_ref="e1"/></SpectrumIdentificationItem>
</SpectrumIdentificationResult>"#,
    );

    let error = import_mzid(&path).unwrap_err();

    let XmlImportError::Invalid {
        element, problem, ..
    } = &error
    else {
        panic!("unexpected error: {error}");
    };
    assert_eq!(element, "SpectrumIdentificationItem");
    assert!(
        problem.contains("MS:1001460") && problem.contains("no monoisotopic mass delta"),
        "{problem}"
    );
}

#[test]
fn evidence_for_another_peptide_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        r#"<DBSequence id="d1" accession="P1"/>
<Peptide id="a"><PeptideSequence>AAAK</PeptideSequence></Peptide><Peptide id="b"><PeptideSequence>BBBK</PeptideSequence></Peptide>
<PeptideEvidence id="eb" peptide_ref="b" dBSequence_ref="d1"/>"#,
        r#"<SpectrumIdentificationResult id="r" spectrumID="index=1" spectraData_ref="sd">
<SpectrumIdentificationItem id="i" peptide_ref="a" rank="1" chargeState="2"><PeptideEvidenceRef peptideEvidence_ref="eb"/></SpectrumIdentificationItem>
</SpectrumIdentificationResult>"#,
    );
    let message = import_mzid(&path).unwrap_err().to_string();
    assert!(
        message.contains("belongs to Peptide 'b', not 'a'"),
        "{message}"
    );
}

#[test]
fn unknown_references_name_what_is_missing() {
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        "",
        r#"<SpectrumIdentificationResult id="r" spectrumID="index=1" spectraData_ref="nowhere"/>"#,
    );
    let message = import_mzid(&path).unwrap_err().to_string();
    assert!(
        message.contains("unknown SpectraData 'nowhere'"),
        "{message}"
    );
}

#[test]
fn substitutions_are_reported_as_unsupported() {
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        r#"<Peptide id="p"><PeptideSequence>PEPK</PeptideSequence><SubstitutionModification originalResidue="E" replacementResidue="Q" location="2"/></Peptide>"#,
        "",
    );
    let message = import_mzid(&path).unwrap_err().to_string();
    assert!(message.contains("not supported yet"), "{message}");
}

#[test]
fn scan_number_params_key_spectra_like_other_importers() {
    // MS-GF+ style: spectrumID is a 0-based index; the scan is a cvParam.
    // Two results on the same scan merge; the later one supplies the RT the
    // first lacked.
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        r#"<DBSequence id="d1" accession="P1"/><Peptide id="pep"><PeptideSequence>PEPTIDEK</PeptideSequence></Peptide>
<PeptideEvidence id="e1" peptide_ref="pep" dBSequence_ref="d1" start="10" end="17"/>
<PeptideEvidence id="e2" peptide_ref="pep" dBSequence_ref="d1" start="40" end="47"/>"#,
        r#"<SpectrumIdentificationResult id="r1" spectrumID="index=1406" spectraData_ref="sd">
<SpectrumIdentificationItem id="i1" peptide_ref="pep" rank="1" chargeState="2"><PeptideEvidenceRef peptideEvidence_ref="e1"/><PeptideEvidenceRef peptideEvidence_ref="e2"/></SpectrumIdentificationItem>
<cvParam accession="MS:1001115" name="scan number(s)" cvRef="PSI-MS" value="4021"/>
</SpectrumIdentificationResult>
<SpectrumIdentificationResult id="r2" spectrumID="index=1406" spectraData_ref="sd">
<SpectrumIdentificationItem id="i2" peptide_ref="pep" rank="2" chargeState="2"><PeptideEvidenceRef peptideEvidence_ref="e1"/></SpectrumIdentificationItem>
<cvParam accession="MS:1001115" name="scan number(s)" cvRef="PSI-MS" value="4021"/>
<cvParam accession="MS:1000016" name="scan start time" cvRef="PSI-MS" value="30" unitName="second"/>
</SpectrumIdentificationResult>"#,
    );

    let batch = import_mzid(&path).unwrap();

    assert_eq!(batch.spectra.len(), 1);
    assert_eq!(batch.spectra[0].key, SpectrumKey::Scan(4021));
    assert_eq!(batch.spectra[0].rt_seconds, Some(30.0));
    assert!(
        batch
            .psms
            .iter()
            .all(|psm| psm.spectrum == SpectrumKey::Scan(4021))
    );
    // PEPTIDEK occurs twice in P1: two loci, both kept.
    let starts: Vec<Option<u32>> = batch.psms[0].evidence.iter().map(|e| e.start).collect();
    assert_eq!(starts, [Some(10), Some(40)]);
}

#[test]
fn average_only_deltas_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let path = mzid(
        dir.path(),
        r#"<DBSequence id="d1" accession="P1"/><Peptide id="pep"><PeptideSequence>PEPTIDEK</PeptideSequence>
<Modification location="2" avgMassDelta="80.0"><cvParam accession="MS:1001460" name="unknown modification" cvRef="PSI-MS"/></Modification></Peptide>
<PeptideEvidence id="e1" peptide_ref="pep" dBSequence_ref="d1"/>"#,
        r#"<SpectrumIdentificationResult id="r" spectrumID="index=1" spectraData_ref="sd">
<SpectrumIdentificationItem id="i" peptide_ref="pep" rank="1" chargeState="2"><PeptideEvidenceRef peptideEvidence_ref="e1"/></SpectrumIdentificationItem>
</SpectrumIdentificationResult>"#,
    );
    let message = import_mzid(&path).unwrap_err().to_string();
    assert!(message.contains("only an average delta, 80"), "{message}");
}
