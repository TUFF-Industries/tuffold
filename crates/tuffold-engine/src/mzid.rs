//! mzIdentML 1.2 export of a validated project: Tuffold's results, for
//! PRIDE submissions and other tools.
//!
//! - One `SpectrumIdentificationResult` per validated spectrum, target and
//!   decoy, holding Tuffold's rank-1 peptide for it (the one the engines'
//!   votes chose), with `passThreshold` at the saved PSM threshold, the
//!   PSM-level q-value (MS:1002354) and local FDR (= PEP, MS:1002351), and
//!   Tuffold's combined score as a userParam. Engine scores stay in the
//!   engines' own files.
//! - Peptides carry their modifications (`location` 0 is the N-terminus,
//!   length + 1 the C-terminus) with the Unimod accession when known, else
//!   "unknown modification" with the mass delta; `PeptideEvidence` marks
//!   decoys (`isDecoy`).
//! - The `ProteinDetectionList` holds every protein group: leader first
//!   (`group representative`, `leading protein`), the q-value (MS:1002373)
//!   and local FDR (MS:1002370) of the group, and `passThreshold` at the
//!   saved protein threshold and peptide minimum.
//! - The database is named but not located (Tuffold keeps only the proteins
//!   it saw): `SearchDatabase/@location` is `not-recorded`, a valid URI.
//! - `SpectraData/@location` is the attached spectra file, percent-encoded,
//!   or the run's name when none was attached; `FileFormat` follows the
//!   extension (the generic "mass spectrometer file format" otherwise).
//!   Checked against the mzIdentML 1.2.0 XSD.

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::io::{self, Write};

use tuffold_core::{ModSite, Modification};
use tuffold_store::{Exchange, ExchangeSpectrum, Project, StoreError};
use tuffold_validate::ThresholdMode;

#[derive(Debug, thiserror::Error)]
pub enum MzidError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error(transparent)]
    Io(#[from] io::Error),
    #[error("this project has no current validation results; validate it first")]
    NotValidated,
}

/// What the file holds. Passing counts are targets, as everywhere in
/// Tuffold; passing decoys are in the file too, marked as passing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MzidReport {
    pub spectra: usize,
    pub passing_spectra: usize,
    pub peptides: usize,
    pub protein_groups: usize,
    pub passing_groups: usize,
}

/// Write the project's validated results as mzIdentML 1.2.
pub fn write_mzid(project: &Project, out: &mut impl Write) -> Result<MzidReport, MzidError> {
    let current = project
        .validation_summary()?
        .is_some_and(|summary| summary.validated_at.is_some());
    if !current {
        return Err(MzidError::NotValidated);
    }
    let data = project.exchange()?;
    let mut writer = Writer {
        out,
        line: String::new(),
        native_runs: native_runs(&data.spectra),
    };
    writer.document(&data)?;
    Ok(MzidReport {
        spectra: data.spectra.len(),
        passing_spectra: data
            .spectra
            .iter()
            .filter(|s| s.passes && !s.is_decoy)
            .count(),
        peptides: data.peptides.len(),
        protein_groups: data.groups.len(),
        passing_groups: data
            .groups
            .iter()
            .filter(|g| g.passes && !g.is_decoy)
            .count(),
    })
}

struct Writer<'a, W: Write> {
    out: &'a mut W,
    line: String,
    /// Runs whose spectra all have the spectra file's own id.
    native_runs: HashSet<i64>,
}

impl<W: Write> Writer<'_, W> {
    /// Write the buffered line(s).
    fn flush_line(&mut self) -> io::Result<()> {
        self.out.write_all(self.line.as_bytes())?;
        self.line.clear();
        Ok(())
    }

    fn put(&mut self, text: &str) -> io::Result<()> {
        self.line.push_str(text);
        self.line.push('\n');
        if self.line.len() > 1 << 16 {
            self.flush_line()?;
        }
        Ok(())
    }

    fn document(&mut self, data: &Exchange) -> io::Result<()> {
        let version = env!("CARGO_PKG_VERSION");
        self.put(r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
        self.put(
            r#"<MzIdentML id="tuffold" version="1.2.0" xmlns="http://psidev.info/psi/pi/mzIdentML/1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psidev.info/psi/pi/mzIdentML/1.2 https://www.psidev.info/files/mzIdentML1.2.0.xsd">"#,
        )?;
        self.put("  <cvList>")?;
        self.put(r#"    <cv id="PSI-MS" fullName="PSI-MS" uri="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>"#)?;
        self.put(
            r#"    <cv id="UNIMOD" fullName="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo"/>"#,
        )?;
        self.put(r#"    <cv id="UO" fullName="UNIT-ONTOLOGY" uri="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>"#)?;
        self.put("  </cvList>")?;
        self.put("  <AnalysisSoftwareList>")?;
        self.put(&format!(
            r#"    <AnalysisSoftware id="tuffold" name="Tuffold" version="{version}" uri="https://github.com/TUFF-Industries/tuffold">"#
        ))?;
        self.put(r#"      <SoftwareName><userParam name="Tuffold"/></SoftwareName>"#)?;
        self.put("    </AnalysisSoftware>")?;
        self.put("  </AnalysisSoftwareList>")?;
        self.sequences(data)?;
        self.put("  <AnalysisCollection>")?;
        self.put(r#"    <SpectrumIdentification id="SI" spectrumIdentificationProtocol_ref="SIP" spectrumIdentificationList_ref="SIL">"#)?;
        for run in &data.runs {
            self.put(&format!(
                r#"      <InputSpectra spectraData_ref="SD_{}"/>"#,
                run.id
            ))?;
        }
        self.put(r#"      <SearchDatabaseRef searchDatabase_ref="SDB"/>"#)?;
        self.put("    </SpectrumIdentification>")?;
        self.put(r#"    <ProteinDetection id="PD" proteinDetectionProtocol_ref="PDP" proteinDetectionList_ref="PDL">"#)?;
        self.put(r#"      <InputSpectrumIdentifications spectrumIdentificationList_ref="SIL"/>"#)?;
        self.put("    </ProteinDetection>")?;
        self.put("  </AnalysisCollection>")?;
        self.protocols(data)?;
        self.put("  <DataCollection>")?;
        self.inputs(data)?;
        self.put("    <AnalysisData>")?;
        self.spectra(data)?;
        self.groups(data)?;
        self.put("    </AnalysisData>")?;
        self.put("  </DataCollection>")?;
        self.put("</MzIdentML>")?;
        self.flush_line()
    }

    fn sequences(&mut self, data: &Exchange) -> io::Result<()> {
        self.put("  <SequenceCollection>")?;
        for protein in &data.proteins {
            self.put(&format!(
                r#"    <DBSequence id="DBSeq_{}" accession="{}" searchDatabase_ref="SDB"{}>"#,
                protein.id,
                escape(&protein.accession),
                protein
                    .sequence
                    .as_ref()
                    .map_or_else(String::new, |s| format!(
                        r#" length="{}""#,
                        s.chars().count()
                    ))
            ))?;
            if let Some(sequence) = &protein.sequence {
                self.put(&format!("      <Seq>{}</Seq>", escape(sequence)))?;
            }
            if let Some(description) = &protein.description {
                self.put(&format!(
                    r#"      <cvParam cvRef="PSI-MS" accession="MS:1001088" name="protein description" value="{}"/>"#,
                    escape(description)
                ))?;
            }
            self.put("    </DBSequence>")?;
        }
        for peptide in &data.peptides {
            let form = &peptide.form;
            self.put(&format!(r#"    <Peptide id="Pep_{}">"#, peptide.id))?;
            self.put(&format!(
                "      <PeptideSequence>{}</PeptideSequence>",
                escape(&form.sequence)
            ))?;
            let length = form.sequence.chars().count();
            for modification in &form.modifications {
                self.put(&modification_xml(modification, length))?;
            }
            self.put("    </Peptide>")?;
        }
        let decoys: HashMap<i64, bool> = data.proteins.iter().map(|p| (p.id, p.is_decoy)).collect();
        for (index, evidence) in data.evidence.iter().enumerate() {
            let mut attributes = format!(
                r#"id="{}" peptide_ref="Pep_{}" dBSequence_ref="DBSeq_{}""#,
                evidence_id(index),
                evidence.peptide_id,
                evidence.protein_id
            );
            if let (Some(start), Some(end)) = (evidence.start, evidence.end) {
                let _ = write!(attributes, r#" start="{start}" end="{end}""#);
            }
            for (name, residue) in [("pre", &evidence.pre), ("post", &evidence.post)] {
                if let Some(residue) = residue.as_deref().filter(|r| r.chars().count() == 1) {
                    let _ = write!(attributes, r#" {name}="{}""#, escape(residue));
                }
            }
            let decoy = decoys.get(&evidence.protein_id).copied().unwrap_or(false);
            self.put(&format!(
                r#"    <PeptideEvidence {attributes} isDecoy="{decoy}"/>"#
            ))?;
        }
        self.put("  </SequenceCollection>")
    }

    fn protocols(&mut self, data: &Exchange) -> io::Result<()> {
        let thresholds = &data.thresholds;
        let (psm, protein) = (thresholds.psm, thresholds.protein);
        let (psm_threshold, protein_threshold) = match thresholds.mode {
            ThresholdMode::Fdr => (
                format!(
                    r#"<cvParam cvRef="PSI-MS" accession="MS:1002350" name="PSM-level global FDR" value="{psm}"/>"#
                ),
                format!(
                    r#"<cvParam cvRef="PSI-MS" accession="MS:1002369" name="protein group-level global FDR" value="{protein}"/>"#
                ),
            ),
            ThresholdMode::Probability => (
                format!(
                    r#"<userParam name="Tuffold:PSM probability (1 - PEP) at least" value="{psm}"/>"#
                ),
                format!(
                    r#"<userParam name="Tuffold:protein group probability (1 - PEP) at least" value="{protein}"/>"#
                ),
            ),
        };
        self.put("  <AnalysisProtocolCollection>")?;
        self.put(
            r#"    <SpectrumIdentificationProtocol id="SIP" analysisSoftware_ref="tuffold">"#,
        )?;
        self.put(r#"      <SearchType><cvParam cvRef="PSI-MS" accession="MS:1001083" name="ms-ms search"/></SearchType>"#)?;
        self.put(&format!("      <Threshold>{psm_threshold}</Threshold>"))?;
        self.put("    </SpectrumIdentificationProtocol>")?;
        self.put(r#"    <ProteinDetectionProtocol id="PDP" analysisSoftware_ref="tuffold">"#)?;
        self.put(&format!(
            "      <AnalysisParams><userParam name=\"Tuffold:minimum unique or razor peptides\" value=\"{}\"/></AnalysisParams>",
            thresholds.min_peptides
        ))?;
        self.put(&format!("      <Threshold>{protein_threshold}</Threshold>"))?;
        self.put("    </ProteinDetectionProtocol>")?;
        self.put("  </AnalysisProtocolCollection>")
    }

    fn inputs(&mut self, data: &Exchange) -> io::Result<()> {
        self.put("    <Inputs>")?;
        // The schema wants a URI; Tuffold doesn't keep the FASTA's path.
        self.put(r#"      <SearchDatabase id="SDB" location="not-recorded">"#)?;
        self.put(r#"        <FileFormat><cvParam cvRef="PSI-MS" accession="MS:1001348" name="FASTA format"/></FileFormat>"#)?;
        self.put(r#"        <DatabaseName><userParam name="protein database (path not recorded; the proteins the searches matched are in SequenceCollection)"/></DatabaseName>"#)?;
        self.put("      </SearchDatabase>")?;
        let first_key: HashMap<i64, &str> = data
            .spectra
            .iter()
            .rev()
            .map(|s| (s.run_id, spectrum_id(&self.native_runs, s)))
            .collect();
        for run in &data.runs {
            // The spectra file when one was attached, else the run's name
            // as the results gave it (never an invented file).
            let location = run.spectra_file.as_deref().unwrap_or(&run.name);
            self.put(&format!(
                r#"      <SpectraData id="SD_{}" location="{}" name="{}">"#,
                run.id,
                escape(&uri(location)),
                escape(&run.name)
            ))?;
            let (format_accession, format_name) = file_format(run.spectra_file.as_deref());
            self.put(&format!(
                r#"        <FileFormat><cvParam cvRef="PSI-MS" accession="{format_accession}" name="{format_name}"/></FileFormat>"#
            ))?;
            let (accession, name) = id_format(first_key.get(&run.id).copied().unwrap_or(""));
            self.put(&format!(
                r#"        <SpectrumIDFormat><cvParam cvRef="PSI-MS" accession="{accession}" name="{name}"/></SpectrumIDFormat>"#
            ))?;
            self.put("      </SpectraData>")?;
        }
        self.put("    </Inputs>")
    }

    fn spectra(&mut self, data: &Exchange) -> io::Result<()> {
        // Evidence ids for each peptide, in file order.
        let mut evidence: HashMap<i64, Vec<usize>> = HashMap::new();
        for (index, item) in data.evidence.iter().enumerate() {
            evidence.entry(item.peptide_id).or_default().push(index);
        }
        self.put(r#"      <SpectrumIdentificationList id="SIL">"#)?;
        for spectrum in &data.spectra {
            self.put(&format!(
                r#"        <SpectrumIdentificationResult id="SIR_{}" spectrumID="{}" spectraData_ref="SD_{}">"#,
                spectrum.id,
                escape(spectrum_id(&self.native_runs, spectrum)),
                spectrum.run_id
            ))?;
            let mut item = format!(
                r#"          <SpectrumIdentificationItem id="SII_{}" chargeState="{}" experimentalMassToCharge="{}""#,
                spectrum.id,
                spectrum.charge.unwrap_or(0),
                spectrum.experimental_mz.unwrap_or(0.0)
            );
            if let Some(calculated) = spectrum.calculated_mz {
                let _ = write!(item, r#" calculatedMassToCharge="{calculated}""#);
            }
            let _ = write!(
                item,
                r#" peptide_ref="Pep_{}" rank="1" passThreshold="{}">"#,
                spectrum.peptide_id, spectrum.passes
            );
            self.put(&item)?;
            for &index in evidence.get(&spectrum.peptide_id).into_iter().flatten() {
                self.put(&format!(
                    r#"            <PeptideEvidenceRef peptideEvidence_ref="{}"/>"#,
                    evidence_id(index)
                ))?;
            }
            self.put(&format!(
                r#"            <cvParam cvRef="PSI-MS" accession="MS:1002354" name="PSM-level q-value" value="{}"/>"#,
                spectrum.q_value
            ))?;
            self.put(&format!(
                r#"            <cvParam cvRef="PSI-MS" accession="MS:1002351" name="PSM-level local FDR" value="{}"/>"#,
                spectrum.pep
            ))?;
            self.put(&format!(
                r#"            <userParam name="Tuffold:score" value="{}"/>"#,
                spectrum.score
            ))?;
            self.put("          </SpectrumIdentificationItem>")?;
            if let Some(rt) = spectrum.rt_seconds {
                self.put(&format!(
                    r#"          <cvParam cvRef="PSI-MS" accession="MS:1000016" name="scan start time" value="{rt}" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>"#
                ))?;
            }
            self.put("        </SpectrumIdentificationResult>")?;
        }
        self.put("      </SpectrumIdentificationList>")
    }

    fn groups(&mut self, data: &Exchange) -> io::Result<()> {
        // For each (peptide, protein): its evidence; for each peptide: its
        // spectra.
        let mut evidence: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (index, item) in data.evidence.iter().enumerate() {
            evidence
                .entry((item.peptide_id, item.protein_id))
                .or_default()
                .push(index);
        }
        let mut spectra: HashMap<i64, Vec<i64>> = HashMap::new();
        for spectrum in &data.spectra {
            spectra
                .entry(spectrum.peptide_id)
                .or_default()
                .push(spectrum.id);
        }
        self.put(r#"      <ProteinDetectionList id="PDL">"#)?;
        for group in &data.groups {
            self.put(&format!(
                r#"        <ProteinAmbiguityGroup id="PAG_{}">"#,
                group.id
            ))?;
            for (position, (protein_id, role)) in group.members.iter().enumerate() {
                self.put(&format!(
                    r#"          <ProteinDetectionHypothesis id="PDH_{}_{protein_id}" dBSequence_ref="DBSeq_{protein_id}" passThreshold="{}">"#,
                    group.id,
                    group.passes && position == 0
                ))?;
                for (peptide_id, _) in &group.peptides {
                    let Some(indices) = evidence.get(&(*peptide_id, *protein_id)) else {
                        continue;
                    };
                    for &index in indices {
                        self.put(&format!(
                            r#"            <PeptideHypothesis peptideEvidence_ref="{}">"#,
                            evidence_id(index)
                        ))?;
                        for spectrum in spectra.get(peptide_id).into_iter().flatten() {
                            self.put(&format!(
                                r#"              <SpectrumIdentificationItemRef spectrumIdentificationItem_ref="SII_{spectrum}"/>"#
                            ))?;
                        }
                        self.put("            </PeptideHypothesis>")?;
                    }
                }
                if role == "leader" {
                    self.put(r#"            <cvParam cvRef="PSI-MS" accession="MS:1002403" name="group representative"/>"#)?;
                    self.put(r#"            <cvParam cvRef="PSI-MS" accession="MS:1002401" name="leading protein"/>"#)?;
                } else {
                    self.put(&format!(
                        r#"            <userParam name="Tuffold:member role" value="{role}"/>"#
                    ))?;
                }
                self.put("          </ProteinDetectionHypothesis>")?;
            }
            if let Some(q) = group.q_value {
                self.put(&format!(
                    r#"          <cvParam cvRef="PSI-MS" accession="MS:1002373" name="protein group-level q-value" value="{q}"/>"#
                ))?;
            }
            if let Some(pep) = group.pep {
                self.put(&format!(
                    r#"          <cvParam cvRef="PSI-MS" accession="MS:1002370" name="protein group-level local FDR" value="{pep}"/>"#
                ))?;
            }
            self.put(&format!(
                r#"          <cvParam cvRef="PSI-MS" accession="MS:1002415" name="protein group passes threshold" value="{}"/>"#,
                group.passes
            ))?;
            self.put("        </ProteinAmbiguityGroup>")?;
        }
        self.put("      </ProteinDetectionList>")
    }
}

fn evidence_id(index: usize) -> String {
    format!("PE_{}", index + 1)
}

/// A spectrum's `spectrumID`: the spectra file's own id when every
/// spectrum of its run has one (so `SpectraData` names one format and the
/// ids resolve in the file, as PRIDE's validator looks them up), else
/// Tuffold's key.
fn spectrum_id<'a>(native_runs: &HashSet<i64>, spectrum: &'a ExchangeSpectrum) -> &'a str {
    match &spectrum.native_id {
        Some(native) if native_runs.contains(&spectrum.run_id) => native,
        _ => &spectrum.key,
    }
}

fn native_runs(spectra: &[ExchangeSpectrum]) -> HashSet<i64> {
    let mut runs: HashSet<i64> = spectra.iter().map(|s| s.run_id).collect();
    for spectrum in spectra.iter().filter(|s| s.native_id.is_none()) {
        runs.remove(&spectrum.run_id);
    }
    runs
}

/// The nativeID format a run's spectrum keys follow.
fn id_format(key: &str) -> (&'static str, &'static str) {
    if key.starts_with("scan=") {
        ("MS:1000776", "scan number only nativeID format")
    } else if key.starts_with("index=") {
        ("MS:1000774", "multiple peak list nativeID format")
    } else if key.starts_with("controllerType=") {
        ("MS:1000768", "Thermo nativeID format")
    } else {
        ("MS:1000824", "no nativeID format")
    }
}

/// A `Modification` element: Unimod when the accession is one, else
/// "unknown modification".
fn modification_xml(modification: &Modification, length: usize) -> String {
    let location = match modification.site {
        ModSite::NTerm => 0,
        ModSite::CTerm => length + 1,
        ModSite::Residue(position) => usize::try_from(position).unwrap_or(0),
    };
    let term = match (&modification.accession, &modification.name) {
        (Some(accession), Some(name)) if accession.starts_with("UNIMOD:") => format!(
            r#"<cvParam cvRef="UNIMOD" accession="{}" name="{}"/>"#,
            escape(accession),
            escape(name)
        ),
        (_, name) => format!(
            r#"<cvParam cvRef="PSI-MS" accession="MS:1001460" name="unknown modification"{}/>"#,
            name.as_ref()
                .map_or_else(String::new, |n| format!(r#" value="{}""#, escape(n)))
        ),
    };
    format!(
        r#"      <Modification location="{location}" monoisotopicMassDelta="{}">{term}</Modification>"#,
        modification.mass_delta
    )
}

/// XML attribute and text escaping.
/// A path as an `xs:anyURI`: bytes outside URI characters (spaces,
/// non-ASCII) percent-encoded.
fn uri(path: &str) -> String {
    let mut encoded = String::with_capacity(path.len());
    for byte in path.bytes() {
        if byte.is_ascii_alphanumeric() || b"-._~/:@!$&'()*+,;=".contains(&byte) {
            encoded.push(byte as char);
        } else {
            let _ = write!(encoded, "%{byte:02X}");
        }
    }
    encoded
}

/// The PSI-MS term for a spectra file's format, by extension; the generic
/// parent term when the file isn't known or recognized.
fn file_format(file: Option<&str>) -> (&'static str, &'static str) {
    let extension = file
        .and_then(|file| std::path::Path::new(file).extension())
        .map(|extension| extension.to_string_lossy().to_ascii_lowercase());
    match extension.as_deref() {
        Some("mzml") => ("MS:1000584", "mzML format"),
        Some("mgf") => ("MS:1001062", "Mascot MGF format"),
        Some("raw") => ("MS:1000563", "Thermo RAW format"),
        _ => ("MS:1000560", "mass spectrometer file format"),
    }
}

fn escape(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&apos;"),
            c if (c as u32) < 0x20 && !matches!(c, '\t' | '\n' | '\r') => {}
            c => escaped.push(c),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modifications_use_mzidentml_locations() {
        let named = Modification {
            site: ModSite::Residue(3),
            mass_delta: 15.994915,
            name: Some("Oxidation".into()),
            accession: Some("UNIMOD:35".into()),
        };
        assert_eq!(
            modification_xml(&named, 8),
            r#"      <Modification location="3" monoisotopicMassDelta="15.994915"><cvParam cvRef="UNIMOD" accession="UNIMOD:35" name="Oxidation"/></Modification>"#
        );
        let c_term = Modification {
            site: ModSite::CTerm,
            mass_delta: -0.984016,
            name: None,
            accession: None,
        };
        assert!(modification_xml(&c_term, 8).contains(r#"location="9""#));
        assert!(modification_xml(&c_term, 8).contains("MS:1001460"));
        assert_eq!(escape(r#"a<b & "c""#), "a&lt;b &amp; &quot;c&quot;");
        assert_eq!(id_format("scan=7").0, "MS:1000776");
        assert_eq!(id_format("index=7").0, "MS:1000774");
    }

    fn spectrum(run_id: i64, key: &str, native_id: Option<&str>) -> ExchangeSpectrum {
        ExchangeSpectrum {
            id: 1,
            run_id,
            key: key.into(),
            native_id: native_id.map(Into::into),
            rt_seconds: None,
            peptide_id: 1,
            charge: Some(2),
            experimental_mz: None,
            calculated_mz: None,
            is_decoy: false,
            score: 1.0,
            q_value: 0.0,
            pep: 0.0,
            passes: true,
        }
    }

    #[test]
    fn spectrum_ids_are_the_files_own_where_every_spectrum_of_the_run_has_one() {
        // Run 1: both spectra carry their mzML ids, so those are written
        // (Thermo format). Run 2: one lacks it, so the whole run keeps
        // Tuffold's keys and one id format.
        let thermo = "controllerType=0 controllerNumber=1 scan=7";
        let spectra = [
            spectrum(1, "scan=7", Some(thermo)),
            spectrum(
                1,
                "scan=8",
                Some("controllerType=0 controllerNumber=1 scan=8"),
            ),
            spectrum(2, "scan=7", Some(thermo)),
            spectrum(2, "scan=9", None),
        ];
        let native = native_runs(&spectra);
        assert_eq!(native, HashSet::from([1]));
        assert_eq!(spectrum_id(&native, &spectra[0]), thermo);
        assert_eq!(id_format(spectrum_id(&native, &spectra[0])).0, "MS:1000768");
        assert_eq!(spectrum_id(&native, &spectra[2]), "scan=7");
        assert_eq!(spectrum_id(&native, &spectra[3]), "scan=9");
    }
}
