//! mzIdentML 1.1 / 1.2 → [`ImportBatch`].
//!
//! The schema orders `SequenceCollection` and `AnalysisCollection` before
//! `DataCollection`, so one streaming pass resolves every reference.
//!
//! Mapping conventions:
//! - Each `SpectraData` used by a result is one MS run named after its
//!   `location`'s file name without a spectra extension; each run is its own
//!   biosample until the user groups them.
//! - Each `SpectrumIdentificationResult` is one spectrum keyed from its
//!   `spectrumID` (`scan=N`, `index=N`, else verbatim); results of different
//!   lists on the same `spectrumID` share the spectrum. A `scan number(s)`
//!   cvParam (MS:1001115, which MS-GF+ writes beside `spectrumID="index=N"`)
//!   takes precedence, so these spectra meet the scan-keyed spectra of
//!   pepXML, FragPipe, and mzML. Retention time comes from a `scan start
//!   time` cvParam with its unit; repeated results fill unknown values.
//! - Each `SpectrumIdentificationItem` is one PSM of its `rank`, with
//!   `passThreshold`, both m/z values, and every numeric cvParam/userParam as
//!   a score (cvParams keep their PSI-MS accession).
//! - Modifications: `location` 0 is the N-terminus, length+1 the C-terminus.
//!   The delta is `monoisotopicMassDelta`, else the Unimod table for the
//!   modification's accession; otherwise (including `avgMassDelta` alone,
//!   which is not a monoisotopic delta) the file is rejected with the peptide
//!   named.
//! - Evidence is one entry per protein locus (accession and start).
//! - Decoys follow `PeptideEvidence/@isDecoy` and are marked explicit;
//!   proteins never flagged either way fall back to the accession convention.
//! - The batch's software lists every search engine used by a
//!   `SpectrumIdentification`, joined with " + ".
//! - `ProteinDetectionList` is not imported: protein inference is Tuffold's.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use quick_xml::XmlVersion;
use quick_xml::events::Event;
use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, Spectrum, SpectrumKey, unimod,
};

use crate::xml::{Attrs, XmlImportError, invalid, malformed, reader, run_name};

const SCAN_START_TIME: &str = "MS:1000016";
const SCAN_NUMBERS: &str = "MS:1001115";

#[derive(Default)]
struct DraftSoftware {
    attribute_name: Option<String>,
    cv_name: Option<String>,
    version: Option<String>,
}

#[derive(Default)]
struct DraftModification {
    location: u32,
    mono: Option<f64>,
    average: Option<f64>,
    name: Option<String>,
    accession: Option<String>,
}

#[derive(Default)]
struct DraftPeptide {
    sequence: String,
    modifications: Vec<DraftModification>,
}

struct DraftEvidence {
    peptide: String,
    protein: String,
    start: Option<u32>,
    end: Option<u32>,
    pre: Option<char>,
    post: Option<char>,
}

struct DraftResult {
    run: String,
    key: SpectrumKey,
    source_id: String,
    rt_seconds: Option<f64>,
    first_mz: Option<f64>,
    first_charge: Option<i16>,
    /// From a `scan number(s)` cvParam; overrides the `spectrumID` key.
    scan: Option<u32>,
    items: Vec<Psm>,
}

#[derive(Default)]
struct State {
    stack: Vec<String>,
    text: String,
    software: HashMap<String, DraftSoftware>,
    current_software: Option<String>,
    proteins: HashMap<String, Protein>,
    current_protein: Option<String>,
    peptides: HashMap<String, DraftPeptide>,
    current_peptide: Option<String>,
    current_modification: Option<DraftModification>,
    evidence: HashMap<String, DraftEvidence>,
    decoy_flags: HashMap<String, bool>,
    protocol_software: HashMap<String, String>,
    list_protocol: HashMap<String, String>,
    spectra_data: HashMap<String, String>,
    result: Option<DraftResult>,
    item: Option<Psm>,
    item_peptide_ref: String,
    engines: Vec<String>,
    batch: ImportBatch,
    /// Position of each (run, key) spectrum in `batch.spectra`.
    spectra_seen: HashMap<(String, SpectrumKey), usize>,
    runs_seen: HashSet<String>,
    protein_order: Vec<String>,
    proteins_used: HashSet<String>,
}

pub fn import_mzid(path: &Path) -> Result<ImportBatch, XmlImportError> {
    let mut reader = reader(path)?;
    let mut state = State::default();
    let mut buffer = Vec::new();
    loop {
        let offset = reader.buffer_position();
        let event = reader
            .read_event_into(&mut buffer)
            .map_err(|error| malformed(path, reader.error_position(), error))?;
        let is_empty = matches!(event, Event::Empty(_));
        match event {
            Event::Start(start) | Event::Empty(start) => {
                let attrs =
                    Attrs::read(&start).map_err(|problem| malformed(path, offset, problem))?;
                state.stack.push(attrs.element.clone());
                state.text.clear();
                state
                    .start(&attrs)
                    .map_err(|problem| invalid(path, offset, &attrs.element, problem))?;
                if is_empty {
                    let element = state.stack.pop().unwrap_or_default();
                    state
                        .end(&element)
                        .map_err(|problem| invalid(path, offset, &element, problem))?;
                }
            }
            Event::Text(text) => state
                .text
                .push_str(&text.xml_content(XmlVersion::Implicit1_0)),
            Event::End(_) => {
                let element = state.stack.pop().unwrap_or_default();
                state
                    .end(&element)
                    .map_err(|problem| invalid(path, offset, &element, problem))?;
            }
            Event::Eof => break,
            _ => {}
        }
        buffer.clear();
    }
    Ok(state.finish())
}

impl State {
    fn parent(&self) -> &str {
        // The stack holds the current element last; its parent is before it.
        self.stack
            .len()
            .checked_sub(2)
            .map_or("", |index| self.stack[index].as_str())
    }

    fn start(&mut self, attrs: &Attrs) -> Result<(), String> {
        match attrs.element.as_str() {
            "AnalysisSoftware" => {
                let id = attrs.required("id")?.to_string();
                self.software.insert(
                    id.clone(),
                    DraftSoftware {
                        attribute_name: attrs.get("name").map(str::to_string),
                        cv_name: None,
                        version: attrs.get("version").map(str::to_string),
                    },
                );
                self.current_software = Some(id);
            }
            "DBSequence" => {
                let id = attrs.required("id")?.to_string();
                let accession = attrs.required("accession")?.to_string();
                self.proteins.insert(
                    id.clone(),
                    Protein {
                        accession,
                        description: None,
                        sequence: None,
                        is_decoy: false,
                        decoy_explicit: false,
                    },
                );
                self.current_protein = Some(id);
            }
            "Peptide" => {
                let id = attrs.required("id")?.to_string();
                self.peptides.insert(id.clone(), DraftPeptide::default());
                self.current_peptide = Some(id);
            }
            "Modification" if self.current_peptide.is_some() => {
                self.current_modification = Some(DraftModification {
                    location: attrs
                        .parse("location")?
                        .ok_or("is missing its required 'location' attribute")?,
                    mono: attrs.number("monoisotopicMassDelta")?,
                    average: attrs.number("avgMassDelta")?,
                    ..DraftModification::default()
                });
            }
            "SubstitutionModification" => {
                return Err("is not supported yet (amino-acid substitutions)".into());
            }
            "PeptideEvidence" => {
                let id = attrs.required("id")?.to_string();
                let protein = attrs.required("dBSequence_ref")?.to_string();
                if let Some(flag) = attrs.get("isDecoy") {
                    let is_decoy = flag.eq_ignore_ascii_case("true") || flag == "1";
                    *self.decoy_flags.entry(protein.clone()).or_default() |= is_decoy;
                }
                let residue = |name: &str| {
                    attrs
                        .get(name)
                        .and_then(|text| text.chars().next())
                        .filter(|ch| *ch != '-' && *ch != '?')
                };
                self.evidence.insert(
                    id,
                    DraftEvidence {
                        peptide: attrs.required("peptide_ref")?.to_string(),
                        protein,
                        start: attrs.parse("start")?,
                        end: attrs.parse("end")?,
                        pre: residue("pre"),
                        post: residue("post"),
                    },
                );
            }
            "SpectrumIdentification" => {
                let list = attrs.required("spectrumIdentificationList_ref")?;
                let protocol = attrs.required("spectrumIdentificationProtocol_ref")?;
                self.list_protocol
                    .insert(list.to_string(), protocol.to_string());
            }
            "SpectrumIdentificationProtocol" => {
                let id = attrs.required("id")?;
                let software = attrs.required("analysisSoftware_ref")?;
                self.protocol_software
                    .insert(id.to_string(), software.to_string());
            }
            "SpectraData" => {
                let id = attrs.required("id")?.to_string();
                let location = attrs.get("location").filter(|value| !value.is_empty());
                let name = match location {
                    Some(location) => run_name(&percent_decoded(location)),
                    None => run_name(attrs.get("name").unwrap_or(&id)),
                };
                self.spectra_data.insert(id, name);
            }
            "SpectrumIdentificationList" => {
                let list = attrs.required("id")?;
                let engine = self.engine_for_list(list);
                if !self.engines.contains(&engine) {
                    self.engines.push(engine);
                }
            }
            "SpectrumIdentificationResult" => self.start_result(attrs)?,
            "SpectrumIdentificationItem" => self.start_item(attrs)?,
            "PeptideEvidenceRef" => self.add_evidence(attrs.required("peptideEvidence_ref")?)?,
            "cvParam" => self.cv_param(attrs)?,
            "userParam" => self.user_param(attrs)?,
            _ => {}
        }
        Ok(())
    }

    fn end(&mut self, element: &str) -> Result<(), String> {
        match element {
            "AnalysisSoftware" => self.current_software = None,
            "DBSequence" => self.current_protein = None,
            "Seq" => {
                if let Some(protein) = self
                    .current_protein
                    .as_ref()
                    .and_then(|id| self.proteins.get_mut(id))
                {
                    let sequence: String = self
                        .text
                        .chars()
                        .filter(|ch| !ch.is_whitespace())
                        .collect::<String>()
                        .to_ascii_uppercase();
                    protein.sequence = (!sequence.is_empty()).then_some(sequence);
                }
            }
            "PeptideSequence" => {
                if let Some(peptide) = self
                    .current_peptide
                    .as_ref()
                    .and_then(|id| self.peptides.get_mut(id))
                {
                    peptide.sequence = self.text.trim().to_ascii_uppercase();
                }
            }
            "Modification" => {
                if let (Some(modification), Some(peptide)) = (
                    self.current_modification.take(),
                    self.current_peptide
                        .as_ref()
                        .and_then(|id| self.peptides.get_mut(id)),
                ) {
                    peptide.modifications.push(modification);
                }
            }
            "Peptide" => self.current_peptide = None,
            "SpectrumIdentificationItem" => {
                if let (Some(item), Some(result)) = (self.item.take(), self.result.as_mut()) {
                    result.first_mz = result.first_mz.or(item.experimental_mz);
                    result.first_charge = result.first_charge.or(Some(item.charge));
                    result.items.push(item);
                }
            }
            "SpectrumIdentificationResult" => self.finish_result(),
            _ => {}
        }
        Ok(())
    }

    fn engine_for_list(&self, list: &str) -> String {
        let software = self
            .list_protocol
            .get(list)
            .and_then(|protocol| self.protocol_software.get(protocol))
            .and_then(|id| self.software.get(id));
        software
            .and_then(|software| {
                software
                    .cv_name
                    .clone()
                    .or_else(|| software.attribute_name.clone())
            })
            .unwrap_or_else(|| "unknown search engine (mzIdentML)".into())
    }

    fn cv_param(&mut self, attrs: &Attrs) -> Result<(), String> {
        let name = attrs.get("name").unwrap_or_default().to_string();
        let accession = attrs.get("accession").map(str::to_string);
        match self.parent() {
            "SoftwareName" => {
                if let Some(software) = self
                    .current_software
                    .as_ref()
                    .and_then(|id| self.software.get_mut(id))
                {
                    software.cv_name = Some(name);
                }
            }
            "DBSequence" if accession.as_deref() == Some("MS:1001088") => {
                if let Some(protein) = self
                    .current_protein
                    .as_ref()
                    .and_then(|id| self.proteins.get_mut(id))
                {
                    let description = attrs.get("value").unwrap_or_default();
                    let description = description.strip_prefix('>').unwrap_or(description);
                    protein.description = Some(description.to_string());
                }
            }
            "Modification" => {
                if let Some(modification) = self.current_modification.as_mut()
                    && modification.accession.is_none()
                {
                    modification.name = Some(name);
                    modification.accession = accession;
                }
            }
            "SpectrumIdentificationResult" if accession.as_deref() == Some(SCAN_NUMBERS) => {
                // "1234" or "1234,1235" (merged scans): the first scan keys it.
                let value = attrs.required("value")?;
                let first = value.split([',', ' ']).next().unwrap_or_default().trim();
                let scan = first
                    .parse()
                    .map_err(|_| format!("has scan number(s) {value:?}, which is not a scan"))?;
                if let Some(result) = self.result.as_mut() {
                    result.scan = Some(scan);
                }
            }
            "SpectrumIdentificationResult" if accession.as_deref() == Some(SCAN_START_TIME) => {
                let value = attrs.required_number("value")?;
                let unit = attrs
                    .get("unitName")
                    .or(attrs.get("unitAccession"))
                    .unwrap_or_default()
                    .to_ascii_lowercase();
                let seconds = match unit.as_str() {
                    "second" | "uo:0000010" => Some(value),
                    "minute" | "uo:0000031" => Some(value * 60.0),
                    _ => None, // no unit stated: leave the time unknown, never guess
                };
                if let Some(result) = self.result.as_mut() {
                    result.rt_seconds = seconds;
                }
            }
            "SpectrumIdentificationItem" => self.add_score(name, accession, attrs)?,
            _ => {}
        }
        Ok(())
    }

    fn user_param(&mut self, attrs: &Attrs) -> Result<(), String> {
        match self.parent() {
            "SoftwareName" => {
                if let Some(software) = self
                    .current_software
                    .as_ref()
                    .and_then(|id| self.software.get_mut(id))
                {
                    software.cv_name = attrs.get("name").map(str::to_string);
                }
            }
            "SpectrumIdentificationItem" => {
                let name = attrs.required("name")?.to_string();
                self.add_score(name, None, attrs)?;
            }
            _ => {}
        }
        Ok(())
    }

    fn add_score(
        &mut self,
        name: String,
        accession: Option<String>,
        attrs: &Attrs,
    ) -> Result<(), String> {
        // Valueless params are flags (e.g. "peptide unique to one protein").
        let Ok(Some(value)) = attrs.number("value") else {
            return Ok(());
        };
        if let Some(item) = self.item.as_mut() {
            item.scores.push(Score {
                name,
                accession,
                value,
            });
        }
        Ok(())
    }

    fn start_result(&mut self, attrs: &Attrs) -> Result<(), String> {
        let spectra_ref = attrs.required("spectraData_ref")?;
        let run = self
            .spectra_data
            .get(spectra_ref)
            .cloned()
            .ok_or_else(|| format!("references unknown SpectraData '{spectra_ref}'"))?;
        let spectrum_id = attrs.required("spectrumID")?;
        self.result = Some(DraftResult {
            run,
            key: SpectrumKey::from_native_id(spectrum_id),
            source_id: spectrum_id.to_string(),
            rt_seconds: None,
            first_mz: None,
            first_charge: None,
            scan: None,
            items: Vec::new(),
        });
        Ok(())
    }

    fn start_item(&mut self, attrs: &Attrs) -> Result<(), String> {
        let result = self
            .result
            .as_ref()
            .ok_or("appears outside any <SpectrumIdentificationResult>")?;
        let peptide_ref = attrs.required("peptide_ref")?;
        let draft = self
            .peptides
            .get(peptide_ref)
            .ok_or_else(|| format!("references unknown Peptide '{peptide_ref}'"))?;
        let peptide = peptide_form(peptide_ref, draft)?;
        self.item_peptide_ref = peptide_ref.to_string();
        let pass_threshold = attrs
            .get("passThreshold")
            .map(|flag| flag.eq_ignore_ascii_case("true") || flag == "1");
        self.item = Some(Psm {
            ms_run: result.run.clone(),
            spectrum: result.key.clone(),
            rank: attrs
                .parse("rank")?
                .ok_or("is missing its required 'rank' attribute")?,
            peptide,
            charge: attrs
                .parse("chargeState")?
                .ok_or("is missing its required 'chargeState' attribute")?,
            calculated_mz: attrs.number("calculatedMassToCharge")?,
            experimental_mz: attrs.number("experimentalMassToCharge")?,
            pass_threshold,
            scores: Vec::new(),
            evidence: Vec::new(),
        });
        Ok(())
    }

    fn add_evidence(&mut self, reference: &str) -> Result<(), String> {
        let draft = self
            .evidence
            .get(reference)
            .ok_or_else(|| format!("references unknown PeptideEvidence '{reference}'"))?;
        if self.item.is_some() && draft.peptide != self.item_peptide_ref {
            return Err(format!(
                "references PeptideEvidence '{reference}', which belongs to Peptide '{}', not '{}'",
                draft.peptide, self.item_peptide_ref
            ));
        }
        let protein = self.proteins.get(&draft.protein).ok_or_else(|| {
            format!(
                "evidence '{reference}' references unknown DBSequence '{}'",
                draft.protein
            )
        })?;
        if let Some(item) = self.item.as_mut() {
            let is_new_locus = !item.evidence.iter().any(|evidence| {
                evidence.accession == protein.accession && evidence.start == draft.start
            });
            if is_new_locus {
                item.evidence.push(Evidence {
                    accession: protein.accession.clone(),
                    start: draft.start,
                    end: draft.end,
                    pre: draft.pre,
                    post: draft.post,
                });
            }
        }
        if self.proteins_used.insert(draft.protein.clone()) {
            self.protein_order.push(draft.protein.clone());
        }
        Ok(())
    }

    fn finish_result(&mut self) {
        let Some(mut result) = self.result.take() else {
            return;
        };
        if let Some(scan) = result.scan {
            result.key = SpectrumKey::Scan(scan);
            for item in &mut result.items {
                item.spectrum = result.key.clone();
            }
        }
        if self.runs_seen.insert(result.run.clone()) {
            self.batch.biosamples.push(BioSample {
                name: result.run.clone(),
                category: None,
            });
            self.batch.ms_runs.push(MsRun {
                name: result.run.clone(),
                biosample: result.run.clone(),
                spectra_file: None,
            });
        }
        let spectrum = Spectrum {
            ms_run: result.run.clone(),
            key: result.key.clone(),
            source_id: Some(result.source_id),
            precursor_mz: result.first_mz,
            precursor_charge: result.first_charge,
            rt_seconds: result.rt_seconds,
        };
        let seen_key = (result.run.clone(), result.key.clone());
        match self.spectra_seen.get(&seen_key) {
            Some(&position) => self.batch.spectra[position].merge(&spectrum),
            None => {
                self.spectra_seen.insert(seen_key, self.batch.spectra.len());
                self.batch.spectra.push(spectrum);
            }
        }
        self.batch.psms.extend(result.items);
    }

    fn finish(mut self) -> ImportBatch {
        let mut batch = std::mem::take(&mut self.batch);
        batch.software = Software {
            version: match self.engines.len() {
                1 => self
                    .software
                    .values()
                    .find(|software| {
                        software
                            .cv_name
                            .as_ref()
                            .or(software.attribute_name.as_ref())
                            == self.engines.first()
                    })
                    .and_then(|software| software.version.clone()),
                _ => None,
            },
            name: if self.engines.is_empty() {
                "unknown search engine (mzIdentML)".into()
            } else {
                self.engines.join(" + ")
            },
        };
        batch.proteins = self
            .protein_order
            .iter()
            .filter_map(|id| {
                let mut protein = self.proteins.remove(id)?;
                match self.decoy_flags.get(id) {
                    Some(flag) => {
                        protein.is_decoy = *flag;
                        protein.decoy_explicit = true;
                    }
                    None => protein.is_decoy = is_decoy_accession(&protein.accession),
                }
                Some(protein)
            })
            .collect();
        batch
    }
}

fn peptide_form(id: &str, draft: &DraftPeptide) -> Result<PeptideForm, String> {
    let length = u32::try_from(draft.sequence.chars().count()).unwrap_or(u32::MAX);
    if length == 0 {
        return Err(format!(
            "references Peptide '{id}', which has no PeptideSequence"
        ));
    }
    let modifications = draft
        .modifications
        .iter()
        .map(|modification| {
            let site = match modification.location {
                0 => ModSite::NTerm,
                location if location == length + 1 => ModSite::CTerm,
                location if location <= length => ModSite::Residue(location),
                location => {
                    return Err(format!(
                        "references Peptide '{id}' with a modification at location {location} \
                         outside its {length} residues"
                    ));
                }
            };
            let mass_delta = modification
                .mono
                .or_else(|| {
                    modification
                        .accession
                        .as_deref()
                        .and_then(unimod::mass_delta)
                })
                .ok_or_else(|| {
                    let average = modification
                        .average
                        .map(|delta| format!(" (only an average delta, {delta}, is given)"))
                        .unwrap_or_default();
                    format!(
                        "references Peptide '{id}' with a modification ({}) that has no \
                         monoisotopic mass delta{average} and is not in Tuffold's Unimod table",
                        modification.accession.as_deref().unwrap_or("no accession")
                    )
                })?;
            Ok(Modification {
                site,
                mass_delta,
                name: modification.name.clone(),
                accession: modification.accession.clone(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    PeptideForm::new(&draft.sequence, modifications)
        .map_err(|error| format!("references Peptide '{id}': {error}"))
}

/// `location` is a URI: `%20` and the like stand for the bytes they
/// encode. A `%` not followed by two hex digits is kept as it is.
fn percent_decoded(location: &str) -> String {
    let bytes = location.as_bytes();
    let mut decoded = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        let hex = bytes
            .get(i + 1..i + 3)
            .and_then(|pair| std::str::from_utf8(pair).ok())
            .and_then(|pair| u8::from_str_radix(pair, 16).ok());
        match (bytes[i], hex) {
            (b'%', Some(byte)) => {
                decoded.push(byte);
                i += 3;
            }
            (byte, _) => {
                decoded.push(byte);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&decoded).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locations_are_percent_decoded() {
        assert_eq!(percent_decoded("/data/run%201.mzML"), "/data/run 1.mzML");
        assert_eq!(percent_decoded("caf%C3%A9.mgf"), "café.mgf");
        assert_eq!(percent_decoded("100%_run"), "100%_run");
        assert_eq!(percent_decoded("x%2"), "x%2");
    }
}
