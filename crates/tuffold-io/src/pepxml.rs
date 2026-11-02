//! pepXML → [`ImportBatch`] (TPP pepXML, any search engine).
//!
//! Mapping conventions:
//! - Each `msms_run_summary` is one MS run named after the file name of its
//!   `base_name` without a spectra extension; its spectra file is that name
//!   plus `raw_data` (e.g. `.mzML`). Each run is its own biosample until the
//!   user groups them.
//! - Each `spectrum_query` is one spectrum keyed by `start_scan` (else the
//!   `spectrum` title verbatim; the query `index` is a 1-based ordinal, not a
//!   spectrum index, so it never keys). Its precursor m/z follows from
//!   `precursor_neutral_mass` and `assumed_charge`. Queries for the same scan
//!   (e.g. charge 2 and 3) share one spectrum whose charge becomes unknown.
//! - Each `search_hit` is one PSM of its `hit_rank`. Every `search_score` is
//!   kept under its verbatim name. `analysis_result`s (PeptideProphet,
//!   iProphet, Percolator) are not imported: they are post-processors' verdicts
//!   trained on the decoys, and Tuffold validates the engine, not them.
//! - Modifications: a `mod_aminoacid_mass` is the *total* residue mass. The
//!   delta is, in order of preference: its `variable`/`static` attributes; the
//!   `massdiff` of the matching `aminoacid_modification` (within 0.01 Da) of
//!   the `search_summary` the hit's `search_result/@search_id` names; otherwise
//!   the mass minus the residue mass of the summary's declared mass type
//!   (`fragment_mass_type`, else `precursor_mass_type`, else monoisotopic).
//!   Terminal masses work the same against the unmodified H / OH groups; null
//!   deltas are dropped, and positions outside the peptide are errors.
//! - Evidence: `protein` plus every `alternative_protein`, with flanking
//!   residues. A protein is a decoy by accession convention, or when its
//!   description ends in `:reversed` (the Inspect convention).
//! - The batch's software is the first `search_summary`'s engine. MSFragger
//!   calls itself `X! Tandem` there (so TPP reads its files) and names
//!   itself in the version (`MSFragger-4.1`); it is imported as MSFragger.

use std::collections::HashMap;
use std::path::Path;

use quick_xml::events::Event;
use tuffold_core::chemistry::{
    C_TERM_AVERAGE, C_TERM_MONO, MassType, N_TERM_AVERAGE, N_TERM_MONO, mz, residue_mass,
};
use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, Spectrum, SpectrumKey,
};

use crate::xml::{Attrs, XmlImportError, invalid, malformed, reader, run_name};

/// Masses this close to a declared modification take its declared delta.
const MASS_MATCH_TOLERANCE: f64 = 0.01;

#[derive(Debug, Clone, Copy)]
struct DeclaredMod {
    residue: char,
    mass: f64,
    massdiff: f64,
}

#[derive(Debug, Clone, Copy)]
struct DeclaredTermMod {
    is_n_term: bool,
    mass: f64,
    massdiff: f64,
}

/// The `spectrum_query` being read.
struct Query {
    run: String,
    key: SpectrumKey,
    charge: i16,
    neutral_mass: Option<f64>,
}

/// The `search_hit` being read.
struct Hit {
    rank: u32,
    sequence: String,
    calculated_neutral: Option<f64>,
    modifications: Vec<Modification>,
    scores: Vec<Score>,
    evidence: Vec<Evidence>,
}

#[derive(Default)]
struct State {
    batch: ImportBatch,
    run: Option<String>,
    /// Search summaries of the current run, by `search_id`.
    summaries: Vec<Summary>,
    /// The `search_id` of the `search_result` being read.
    search_id: Option<u32>,
    query: Option<Query>,
    hit: Option<Hit>,
    analysis: Option<String>,
    proteins: HashMap<String, Protein>,
    protein_order: Vec<String>,
    /// Position of each (run, key) spectrum in `batch.spectra`.
    spectra_seen: HashMap<(String, SpectrumKey), usize>,
}

/// One `search_summary`: its declared modifications and mass type.
struct Summary {
    search_id: u32,
    mass_type: MassType,
    residue_mods: Vec<DeclaredMod>,
    terminal_mods: Vec<DeclaredTermMod>,
}

/// The engine a `search_summary` names, seeing through MSFragger's
/// `search_engine="X! Tandem"`.
fn software(engine: &str, version: Option<&str>) -> Software {
    const MSFRAGGER: &str = "MSFragger";
    match version.and_then(|version| version.strip_prefix(MSFRAGGER)) {
        Some(rest) => Software {
            name: MSFRAGGER.into(),
            version: Some(rest.trim_start_matches(['-', ' ', '_']).to_string())
                .filter(|version| !version.is_empty()),
        },
        None => Software {
            name: engine.to_string(),
            version: version.map(str::to_string),
        },
    }
}

pub fn import_pepxml(path: &Path) -> Result<ImportBatch, XmlImportError> {
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
                state
                    .start(&attrs)
                    .map_err(|problem| invalid(path, offset, &attrs.element, problem))?;
                if is_empty {
                    state
                        .end(&attrs.element)
                        .map_err(|problem| invalid(path, offset, &attrs.element, problem))?;
                }
            }
            Event::End(end) => {
                let element = end.local_name().as_ref().to_string();
                state
                    .end(&element)
                    .map_err(|problem| invalid(path, offset, &element, problem))?;
            }
            Event::Eof => break,
            _ => {}
        }
        buffer.clear();
    }
    let mut batch = state.batch;
    batch.proteins = state
        .protein_order
        .iter()
        .filter_map(|accession| state.proteins.remove(accession))
        .collect();
    if batch.software.name.is_empty() {
        batch.software.name = "unknown search engine (pepXML)".into();
    }
    Ok(batch)
}

impl State {
    fn start(&mut self, attrs: &Attrs) -> Result<(), String> {
        match attrs.element.as_str() {
            "msms_run_summary" => self.start_run(attrs),
            "search_summary" => {
                if self.batch.software.name.is_empty() {
                    self.batch.software = software(
                        attrs.get("search_engine").unwrap_or_default(),
                        attrs.get("search_engine_version"),
                    );
                }
                let mass_type = attrs
                    .get("fragment_mass_type")
                    .or(attrs.get("precursor_mass_type"))
                    .map_or(MassType::Monoisotopic, |kind| {
                        if kind.eq_ignore_ascii_case("average") {
                            MassType::Average
                        } else {
                            MassType::Monoisotopic
                        }
                    });
                let search_id = match attrs.parse::<u32>("search_id")? {
                    Some(id) => id,
                    None => u32::try_from(self.summaries.len() + 1).unwrap_or(u32::MAX),
                };
                self.summaries.push(Summary {
                    search_id,
                    mass_type,
                    residue_mods: Vec::new(),
                    terminal_mods: Vec::new(),
                });
                Ok(())
            }
            "search_result" => {
                self.search_id = attrs.parse("search_id")?;
                Ok(())
            }
            "aminoacid_modification" => {
                let residue = attrs.required("aminoacid")?.chars().next().unwrap_or('X');
                let summary = self
                    .summaries
                    .last_mut()
                    .ok_or("appears outside any <search_summary>")?;
                summary.residue_mods.push(DeclaredMod {
                    residue: residue.to_ascii_uppercase(),
                    mass: attrs.required_number("mass")?,
                    massdiff: attrs.required_number("massdiff")?,
                });
                Ok(())
            }
            "terminal_modification" => {
                let summary = self
                    .summaries
                    .last_mut()
                    .ok_or("appears outside any <search_summary>")?;
                summary.terminal_mods.push(DeclaredTermMod {
                    is_n_term: attrs.required("terminus")?.eq_ignore_ascii_case("n"),
                    mass: attrs.required_number("mass")?,
                    massdiff: attrs.required_number("massdiff")?,
                });
                Ok(())
            }
            "spectrum_query" => self.start_query(attrs),
            "search_hit" => self.start_hit(attrs),
            "alternative_protein" => {
                let accession = attrs.required("protein")?;
                self.add_evidence(accession, attrs.get("protein_descr"), attrs);
                Ok(())
            }
            "modification_info" => self.terminal_modifications(attrs),
            "mod_aminoacid_mass" => self.residue_modification(attrs),
            // Scores inside an analysis_result belong to the post-processor.
            "search_score" if self.analysis.is_none() => {
                let name = attrs.required("name")?;
                // Non-numeric scores (rare, e.g. text flags) carry no value to rank by.
                if let Ok(Some(value)) = attrs.number("value") {
                    self.hit_mut()?.scores.push(Score {
                        name: name.to_string(),
                        accession: None,
                        value,
                    });
                }
                Ok(())
            }
            "analysis_result" => {
                self.analysis = attrs.get("analysis").map(str::to_string);
                Ok(())
            }
            _ => Ok(()),
        }
    }

    fn end(&mut self, element: &str) -> Result<(), String> {
        match element {
            "search_hit" => self.finish_hit()?,
            "search_result" => self.search_id = None,
            "spectrum_query" => self.query = None,
            "analysis_result" => self.analysis = None,
            "msms_run_summary" => {
                self.run = None;
                self.summaries.clear();
            }
            _ => {}
        }
        Ok(())
    }

    /// The summary the current hit belongs to: by `search_id`, else the first.
    fn summary(&self) -> Option<&Summary> {
        self.search_id
            .and_then(|id| {
                self.summaries
                    .iter()
                    .find(|summary| summary.search_id == id)
            })
            .or_else(|| self.summaries.first())
    }

    fn start_run(&mut self, attrs: &Attrs) -> Result<(), String> {
        let name = run_name(attrs.required("base_name")?);
        let spectra_file = attrs
            .get("raw_data")
            .filter(|extension| !extension.is_empty())
            .map(|extension| format!("{name}{extension}"));
        if !self.batch.ms_runs.iter().any(|run| run.name == name) {
            self.batch.biosamples.push(BioSample {
                name: name.clone(),
                category: None,
            });
            self.batch.ms_runs.push(MsRun {
                name: name.clone(),
                biosample: name.clone(),
                spectra_file,
            });
        }
        self.run = Some(name);
        Ok(())
    }

    fn start_query(&mut self, attrs: &Attrs) -> Result<(), String> {
        let run = self
            .run
            .clone()
            .ok_or("appears outside any <msms_run_summary>")?;
        let title = attrs.get("spectrum").filter(|title| !title.is_empty());
        let key = match (attrs.parse::<u32>("start_scan")?, title) {
            (Some(scan), _) => SpectrumKey::Scan(scan),
            (None, Some(title)) => SpectrumKey::Native(title.to_string()),
            (None, None) => return Err("has neither start_scan nor spectrum".into()),
        };
        let charge: i16 = attrs
            .parse("assumed_charge")?
            .ok_or("is missing its required 'assumed_charge' attribute")?;
        let neutral_mass = attrs.number("precursor_neutral_mass")?;
        let charge_for_mz = u8::try_from(charge).ok().filter(|charge| *charge > 0);
        let spectrum = Spectrum {
            ms_run: run.clone(),
            key: key.clone(),
            source_id: title.map(str::to_string),
            precursor_mz: neutral_mass.zip(charge_for_mz).map(|(mass, z)| mz(mass, z)),
            precursor_charge: Some(charge),
            rt_seconds: attrs.number("retention_time_sec")?,
        };
        match self.spectra_seen.get(&(run.clone(), key.clone())) {
            Some(&position) => self.batch.spectra[position].merge(&spectrum),
            None => {
                self.spectra_seen
                    .insert((run.clone(), key.clone()), self.batch.spectra.len());
                self.batch.spectra.push(spectrum);
            }
        }
        self.query = Some(Query {
            run,
            key,
            charge,
            neutral_mass,
        });
        Ok(())
    }

    fn start_hit(&mut self, attrs: &Attrs) -> Result<(), String> {
        if self.query.is_none() {
            return Err("appears outside any <spectrum_query>".into());
        }
        let sequence = attrs.required("peptide")?.to_ascii_uppercase();
        self.hit = Some(Hit {
            rank: attrs
                .parse("hit_rank")?
                .ok_or("is missing its required 'hit_rank' attribute")?,
            sequence,
            calculated_neutral: attrs.number("calc_neutral_pep_mass")?,
            modifications: Vec::new(),
            scores: Vec::new(),
            evidence: Vec::new(),
        });
        let accession = attrs.required("protein")?;
        self.add_evidence(accession, attrs.get("protein_descr"), attrs);
        Ok(())
    }

    fn hit_mut(&mut self) -> Result<&mut Hit, String> {
        self.hit
            .as_mut()
            .ok_or_else(|| "appears outside any <search_hit>".to_string())
    }

    fn add_evidence(&mut self, accession: &str, description: Option<&str>, attrs: &Attrs) {
        let description = description.filter(|text| !text.is_empty());
        let is_decoy = is_decoy_accession(accession)
            || description.is_some_and(|text| text.ends_with(":reversed"));
        let residue = |name: &str| attrs.get(name).and_then(|text| text.chars().next());
        if let Some(hit) = self.hit.as_mut()
            && !hit
                .evidence
                .iter()
                .any(|evidence| evidence.accession == accession)
        {
            hit.evidence.push(Evidence {
                accession: accession.to_string(),
                start: None,
                end: None,
                pre: residue("peptide_prev_aa").filter(|ch| *ch != '-'),
                post: residue("peptide_next_aa").filter(|ch| *ch != '-'),
            });
        }
        let protein = self
            .proteins
            .entry(accession.to_string())
            .or_insert_with(|| {
                self.protein_order.push(accession.to_string());
                Protein {
                    accession: accession.to_string(),
                    description: None,
                    sequence: None,
                    is_decoy: false,
                    decoy_explicit: false,
                }
            });
        protein.is_decoy |= is_decoy;
        if protein.description.is_none() {
            protein.description = description.map(str::to_string);
        }
    }

    fn terminal_modifications(&mut self, attrs: &Attrs) -> Result<(), String> {
        let (terminal_mods, mass_type) = self
            .summary()
            .map_or((Vec::new(), MassType::Monoisotopic), |summary| {
                (summary.terminal_mods.clone(), summary.mass_type)
            });
        let (n_term, c_term) = match mass_type {
            MassType::Monoisotopic => (N_TERM_MONO, C_TERM_MONO),
            MassType::Average => (N_TERM_AVERAGE, C_TERM_AVERAGE),
        };
        for (attribute, is_n_term, unmodified) in [
            ("mod_nterm_mass", true, n_term),
            ("mod_cterm_mass", false, c_term),
        ] {
            let Some(mass) = attrs.number(attribute)? else {
                continue;
            };
            let declared = terminal_mods.iter().find(|declared| {
                declared.is_n_term == is_n_term
                    && (declared.mass - mass).abs() < MASS_MATCH_TOLERANCE
            });
            let mass_delta = declared.map_or(mass - unmodified, |declared| declared.massdiff);
            self.hit_mut()?.modifications.push(Modification {
                site: if is_n_term {
                    ModSite::NTerm
                } else {
                    ModSite::CTerm
                },
                mass_delta,
                name: None,
                accession: None,
            });
        }
        Ok(())
    }

    fn residue_modification(&mut self, attrs: &Attrs) -> Result<(), String> {
        let position: u32 = attrs
            .parse("position")?
            .ok_or("is missing its required 'position' attribute")?;
        let mass = attrs.required_number("mass")?;
        let explicit = match (attrs.number("variable")?, attrs.number("static")?) {
            (None, None) => None,
            (variable, fixed) => Some(variable.unwrap_or(0.0) + fixed.unwrap_or(0.0)),
        };
        let (residue_mods, mass_type) = self
            .summary()
            .map_or((Vec::new(), MassType::Monoisotopic), |summary| {
                (summary.residue_mods.clone(), summary.mass_type)
            });
        let hit = self.hit_mut()?;
        let index = usize::try_from(position).unwrap_or(usize::MAX);
        let residue = position
            .checked_sub(1)
            .and_then(|_| hit.sequence.chars().nth(index - 1))
            .ok_or_else(|| format!("has position {position} outside peptide {}", hit.sequence))?;
        let mass_delta = match explicit {
            Some(delta) => delta,
            None => residue_mods
                .iter()
                .find(|declared| {
                    declared.residue == residue
                        && (declared.mass - mass).abs() < MASS_MATCH_TOLERANCE
                })
                .map(|declared| declared.massdiff)
                .or_else(|| residue_mass(residue, mass_type).map(|base| mass - base))
                .ok_or_else(|| {
                    format!("modifies residue '{residue}', which has no defined mass")
                })?,
        };
        hit.modifications.push(Modification {
            site: ModSite::Residue(position),
            mass_delta,
            name: None,
            accession: None,
        });
        Ok(())
    }

    fn finish_hit(&mut self) -> Result<(), String> {
        let (Some(hit), Some(query)) = (self.hit.take(), self.query.as_ref()) else {
            return Ok(());
        };
        let peptide = PeptideForm::new(&hit.sequence, hit.modifications)
            .map_err(|error| error.to_string())?;
        let charge_for_mz = u8::try_from(query.charge).ok().filter(|charge| *charge > 0);
        self.batch.psms.push(Psm {
            ms_run: query.run.clone(),
            spectrum: query.key.clone(),
            rank: hit.rank,
            peptide,
            charge: query.charge,
            calculated_mz: hit
                .calculated_neutral
                .zip(charge_for_mz)
                .map(|(m, z)| mz(m, z)),
            experimental_mz: query.neutral_mass.zip(charge_for_mz).map(|(m, z)| mz(m, z)),
            pass_threshold: None,
            scores: hit.scores,
            evidence: hit.evidence,
        });
        Ok(())
    }
}
