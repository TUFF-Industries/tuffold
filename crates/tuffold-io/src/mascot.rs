//! Mascot results files (`.dat`, MIME multipart text) → [`ImportBatch`].
//!
//! Sections read: `parameters` (searched file, mass type), `masses`
//! (variable `deltaN` and fixed `FixedModN` modifications), `summary` and
//! `decoy_summary` (`qexpN` precursor, `qmatchN` candidate count,
//! `qplugholeN` homology threshold), `peptides` and `decoy_peptides`
//! (`qN_pM` matches, rank M), `proteins` (descriptions), and `queryN`
//! (title, retention time, scans). Mapping conventions:
//! - The searched file (`FILE=`) is the MS run and its biosample.
//! - A match's modifications: fixed ones on every listed residue, variable
//!   ones from its position string (N-terminus, each residue, C-terminus).
//!   Its mass must equal Mascot's peptide Mr; more than
//!   [`MAX_MASS_FAILURES`] of matches failing refuses the file, fewer are
//!   left out and counted.
//! - Scores: ion score, homology threshold, and the expectation value
//!   Mascot reports (candidate count × 10^(−score/10)). Decoy matches are
//!   imported against `DECOY_<accession>`, flagged explicitly: Mascot lists
//!   them under their target's accession.
//! - Spectrum keys: `scans=` when given, else a scan number in the title,
//!   else the title.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use tuffold_core::chemistry::{PROTON, peptide_mass};
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, Spectrum, SpectrumKey,
};

use crate::xml::run_name;

const MASS_TOLERANCE: f64 = 0.01;
/// Share of matches allowed to fail the mass check before the file is refused.
pub const MAX_MASS_FAILURES: f64 = 0.05;
const DECOY_PREFIX: &str = "DECOY_";

#[derive(Debug, thiserror::Error)]
pub enum MascotError {
    #[error("{file}: {source}")]
    Io {
        file: String,
        #[source]
        source: io::Error,
    },
    #[error("{file}: not a Mascot results file ({problem})")]
    NotMascot { file: String, problem: String },
    #[error("{file}: line {line}: {problem}")]
    Invalid {
        file: String,
        line: usize,
        problem: String,
    },
    #[error(
        "{file}: {failures} of {checked} matches' masses don't match Mascot's (e.g. {example}), \
         so the file was not imported rather than misread"
    )]
    Masses {
        file: String,
        failures: usize,
        checked: usize,
        example: String,
    },
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MascotImport {
    pub batch: ImportBatch,
    pub notes: Vec<String>,
}

/// The residues before and after a peptide in one protein.
type Flanks = (Option<char>, Option<char>);

#[derive(Default)]
struct Query {
    title: Option<String>,
    rt_seconds: Option<f64>,
    scans: Option<String>,
    mz: Option<f64>,
    charge: Option<i16>,
    candidates: Option<f64>,
    homology: Option<f64>,
    decoy_candidates: Option<f64>,
    decoy_homology: Option<f64>,
}

struct Draft {
    query: u32,
    rank: u32,
    decoy: bool,
    peptide: PeptideForm,
    mass: f64,
    score: f64,
    evidence: Vec<Evidence>,
}

#[derive(Default)]
struct State {
    file: Option<String>,
    version: Option<String>,
    average: bool,
    variable: HashMap<u32, (String, f64)>,
    fixed_mass: HashMap<u32, (String, f64)>,
    fixed_residues: HashMap<u32, String>,
    queries: BTreeMap<u32, Query>,
    descriptions: HashMap<String, String>,
    drafts: Vec<Draft>,
    /// (decoy?, query, rank) → residues before and after, per protein.
    terms: HashMap<(bool, u32, u32), Vec<Flanks>>,
    notes: Notes,
}

pub fn import_mascot(path: &Path) -> Result<MascotImport, MascotError> {
    let file = path.display().to_string();
    let io_error = |source| MascotError::Io {
        file: file.clone(),
        source,
    };
    let reader = BufReader::new(File::open(path).map_err(io_error)?);
    let mut state = State::default();
    let mut boundary: Option<String> = None;
    let mut section = String::new();
    for (index, line) in reader.lines().enumerate() {
        let line = line.map_err(io_error)?;
        let line = line.trim_end_matches('\r');
        let number = index + 1;
        if boundary.is_none() {
            if let Some(found) = line.split("boundary=").nth(1) {
                boundary = Some(format!("--{}", found.trim().trim_matches('"')));
            } else if number > 20 {
                break;
            }
            continue;
        }
        let marker = boundary.as_deref().expect("set above");
        if line.starts_with(marker) {
            section.clear();
            continue;
        }
        if let Some(name) = line
            .strip_prefix("Content-Type: application/x-Mascot; name=")
            .or_else(|| line.strip_prefix("Content-Type: application/x-mascot; name="))
        {
            section = name.trim().trim_matches('"').to_string();
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        state
            .line(&section, key, value)
            .map_err(|problem| MascotError::Invalid {
                file: file.clone(),
                line: number,
                problem,
            })?;
    }
    if boundary.is_none() || state.queries.is_empty() && state.drafts.is_empty() {
        return Err(MascotError::NotMascot {
            file,
            problem: "no MIME sections or queries".into(),
        });
    }
    if state.average {
        return Err(MascotError::NotMascot {
            file,
            problem: "searched with average masses, which Tuffold doesn't import".into(),
        });
    }
    let checked = state.notes.checked;
    let failures = state.notes.mass_failures.len();
    if checked > 0 && failures as f64 > MAX_MASS_FAILURES * checked as f64 {
        return Err(MascotError::Masses {
            file,
            failures,
            checked,
            example: state.notes.mass_failures[0].clone(),
        });
    }
    Ok(state.finish(&file))
}

impl State {
    fn line(&mut self, section: &str, key: &str, value: &str) -> Result<(), String> {
        match section {
            "parameters" => match key {
                "FILE" => self.file = Some(value.trim().to_string()),
                "MASS" => self.average = value.trim().eq_ignore_ascii_case("average"),
                _ => {}
            },
            "header" if key == "version" => self.version = Some(value.trim().to_string()),
            "masses" => self.mass_line(key, value)?,
            "summary" | "decoy_summary" => {
                let decoy = section == "decoy_summary";
                self.summary_line(key, value, decoy);
            }
            "peptides" | "decoy_peptides" => {
                self.peptide_line(key, value, section == "decoy_peptides")?
            }
            "proteins" => {
                let accession = key.trim().trim_matches('"').to_string();
                if let Some((_, description)) = value.split_once(',') {
                    self.descriptions
                        .insert(accession, description.trim().trim_matches('"').to_string());
                }
            }
            _ => {
                if let Some(query) = section.strip_prefix("query").and_then(|n| n.parse().ok()) {
                    let entry = self.queries.entry(query).or_default();
                    match key {
                        "title" => entry.title = Some(url_decode(value)),
                        "rtinseconds" => {
                            entry.rt_seconds =
                                value.split(',').next().and_then(|v| v.trim().parse().ok())
                        }
                        "scans" => entry.scans = Some(value.trim().to_string()),
                        _ => {}
                    }
                }
            }
        }
        Ok(())
    }

    fn mass_line(&mut self, key: &str, value: &str) -> Result<(), String> {
        let numbered = |prefix: &str| key.strip_prefix(prefix).and_then(|n| n.parse::<u32>().ok());
        let named = |value: &str| -> Result<(String, f64), String> {
            let (mass, name) = value.split_once(',').unwrap_or((value, ""));
            let mass: f64 = mass
                .trim()
                .parse()
                .map_err(|_| format!("{key}: mass {mass:?}"))?;
            Ok((name.trim().to_string(), mass))
        };
        if let Some(n) = numbered("delta") {
            self.variable.insert(n, named(value)?);
        } else if let Some(n) = numbered("FixedModResidues") {
            self.fixed_residues.insert(n, value.trim().to_string());
        } else if let Some(n) = numbered("FixedMod") {
            self.fixed_mass.insert(n, named(value)?);
        }
        Ok(())
    }

    fn summary_line(&mut self, key: &str, value: &str, decoy: bool) {
        let field = |prefix: &str| key.strip_prefix(prefix).and_then(|n| n.parse::<u32>().ok());
        if let Some(n) = field("qexp") {
            let entry = self.queries.entry(n).or_default();
            let mut parts = value.split(',');
            entry.mz = parts.next().and_then(|mz| mz.trim().parse().ok());
            entry.charge = parts.next().and_then(parse_charge);
        } else if let Some(n) = field("qmatch") {
            let candidates = value.trim().parse().ok();
            let entry = self.queries.entry(n).or_default();
            if decoy {
                entry.decoy_candidates = candidates;
            } else {
                entry.candidates = candidates;
            }
        } else if let Some(n) = field("qplughole") {
            let homology = value.trim().parse().ok();
            let entry = self.queries.entry(n).or_default();
            if decoy {
                entry.decoy_homology = homology;
            } else {
                entry.homology = homology;
            }
        }
    }

    fn peptide_line(&mut self, key: &str, value: &str, decoy: bool) -> Result<(), String> {
        // qN_pM, qN_pM_terms; other qN_pM_* lines (error-tolerant extras) are counted.
        let Some(rest) = key.strip_prefix('q') else {
            return Ok(());
        };
        let mut parts = rest.split('_');
        let (Some(query), Some(rank)) = (
            parts.next().and_then(|n| n.parse::<u32>().ok()),
            parts
                .next()
                .and_then(|p| p.strip_prefix('p'))
                .and_then(|n| n.parse::<u32>().ok()),
        ) else {
            return Ok(());
        };
        match parts.next() {
            None => self.match_line(query, rank, value, decoy),
            Some("terms") => {
                let flanks = value
                    .split(':')
                    .map(|pair| {
                        let mut ends = pair.split(',');
                        let flank = |end: Option<&str>| {
                            end.and_then(|e| e.trim().chars().next())
                                .filter(|ch| ch.is_ascii_alphabetic() || *ch == '-')
                        };
                        (flank(ends.next()), flank(ends.next()))
                    })
                    .collect();
                self.terms.insert((decoy, query, rank), flanks);
                Ok(())
            }
            Some("et") | Some("subst") => {
                self.notes.error_tolerant += 1;
                Ok(())
            }
            Some(_) => Ok(()),
        }
    }

    fn match_line(
        &mut self,
        query: u32,
        rank: u32,
        value: &str,
        decoy: bool,
    ) -> Result<(), String> {
        if value.trim() == "-1" {
            return Ok(()); // no match
        }
        let (fields, proteins) = value.split_once(';').unwrap_or((value, ""));
        let fields: Vec<&str> = fields.split(',').collect();
        if fields.len() < 8 {
            return Err(format!(
                "q{query}_p{rank}: {} fields, expected 11",
                fields.len()
            ));
        }
        let number = |index: usize| -> Result<f64, String> {
            fields[index]
                .trim()
                .parse()
                .map_err(|_| format!("q{query}_p{rank}: field {index} is {:?}", fields[index]))
        };
        let mascot_mass = number(1)?;
        let sequence = fields[4].trim().to_ascii_uppercase();
        let score = number(7)?;
        let peptide = self
            .peptide(&sequence, fields[6].trim())
            .map_err(|problem| format!("q{query}_p{rank}: {problem}"))?;
        let Some(mass) = peptide_mass(&peptide.sequence).map(|mass| {
            mass + peptide
                .modifications
                .iter()
                .map(|m| m.mass_delta)
                .sum::<f64>()
        }) else {
            self.notes.unchecked += 1;
            return Ok(());
        };
        self.notes.checked += 1;
        if (mass - mascot_mass).abs() > MASS_TOLERANCE {
            self.notes.mass_failures.push(format!(
                "{} ({:+.4} Da)",
                peptide.proforma(),
                mascot_mass - mass
            ));
            return Ok(());
        }
        let evidence = proteins
            .split(',')
            .filter_map(|entry| {
                let accession = entry.split(':').next()?.trim().trim_matches('"');
                (!accession.is_empty()).then(|| Evidence {
                    accession: if decoy {
                        format!("{DECOY_PREFIX}{accession}")
                    } else {
                        accession.to_string()
                    },
                    start: None,
                    end: None,
                    pre: None,
                    post: None,
                })
            })
            .collect();
        self.drafts.push(Draft {
            query,
            rank,
            decoy,
            peptide,
            mass,
            score,
            evidence,
        });
        Ok(())
    }

    /// Fixed modifications on their residues, variable ones from the
    /// position string: N-terminus, each residue, C-terminus.
    fn peptide(&self, sequence: &str, positions: &str) -> Result<PeptideForm, String> {
        let length = sequence.len();
        let mut modifications = Vec::new();
        let digits: Vec<char> = positions.chars().collect();
        if !digits.is_empty() && digits.len() != length + 2 {
            return Err(format!(
                "modification string {positions:?} doesn't fit {sequence}"
            ));
        }
        for (index, digit) in digits.iter().enumerate() {
            let id = match digit {
                '0' => continue,
                '1'..='9' => u32::from(*digit) - u32::from('0'),
                'A'..='W' => u32::from(*digit) - u32::from('A') + 10,
                other => return Err(format!("modification code {other:?}")),
            };
            let (name, delta) = self
                .variable
                .get(&id)
                .ok_or_else(|| format!("variable modification {id} isn't listed"))?;
            let site = match index {
                0 => ModSite::NTerm,
                i if i == length + 1 => ModSite::CTerm,
                i => ModSite::Residue(u32::try_from(i).map_err(|_| "long peptide".to_string())?),
            };
            modifications.push(Modification {
                site,
                mass_delta: *delta,
                name: Some(short_name(name)).filter(|n| !n.is_empty()),
                accession: None,
            });
        }
        for (n, (name, delta)) in &self.fixed_mass {
            let residues = self.fixed_residues.get(n).map_or("", String::as_str);
            let mut sites = Vec::new();
            match residues {
                "N_term" | "N-term" => sites.push(ModSite::NTerm),
                "C_term" | "C-term" => sites.push(ModSite::CTerm),
                _ => {
                    for (index, residue) in sequence.chars().enumerate() {
                        if residues.contains(residue) {
                            sites.push(ModSite::Residue(
                                u32::try_from(index + 1).map_err(|_| "long peptide".to_string())?,
                            ));
                        }
                    }
                }
            }
            for site in sites {
                modifications.push(Modification {
                    site,
                    mass_delta: *delta,
                    name: Some(short_name(name)).filter(|n| !n.is_empty()),
                    accession: None,
                });
            }
        }
        PeptideForm::new(sequence, modifications).map_err(|error| error.to_string())
    }

    fn finish(mut self, file: &str) -> MascotImport {
        let run = self
            .file
            .as_deref()
            .map(run_name)
            .filter(|name| !name.is_empty())
            .unwrap_or_else(|| run_name(file));
        let mut batch = ImportBatch {
            software: Software {
                name: "Mascot".into(),
                version: self.version.clone(),
            },
            biosamples: vec![BioSample {
                name: run.clone(),
                category: None,
            }],
            ms_runs: vec![MsRun {
                name: run.clone(),
                biosample: run.clone(),
                spectra_file: self.file.clone(),
            }],
            ..ImportBatch::default()
        };
        let mut spectra: BTreeMap<u32, Spectrum> = BTreeMap::new();
        let mut proteins: BTreeMap<String, Protein> = BTreeMap::new();
        for mut draft in std::mem::take(&mut self.drafts) {
            let query = self.queries.get(&draft.query);
            let Some(charge) = query.and_then(|q| q.charge) else {
                self.notes.without_charge += 1;
                continue;
            };
            let key = spectrum_key(draft.query, query);
            spectra.entry(draft.query).or_insert_with(|| Spectrum {
                ms_run: run.clone(),
                key: key.clone(),
                source_id: query.and_then(|q| q.title.clone()),
                precursor_mz: query.and_then(|q| q.mz),
                precursor_charge: Some(charge),
                rt_seconds: query.and_then(|q| q.rt_seconds),
            });
            if let Some(flanks) = self.terms.get(&(draft.decoy, draft.query, draft.rank)) {
                for (evidence, (pre, post)) in draft.evidence.iter_mut().zip(flanks) {
                    (evidence.pre, evidence.post) = (*pre, *post);
                }
            }
            for evidence in &draft.evidence {
                let target = evidence
                    .accession
                    .strip_prefix(DECOY_PREFIX)
                    .filter(|_| draft.decoy);
                proteins
                    .entry(evidence.accession.clone())
                    .or_insert_with(|| Protein {
                        accession: evidence.accession.clone(),
                        description: self
                            .descriptions
                            .get(target.unwrap_or(&evidence.accession))
                            .filter(|_| !draft.decoy)
                            .cloned(),
                        sequence: None,
                        is_decoy: draft.decoy,
                        decoy_explicit: true,
                    });
            }
            let mut scores = vec![Score {
                name: "Mascot:score".into(),
                accession: None,
                value: draft.score,
            }];
            let (candidates, homology) = match query {
                Some(q) if draft.decoy => (q.decoy_candidates, q.decoy_homology),
                Some(q) => (q.candidates, q.homology),
                None => (None, None),
            };
            if let Some(candidates) = candidates.filter(|c| *c >= 1.0) {
                scores.push(Score {
                    name: "Mascot:expectation value".into(),
                    accession: None,
                    value: candidates * 10f64.powf(-draft.score / 10.0),
                });
            }
            if let Some(homology) = homology {
                scores.push(Score {
                    name: "Mascot:homology threshold".into(),
                    accession: None,
                    value: homology,
                });
            }
            let z = f64::from(charge);
            batch.psms.push(Psm {
                ms_run: run.clone(),
                spectrum: key,
                rank: draft.rank,
                charge,
                calculated_mz: Some((draft.mass + z * PROTON) / z),
                experimental_mz: query.and_then(|q| q.mz),
                pass_threshold: None,
                scores,
                evidence: draft.evidence,
                peptide: draft.peptide,
            });
        }
        batch.spectra = spectra.into_values().collect();
        batch.proteins = proteins.into_values().collect();
        MascotImport {
            batch,
            notes: self.notes.lines(),
        }
    }
}

/// `scans=` (the first of a range), else a scan in the title, else the title.
fn spectrum_key(number: u32, query: Option<&Query>) -> SpectrumKey {
    let scans = query.and_then(|q| q.scans.as_deref());
    if let Some(scan) = scans.and_then(|scans| scans.split(['-', ',']).next()?.trim().parse().ok())
    {
        return SpectrumKey::Scan(scan);
    }
    let title = query
        .and_then(|q| q.title.as_deref())
        .map(str::trim)
        .filter(|t| !t.is_empty());
    let Some(title) = title else {
        return SpectrumKey::Native(format!("query {number}"));
    };
    if let Some((_, rest)) = title.split_once("scan=") {
        let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
        if let Ok(scan) = digits.parse() {
            return SpectrumKey::Scan(scan);
        }
    }
    crate::sf3::key_from_title(Some(title), 0)
}

fn parse_charge(text: &str) -> Option<i16> {
    let text = text.trim();
    let (digits, sign) = match text.strip_suffix('+') {
        Some(digits) => (digits, 1),
        None => (text.strip_suffix('-')?, -1),
    };
    digits
        .parse::<i16>()
        .ok()
        .map(|charge| charge * sign)
        .filter(|&charge| charge > 0)
}

/// `Oxidation (M)` → `Oxidation`.
fn short_name(name: &str) -> String {
    name.split(" (").next().unwrap_or(name).trim().to_string()
}

fn url_decode(text: &str) -> String {
    let bytes = text.trim().as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index] == b'%'
            && let Some(byte) = bytes
                .get(index + 1..index + 3)
                .and_then(|hex| std::str::from_utf8(hex).ok())
                .and_then(|hex| u8::from_str_radix(hex, 16).ok())
        {
            out.push(byte);
            index += 3;
        } else {
            out.push(bytes[index]);
            index += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

#[derive(Default)]
struct Notes {
    checked: usize,
    unchecked: usize,
    mass_failures: Vec<String>,
    without_charge: usize,
    error_tolerant: usize,
}

impl Notes {
    fn lines(&self) -> Vec<String> {
        let mut lines = vec![
            "Mascot's significance thresholds were not applied: Tuffold validates the ion \
             scores itself."
                .to_string(),
        ];
        if !self.mass_failures.is_empty() {
            lines.push(format!(
                "{} matches whose mass disagrees with Mascot's were left out, e.g. {}.",
                self.mass_failures.len(),
                self.mass_failures
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        if self.unchecked > 0 {
            lines.push(format!(
                "{} matches with non-standard residues couldn't be mass-checked and were left out.",
                self.unchecked
            ));
        }
        if self.without_charge > 0 {
            lines.push(format!(
                "{} matches on queries without a charge were left out.",
                self.without_charge
            ));
        }
        if self.error_tolerant > 0 {
            lines.push(format!(
                "{} error-tolerant or substitution details were ignored.",
                self.error_tolerant
            ));
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn charges_titles_and_keys_decode() {
        assert_eq!(parse_charge("2+"), Some(2));
        assert_eq!(parse_charge("1-"), None);
        assert_eq!(url_decode("run%201%2Escan%3D7"), "run 1.scan=7");
        let query = Query {
            title: Some("run.raw NativeID:\"controllerType=0 scan=1234\"".into()),
            ..Query::default()
        };
        assert_eq!(spectrum_key(3, Some(&query)), SpectrumKey::Scan(1234));
        assert_eq!(spectrum_key(3, None), SpectrumKey::Native("query 3".into()));
    }
}
