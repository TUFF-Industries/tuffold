//! MaxQuant `txt` output → [`ImportBatch`].
//!
//! Reads `msms.txt` (one row per identified MS/MS scan), or `evidence.txt`
//! when that's all there is (its rows without an MS/MS scan, transferred
//! by match-between-runs, are skipped and counted), with
//! `parameters.txt` for the version and fixed modifications and
//! `summary.txt` for which experiment each raw file belongs to, when they
//! sit beside it. Mapping conventions:
//! - Each raw file is an MS run; its MaxQuant experiment is its biosample
//!   (the raw file itself when `summary.txt` is absent).
//! - Rows of one scan are ranked by Andromeda score (ties share a rank);
//!   MaxQuant's PEP and q-values are not imported: Tuffold validates.
//! - Modifications come from `Modified sequence` (`_(Acetyl (Protein
//!   N-term))M(Oxidation (M))PEPTIDE_`, or the old `(ox)` style) through a
//!   table of MaxQuant's names; fixed modifications from `parameters.txt`
//!   apply to every matching residue. Each PSM's mass must then equal
//!   MaxQuant's `Mass`; one that doesn't (an unknown label, say) is left
//!   out and counted, and more than [`MAX_MASS_FAILURES`] of them refuses
//!   the file.
//! - Reverse hits carry MaxQuant's `REV__` accessions, contaminants
//!   `CON__` (targets).

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use tuffold_core::chemistry::{PROTON, peptide_mass};
use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, Spectrum, SpectrumKey,
};

/// MaxQuant reports `Mass` to 5+ decimals; a missing modification is off
/// by at least a dalton.
const MASS_TOLERANCE: f64 = 0.01;
/// Share of PSMs allowed to fail the mass check before the file is refused.
pub const MAX_MASS_FAILURES: f64 = 0.05;

#[derive(Debug, thiserror::Error)]
pub enum MaxQuantError {
    #[error("{file}: {source}")]
    Io {
        file: String,
        #[source]
        source: io::Error,
    },
    #[error("{0}: no msms.txt or evidence.txt here (a MaxQuant txt folder holds them)")]
    NoMsms(String),
    #[error("{file}: missing column {column:?} (columns: {found})")]
    MissingColumn {
        file: String,
        column: &'static str,
        found: String,
    },
    #[error("{file}: line {line}: {problem}")]
    Invalid {
        file: String,
        line: usize,
        problem: String,
    },
    #[error(
        "{file}: {failures} of {checked} PSMs' masses don't match MaxQuant's (e.g. {example}); \
         labels or modifications Tuffold doesn't know are in use, so the file was not imported"
    )]
    Masses {
        file: String,
        failures: usize,
        checked: usize,
        example: String,
    },
}

/// The batch and notes on what was left out.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MaxQuantImport {
    pub batch: ImportBatch,
    pub notes: Vec<String>,
}

/// A MaxQuant modification: mass and where it sits.
#[derive(Debug, Clone, Copy, PartialEq)]
struct KnownModification {
    delta: f64,
    site: SiteKind,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum SiteKind {
    Residue,
    NTerm,
    CTerm,
}

/// MaxQuant's names (and old two-letter abbreviations) for common
/// modifications; monoisotopic deltas from Unimod.
fn known_modification(name: &str) -> Option<KnownModification> {
    let residue = |delta| KnownModification {
        delta,
        site: SiteKind::Residue,
    };
    let n_term = |delta| KnownModification {
        delta,
        site: SiteKind::NTerm,
    };
    let c_term = |delta| KnownModification {
        delta,
        site: SiteKind::CTerm,
    };
    let base = name.split(" (").next().unwrap_or(name).trim();
    let terminal = name.contains("N-term");
    Some(match (base, terminal) {
        ("Oxidation" | "ox", _) => residue(15.994_915),
        ("Acetyl" | "ac", true) => n_term(42.010_565),
        ("ac", false) => n_term(42.010_565),
        ("Acetyl", false) => residue(42.010_565),
        ("Phospho" | "ph", _) => residue(79.966_331),
        ("Carbamidomethyl" | "cam", _) => residue(57.021_464),
        ("Deamidation" | "de", _) => residue(0.984_016),
        ("Gln->pyro-Glu" | "gl" | "pq", _) => n_term(-17.026_549),
        ("Glu->pyro-Glu" | "pe", _) => n_term(-18.010_565),
        ("GlyGly" | "gg", _) => residue(114.042_927),
        ("Methyl" | "me", _) => residue(14.015_650),
        ("Dimethyl" | "di", _) => residue(28.031_300),
        ("Trimethyl" | "tr", _) => residue(42.046_950),
        ("Carbamyl", true) => n_term(43.005_814),
        ("Carbamyl", false) => residue(43.005_814),
        ("Nitro", _) => residue(44.985_078),
        ("Amidated", _) if name.contains("C-term") => c_term(-0.984_016),
        _ => return None,
    })
}

/// The `txt` folder and its PSM table (`msms.txt`, else `evidence.txt`).
fn locate(path: &Path) -> Result<(PathBuf, PathBuf), MaxQuantError> {
    let msms = if path.is_dir() {
        let msms = path.join("msms.txt");
        if msms.is_file() {
            msms
        } else {
            path.join("evidence.txt")
        }
    } else {
        path.to_path_buf()
    };
    if !msms.is_file() {
        return Err(MaxQuantError::NoMsms(path.display().to_string()));
    }
    let folder = msms.parent().map(Path::to_path_buf).unwrap_or_default();
    Ok((folder, msms))
}

/// Whether a folder looks like MaxQuant output.
pub fn is_maxquant_folder(path: &Path) -> bool {
    path.is_dir() && (path.join("msms.txt").is_file() || path.join("evidence.txt").is_file())
}

pub fn import_maxquant(path: &Path) -> Result<MaxQuantImport, MaxQuantError> {
    let (folder, msms) = locate(path)?;
    let parameters = key_values(&folder.join("parameters.txt"))?;
    let experiments = experiments(&folder.join("summary.txt"))?;
    let fixed = fixed_modifications(&parameters);
    let mut notes = Notes::default();
    for name in &fixed.unknown {
        notes.unknown_fixed.push(name.clone());
    }

    let file = msms.display().to_string();
    let io_error = |source| MaxQuantError::Io {
        file: file.clone(),
        source,
    };
    let reader = BufReader::new(File::open(&msms).map_err(io_error)?);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .transpose()
        .map_err(io_error)?
        .unwrap_or_default();
    let columns = Columns::new(&header, &file)?;

    // (run, scan) → rows, so each scan's rows can be ranked.
    let mut scans: BTreeMap<(String, u32), Vec<Row>> = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let line = line.map_err(io_error)?;
        if line.trim().is_empty() {
            continue;
        }
        let number = index + 2;
        let row = columns.row(&line, &file, number, &fixed, &mut notes)?;
        if let Some(row) = row {
            scans
                .entry((row.run.clone(), row.scan))
                .or_default()
                .push(row);
        }
    }
    if notes.checked > 0
        && notes.mass_failures.len() as f64 > MAX_MASS_FAILURES * notes.checked as f64
    {
        return Err(MaxQuantError::Masses {
            file,
            failures: notes.mass_failures.len(),
            checked: notes.checked,
            example: notes.mass_failures.first().cloned().unwrap_or_default(),
        });
    }
    let version = parameters.get("Version").cloned();
    Ok(MaxQuantImport {
        batch: batch(scans, &experiments, version),
        notes: notes.lines(),
    })
}

/// `parameters.txt`: two tab-separated columns.
fn key_values(path: &Path) -> Result<HashMap<String, String>, MaxQuantError> {
    let Ok(file) = File::open(path) else {
        return Ok(HashMap::new());
    };
    let mut values = HashMap::new();
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|source| MaxQuantError::Io {
            file: path.display().to_string(),
            source,
        })?;
        if let Some((key, value)) = line.split_once('\t') {
            values.insert(key.trim().to_string(), value.trim().to_string());
        }
    }
    Ok(values)
}

/// `summary.txt`: raw file → experiment.
fn experiments(path: &Path) -> Result<HashMap<String, String>, MaxQuantError> {
    let Ok(file) = File::open(path) else {
        return Ok(HashMap::new());
    };
    let file_name = path.display().to_string();
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()
        .map_err(|source| MaxQuantError::Io {
            file: file_name.clone(),
            source,
        })?
        .unwrap_or_default();
    let names: Vec<&str> = header.split('\t').map(str::trim).collect();
    let (Some(raw), Some(experiment)) = (
        names.iter().position(|name| *name == "Raw file"),
        names.iter().position(|name| *name == "Experiment"),
    ) else {
        return Ok(HashMap::new());
    };
    let mut out = HashMap::new();
    for line in lines {
        let line = line.map_err(|source| MaxQuantError::Io {
            file: file_name.clone(),
            source,
        })?;
        let fields: Vec<&str> = line.split('\t').collect();
        if let (Some(raw), Some(experiment)) = (fields.get(raw), fields.get(experiment)) {
            let (raw, experiment) = (raw.trim(), experiment.trim());
            if !raw.is_empty() && !experiment.is_empty() && raw != "Total" {
                out.insert(raw.to_string(), experiment.to_string());
            }
        }
    }
    Ok(out)
}

struct Fixed {
    /// Residue → (name, delta).
    residues: Vec<(char, String, f64)>,
    unknown: Vec<String>,
}

/// Fixed modifications from `parameters.txt`, e.g. `Carbamidomethyl (C)`.
fn fixed_modifications(parameters: &HashMap<String, String>) -> Fixed {
    let mut fixed = Fixed {
        residues: Vec::new(),
        unknown: Vec::new(),
    };
    let Some(value) = parameters.get("Fixed modifications") else {
        return fixed;
    };
    for name in value
        .split(';')
        .map(str::trim)
        .filter(|name| !name.is_empty())
    {
        let residues: Vec<char> = name
            .rsplit_once(" (")
            .map(|(_, inside)| inside.trim_end_matches(')'))
            .filter(|inside| inside.chars().all(|ch| ch.is_ascii_uppercase()))
            .map(|inside| inside.chars().collect())
            .unwrap_or_default();
        match known_modification(name) {
            Some(known) if known.site == SiteKind::Residue && !residues.is_empty() => {
                let base = name.split(" (").next().unwrap_or(name).to_string();
                fixed.residues.extend(
                    residues
                        .into_iter()
                        .map(|residue| (residue, base.clone(), known.delta)),
                );
            }
            _ => fixed.unknown.push(name.to_string()),
        }
    }
    fixed
}

struct Row {
    run: String,
    /// The `Experiment` column of evidence.txt.
    experiment: Option<String>,
    scan: u32,
    score: f64,
    psm: Psm,
    rt_seconds: Option<f64>,
}

struct Columns {
    raw: usize,
    scan: usize,
    sequence: usize,
    modified: usize,
    charge: usize,
    score: usize,
    proteins: Option<usize>,
    /// evidence.txt leaves `Proteins` empty on reverse hits; this has them.
    leading: Option<usize>,
    mz: Option<usize>,
    mass: Option<usize>,
    rt: Option<usize>,
    delta: Option<usize>,
    reverse: Option<usize>,
    experiment: Option<usize>,
}

impl Columns {
    fn new(header: &str, file: &str) -> Result<Self, MaxQuantError> {
        let names: Vec<&str> = header
            .trim_start_matches('\u{feff}')
            .split('\t')
            .map(str::trim)
            .collect();
        let find = |name: &str| {
            names
                .iter()
                .position(|known| known.eq_ignore_ascii_case(name))
        };
        let required = |column: &'static str| {
            find(column).ok_or_else(|| MaxQuantError::MissingColumn {
                file: file.to_string(),
                column,
                found: names.join(", "),
            })
        };
        Ok(Self {
            raw: required("Raw file")?,
            scan: find("Scan number")
                .or_else(|| find("MS/MS scan number"))
                .ok_or_else(|| MaxQuantError::MissingColumn {
                    file: file.to_string(),
                    column: "Scan number",
                    found: names.join(", "),
                })?,
            sequence: required("Sequence")?,
            modified: required("Modified sequence")?,
            charge: required("Charge")?,
            score: required("Score")?,
            proteins: find("Proteins"),
            leading: find("Leading proteins"),
            mz: find("m/z"),
            mass: find("Mass"),
            rt: find("Retention time"),
            delta: find("Delta score"),
            reverse: find("Reverse"),
            experiment: find("Experiment"),
        })
    }

    fn row(
        &self,
        line: &str,
        file: &str,
        number: usize,
        fixed: &Fixed,
        notes: &mut Notes,
    ) -> Result<Option<Row>, MaxQuantError> {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |index: usize| fields.get(index).map_or("", |field| field.trim());
        let optional = |index: Option<usize>| index.map(get).filter(|value| !value.is_empty());
        let number_of = |value: &str| value.parse::<f64>().ok().filter(|value| value.is_finite());
        let invalid = |problem: String| MaxQuantError::Invalid {
            file: file.to_string(),
            line: number,
            problem,
        };
        let run = get(self.raw).to_string();
        if get(self.scan).is_empty() {
            notes.without_scan += 1; // a match-between-runs transfer
            return Ok(None);
        }
        let scan: u32 = get(self.scan)
            .parse()
            .map_err(|_| invalid(format!("scan number {:?}", get(self.scan))))?;
        let charge: i16 = get(self.charge)
            .parse()
            .ok()
            .filter(|&charge| charge > 0)
            .ok_or_else(|| invalid(format!("charge {:?}", get(self.charge))))?;
        let Some(score) = number_of(get(self.score)) else {
            notes.unscored += 1;
            return Ok(None);
        };
        let sequence = get(self.sequence).to_ascii_uppercase();
        let peptide = match parse_modified(get(self.modified), &sequence, fixed) {
            Ok(peptide) => peptide,
            Err(Unparsed::Unknown(name)) => {
                *notes.unknown_modifications.entry(name).or_default() += 1;
                return Ok(None);
            }
            Err(Unparsed::Malformed(problem)) => return Err(invalid(problem)),
        };
        let mass = peptide_mass(&peptide.sequence).map(|mass| {
            mass + peptide
                .modifications
                .iter()
                .map(|modification| modification.mass_delta)
                .sum::<f64>()
        });
        if let (Some(mass), Some(reported)) = (mass, optional(self.mass).and_then(number_of)) {
            notes.checked += 1;
            if (mass - reported).abs() > MASS_TOLERANCE {
                notes.mass_failures.push(format!(
                    "{} ({:+.4} Da)",
                    peptide.proforma(),
                    reported - mass
                ));
                return Ok(None);
            }
        }
        let accessions: Vec<String> = optional(self.proteins)
            .or_else(|| optional(self.leading))
            .map(|proteins| {
                proteins
                    .split(';')
                    .map(str::trim)
                    .filter(|accession| !accession.is_empty())
                    .map(str::to_string)
                    .collect()
            })
            .unwrap_or_default();
        if optional(self.reverse) == Some("+") && !accessions.iter().any(|a| is_decoy_accession(a))
        {
            notes.reverse_without_label += 1;
        }
        let mut scores = vec![Score {
            name: "Andromeda:score".into(),
            accession: None,
            value: score,
        }];
        if let Some(delta) = optional(self.delta).and_then(number_of) {
            scores.push(Score {
                name: "MaxQuant:delta score".into(),
                accession: None,
                value: delta,
            });
        }
        let z = f64::from(charge);
        Ok(Some(Row {
            experiment: optional(self.experiment).map(str::to_string),
            run: run.clone(),
            scan,
            score,
            rt_seconds: optional(self.rt)
                .and_then(number_of)
                .map(|minutes| minutes * 60.0),
            psm: Psm {
                ms_run: run,
                spectrum: SpectrumKey::Scan(scan),
                rank: 1,
                charge,
                calculated_mz: mass.map(|mass| (mass + z * PROTON) / z),
                experimental_mz: optional(self.mz).and_then(number_of),
                pass_threshold: None,
                scores,
                evidence: accessions
                    .into_iter()
                    .map(|accession| Evidence {
                        accession,
                        start: None,
                        end: None,
                        pre: None,
                        post: None,
                    })
                    .collect(),
                peptide,
            },
        }))
    }
}

enum Unparsed {
    Unknown(String),
    Malformed(String),
}

/// `_(Acetyl (Protein N-term))M(Oxidation (M))PEPTIDE_` → a peptide form,
/// with fixed modifications added on their residues.
fn parse_modified(modified: &str, sequence: &str, fixed: &Fixed) -> Result<PeptideForm, Unparsed> {
    let body = modified.trim().trim_matches('_');
    let mut residues = String::new();
    let mut modifications = Vec::new();
    let mut chars = body.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '(' {
            let mut depth = 1;
            let mut name = String::new();
            for inner in chars.by_ref() {
                match inner {
                    '(' => depth += 1,
                    ')' => depth -= 1,
                    _ => {}
                }
                if depth == 0 {
                    break;
                }
                name.push(inner);
            }
            if depth != 0 {
                return Err(Unparsed::Malformed(format!(
                    "unbalanced parentheses in {modified:?}"
                )));
            }
            let known = known_modification(&name).ok_or_else(|| Unparsed::Unknown(name.clone()))?;
            let position = u32::try_from(residues.len()).unwrap_or(u32::MAX);
            let site = match known.site {
                SiteKind::NTerm => ModSite::NTerm,
                SiteKind::CTerm => ModSite::CTerm,
                SiteKind::Residue if position == 0 => ModSite::NTerm,
                SiteKind::Residue => ModSite::Residue(position),
            };
            modifications.push(Modification {
                site,
                mass_delta: known.delta,
                name: Some(name.split(" (").next().unwrap_or(&name).to_string()),
                accession: None,
            });
        } else if ch.is_ascii_alphabetic() {
            residues.push(ch.to_ascii_uppercase());
        } else {
            return Err(Unparsed::Malformed(format!(
                "unexpected {ch:?} in {modified:?}"
            )));
        }
    }
    if residues != sequence {
        return Err(Unparsed::Malformed(format!(
            "modified sequence {modified:?} doesn't spell {sequence}"
        )));
    }
    for (index, residue) in residues.chars().enumerate() {
        for (fixed_residue, name, delta) in &fixed.residues {
            let position = u32::try_from(index + 1).unwrap_or(u32::MAX);
            let already = modifications
                .iter()
                .any(|m: &Modification| m.site == ModSite::Residue(position));
            if residue == *fixed_residue && !already {
                modifications.push(Modification {
                    site: ModSite::Residue(position),
                    mass_delta: *delta,
                    name: Some(name.clone()),
                    accession: None,
                });
            }
        }
    }
    PeptideForm::new(&residues, modifications)
        .map_err(|error| Unparsed::Malformed(error.to_string()))
}

fn batch(
    scans: BTreeMap<(String, u32), Vec<Row>>,
    experiments: &HashMap<String, String>,
    version: Option<String>,
) -> ImportBatch {
    let mut batch = ImportBatch {
        software: Software {
            name: "MaxQuant".into(),
            version,
        },
        ..ImportBatch::default()
    };
    let mut runs: BTreeMap<String, String> = BTreeMap::new();
    let mut proteins: BTreeMap<String, Protein> = BTreeMap::new();
    for ((run, scan), rows) in scans {
        let biosample = rows
            .iter()
            .find_map(|row| row.experiment.clone())
            .or_else(|| experiments.get(&run).cloned())
            .unwrap_or_else(|| run.clone());
        runs.entry(run.clone()).or_insert(biosample);
        let best = rows
            .iter()
            .max_by(|a, b| a.score.total_cmp(&b.score))
            .expect("a row");
        batch.spectra.push(Spectrum {
            ms_run: run.clone(),
            key: SpectrumKey::Scan(scan),
            source_id: None,
            precursor_mz: best.psm.experimental_mz,
            precursor_charge: Some(best.psm.charge),
            rt_seconds: best.rt_seconds,
        });
        let scores: Vec<f64> = rows.iter().map(|row| row.score).collect();
        for mut row in rows {
            let better = scores.iter().filter(|&&score| score > row.score).count();
            row.psm.rank = u32::try_from(better + 1).unwrap_or(u32::MAX);
            for evidence in &row.psm.evidence {
                proteins
                    .entry(evidence.accession.clone())
                    .or_insert_with(|| Protein {
                        accession: evidence.accession.clone(),
                        description: None,
                        sequence: None,
                        is_decoy: is_decoy_accession(&evidence.accession),
                        decoy_explicit: false,
                    });
            }
            batch.psms.push(row.psm);
        }
    }
    let mut biosamples: Vec<BioSample> = Vec::new();
    for (run, biosample) in &runs {
        if !biosamples.iter().any(|known| &known.name == biosample) {
            biosamples.push(BioSample {
                name: biosample.clone(),
                category: None,
            });
        }
        batch.ms_runs.push(MsRun {
            name: run.clone(),
            biosample: biosample.clone(),
            spectra_file: None,
        });
    }
    batch.biosamples = biosamples;
    batch.proteins = proteins.into_values().collect();
    batch
}

#[derive(Default)]
struct Notes {
    checked: usize,
    mass_failures: Vec<String>,
    unknown_modifications: BTreeMap<String, usize>,
    unknown_fixed: Vec<String>,
    unscored: usize,
    reverse_without_label: usize,
    without_scan: usize,
}

impl Notes {
    fn lines(&self) -> Vec<String> {
        let mut lines = vec![
            "MaxQuant's PEP and FDR filtering were not imported: Tuffold validates the \
             Andromeda scores itself. msms.txt usually holds only PSMs passing MaxQuant's \
             PSM FDR; set it to 100% in MaxQuant to give Tuffold every PSM."
                .to_string(),
        ];
        for (name, count) in &self.unknown_modifications {
            lines.push(format!(
                "{count} PSMs carry {name:?}, a modification Tuffold doesn't know yet; left out."
            ));
        }
        if !self.unknown_fixed.is_empty() {
            lines.push(format!(
                "Fixed modifications Tuffold doesn't know: {}; PSMs whose mass then disagrees \
                 are left out.",
                self.unknown_fixed.join(", ")
            ));
        }
        if !self.mass_failures.is_empty() {
            lines.push(format!(
                "{} PSMs whose mass disagrees with MaxQuant's were left out, e.g. {}.",
                self.mass_failures.len(),
                self.mass_failures
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        if self.without_scan > 0 {
            lines.push(format!(
                "{} evidence rows without an MS/MS scan (match-between-runs transfers) were left \
                 out: they have no spectrum to validate.",
                self.without_scan
            ));
        }
        if self.unscored > 0 {
            lines.push(format!(
                "{} PSMs without a score were left out.",
                self.unscored
            ));
        }
        if self.reverse_without_label > 0 {
            lines.push(format!(
                "{} reverse hits have no REV__ accession, so they import as targets.",
                self.reverse_without_label
            ));
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_fixed() -> Fixed {
        Fixed {
            residues: Vec::new(),
            unknown: Vec::new(),
        }
    }

    #[test]
    fn modified_sequences_parse_in_both_styles() {
        let new = parse_modified(
            "_(Acetyl (Protein N-term))M(Oxidation (M))PEPTIDE_",
            "MPEPTIDE",
            &no_fixed(),
        )
        .ok()
        .unwrap();
        assert_eq!(new.proforma(), "[+42.0106]-M[+15.9949]PEPTIDE");
        let old = parse_modified("_(ac)M(ox)PEPTIDE_", "MPEPTIDE", &no_fixed())
            .ok()
            .unwrap();
        assert_eq!(old.proforma(), new.proforma());
    }

    #[test]
    fn fixed_modifications_land_on_every_matching_residue() {
        let parameters = HashMap::from([(
            "Fixed modifications".to_string(),
            "Carbamidomethyl (C); Mystery (W)".to_string(),
        )]);
        let fixed = fixed_modifications(&parameters);
        assert_eq!(fixed.unknown, ["Mystery (W)"]);
        let peptide = parse_modified("_ACDCK_", "ACDCK", &fixed).ok().unwrap();
        assert_eq!(peptide.proforma(), "AC[+57.0215]DC[+57.0215]K");
    }

    #[test]
    fn unknown_and_malformed_modifications_are_told_apart() {
        assert!(matches!(
            parse_modified("_PEP(Sulfo (Y))TIDE_", "PEPTIDE", &no_fixed()),
            Err(Unparsed::Unknown(name)) if name == "Sulfo (Y)"
        ));
        assert!(matches!(
            parse_modified("_PEP(ox_", "PEPTIDE", &no_fixed()),
            Err(Unparsed::Malformed(_))
        ));
        assert!(matches!(
            parse_modified("_PEPTIDEK_", "PEPTIDE", &no_fixed()),
            Err(Unparsed::Malformed(_))
        ));
    }
}
