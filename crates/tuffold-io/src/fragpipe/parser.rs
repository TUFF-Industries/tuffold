//! FragPipe bundle → [`ImportBatch`].
//!
//! Mapping conventions (no statistics involved):
//! - Each `psm.tsv` row is one rank-1 PSM. Its spectrum is the `Spectrum`
//!   value `<run>.<start scan>.<end scan>.<charge>`: the run is the raw file's
//!   name, the spectrum key is the start scan.
//! - The biosample is the FragPipe experiment, derived from `Spectrum File`
//!   (`interact-` prefix and identification-file extension removed); every
//!   experiment subfolder is a biosample even when no PSM names it. A
//!   biosample's category is its name prefix before the first `-`/`_`
//!   (`control-1` -> `control`), or unassigned when there is none.
//! - The peptide form is `Peptide` plus `Assigned Modifications`
//!   (`5M(15.9949), N-term(42.0106)`, mass deltas).
//! - The engine's score columns are kept under their verbatim header names:
//!   the primary score, `Nextscore`, `Expectation`; a blank cell stays
//!   absent. PeptideProphet's or Percolator's probability and q-value are
//!   not imported (post-processors' verdicts), and since `psm.tsv` is
//!   already filtered, the inspection warns that Tuffold's FDR on it is not
//!   independent; a folder that still holds MSFragger's pepXML files is
//!   imported from those instead (`msfragger.rs`).
//! - Evidence: `Protein ID` with `Protein Start`/`End` and flanking residues,
//!   plus every `Mapped Proteins` entry without positions.
//! - Proteins take their sequence and decoy flag from the bundle FASTA when it
//!   holds them (forward entries win shared identifiers); otherwise the decoy
//!   flag follows the accession convention and the sequence stays unknown.
//!
//! `combined_peptide.tsv` and `combined_protein.tsv` are validated by the
//! inspection but not imported: peptides follow from PSMs, and protein groups
//! are Tuffold's own inference (Phase 2), not FragPipe's.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io;
use std::path::Path;
use std::str::FromStr;

use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeptideForm, Protein, Psm,
    Score, Software, Spectrum, SpectrumKey,
};

use super::inspection::{BundleInspection, inspect_bundle};
use super::{PSM_FILE, biosample, sample_id_from_spectrum_file, split_protein_list};
use crate::fasta::{FastaRecord, build_identifier_index, parse_fasta};
use crate::tsv::{Columns, read_columns};

const PIPELINE: &str = "FragPipe";
/// Score columns imported verbatim, in the order they are reported. The
/// PeptideProphet/Percolator probability and q-value are not: they are
/// verdicts trained on the decoys (see `pepxml.rs`).
const SCORE_COLUMNS: &[&str] = &["score", "nextscore", "expectation"];

#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    /// The bundle failed inspection; the message lists every error.
    #[error("{}", .0.format_errors())]
    InvalidBundle(Box<BundleInspection>),
    #[error("{file}: {source}")]
    Io {
        file: String,
        #[source]
        source: io::Error,
    },
    #[error("{file}: line {line}: column '{column}' holds {value:?}, expected {expected}")]
    InvalidValue {
        file: &'static str,
        line: usize,
        column: &'static str,
        value: String,
        expected: &'static str,
    },
}

/// Inspect and import the bundle at `path`.
pub fn import_bundle(path: impl AsRef<Path>) -> Result<ImportBatch, ParseError> {
    import_inspected(&inspect_bundle(path))
}

/// Import a bundle that has already been inspected.
///
/// Errors in the inspection reject the bundle; warnings never block import.
pub fn import_inspected(inspection: &BundleInspection) -> Result<ImportBatch, ParseError> {
    if !inspection.ok() {
        return Err(ParseError::InvalidBundle(Box::new(inspection.clone())));
    }
    // A clean inspection always resolved psm.tsv; a hand-built one may not.
    let resolved = inspection
        .resolved
        .get(PSM_FILE)
        .ok_or_else(|| ParseError::Io {
            file: PSM_FILE.to_string(),
            source: io::Error::new(
                io::ErrorKind::InvalidInput,
                "the inspection did not resolve psm.tsv; inspect the bundle first",
            ),
        })?;
    let columns = read_columns(&inspection.path.join(PSM_FILE), resolved).map_err(|source| {
        ParseError::Io {
            file: PSM_FILE.to_string(),
            source,
        }
    })?;
    let fasta = match &inspection.selected_fasta {
        Some(path) => Some(parse_fasta(path).map_err(|source| ParseError::Io {
            file: crate::tsv::file_name(path),
            source,
        })?),
        None => None,
    };
    build_batch(inspection, &columns, fasta.as_deref())
}

/// Typed access to `psm.tsv` cells; every bad cell names its line and column.
struct Report<'a> {
    columns: &'a Columns,
}

impl Report<'_> {
    fn invalid(&self, column: &'static str, row: usize, expected: &'static str) -> ParseError {
        ParseError::InvalidValue {
            file: PSM_FILE,
            line: row + 2,
            column,
            value: self.columns.value(column, row).to_string(),
            expected,
        }
    }

    fn text(&self, column: &str, row: usize) -> &str {
        self.columns.value(column, row)
    }

    /// A key cell: blank values would create anonymous or colliding records.
    fn key(&self, column: &'static str, row: usize) -> Result<&str, ParseError> {
        let value = self.text(column, row);
        if value.is_empty() {
            return Err(self.invalid(column, row, "a non-empty value"));
        }
        Ok(value)
    }

    fn integer<T: FromStr>(&self, column: &'static str, row: usize) -> Result<T, ParseError> {
        self.text(column, row)
            .parse()
            .map_err(|_| self.invalid(column, row, "an integer"))
    }

    /// A finite number; NaN and infinities are rejected rather than stored.
    fn number(&self, column: &'static str, row: usize) -> Result<f64, ParseError> {
        self.text(column, row)
            .parse::<f64>()
            .ok()
            .filter(|value| value.is_finite())
            .ok_or_else(|| self.invalid(column, row, "a finite number"))
    }

    /// An optional number: blank or absent is `None`, anything else must parse.
    fn optional_number(&self, column: &'static str, row: usize) -> Result<Option<f64>, ParseError> {
        if self.text(column, row).is_empty() {
            Ok(None)
        } else {
            self.number(column, row).map(Some)
        }
    }

    fn optional_integer<T: FromStr>(
        &self,
        column: &'static str,
        row: usize,
    ) -> Result<Option<T>, ParseError> {
        if self.text(column, row).is_empty() {
            Ok(None)
        } else {
            self.integer(column, row).map(Some)
        }
    }

    fn residue(&self, column: &str, row: usize) -> Option<char> {
        let mut chars = self.text(column, row).chars();
        match (chars.next(), chars.next()) {
            (Some(residue), None) => Some(residue),
            _ => None,
        }
    }
}

fn build_batch(
    inspection: &BundleInspection,
    columns: &Columns,
    fasta: Option<&[FastaRecord]>,
) -> Result<ImportBatch, ParseError> {
    let report = Report { columns };
    let mut biosamples: Vec<BioSample> = inspection
        .experiment_folders
        .iter()
        .map(|name| biosample(name))
        .collect();
    let mut ms_runs: Vec<MsRun> = Vec::new();
    let mut known_runs: HashSet<String> = HashSet::new();
    let mut spectra = Vec::with_capacity(columns.rows());
    let mut psms = Vec::with_capacity(columns.rows());
    let mut descriptions: BTreeMap<String, String> = BTreeMap::new();
    let mut accessions: Vec<String> = Vec::new();
    let mut seen_accessions: HashSet<String> = HashSet::new();

    for row in 0..columns.rows() {
        let spectrum_id = report.key("psm_id", row)?;
        let experiment = sample_id_from_spectrum_file(report.key("sample_id", row)?);
        let (run, key) = split_spectrum(spectrum_id, &experiment);
        if known_runs.insert(run.clone()) {
            if !biosamples.iter().any(|sample| sample.name == experiment) {
                biosamples.push(biosample(&experiment));
            }
            ms_runs.push(MsRun {
                name: run.clone(),
                biosample: experiment.clone(),
                spectra_file: None,
            });
        }
        let charge: i16 = report.integer("charge", row)?;
        let experimental_mz = report.optional_number("observed_mz", row)?;
        spectra.push(Spectrum {
            ms_run: run.clone(),
            key: key.clone(),
            source_id: Some(spectrum_id.to_string()),
            precursor_mz: experimental_mz,
            precursor_charge: Some(charge),
            rt_seconds: report.optional_number("retention", row)?,
        });

        let sequence = report.key("sequence", row)?;
        let modifications = parse_assigned_modifications(
            report.text("assigned_modifications", row),
        )
        .ok_or_else(|| report.invalid("assigned_modifications", row, "items like 5M(15.9949)"))?;
        let mut evidence = Vec::new();
        let mut listed = HashSet::new();
        let primary = split_protein_list(report.text("protein_id", row));
        let mapped = split_protein_list(report.text("mapped_proteins", row));
        for (index, accession) in primary.chain(mapped).enumerate() {
            if !listed.insert(accession) {
                continue;
            }
            let is_primary = index == 0;
            let (start, end, pre, post) = if is_primary {
                (
                    report.optional_integer("protein_start", row)?,
                    report.optional_integer("protein_end", row)?,
                    report.residue("prev_aa", row),
                    report.residue("next_aa", row),
                )
            } else {
                (None, None, None, None)
            };
            evidence.push(Evidence {
                accession: accession.to_string(),
                start,
                end,
                pre,
                post,
            });
            if seen_accessions.insert(accession.to_string()) {
                accessions.push(accession.to_string());
            }
            let description = report.text("protein_description", row);
            if is_primary && !description.is_empty() {
                descriptions
                    .entry(accession.to_string())
                    .or_insert_with(|| description.to_string());
            }
        }

        let mut scores = Vec::new();
        for &logical in SCORE_COLUMNS {
            let Some(name) = inspection.header(PSM_FILE, logical) else {
                continue;
            };
            if let Some(value) = report.optional_number(logical, row)? {
                scores.push(Score {
                    name: name.to_string(),
                    accession: None,
                    value,
                });
            }
        }

        psms.push(Psm {
            ms_run: run,
            spectrum: key,
            rank: 1,
            peptide: PeptideForm::new(sequence, modifications).map_err(|_| {
                report.invalid(
                    "assigned_modifications",
                    row,
                    "modification positions inside the peptide",
                )
            })?,
            charge,
            calculated_mz: report.optional_number("calculated_mz", row)?,
            experimental_mz,
            pass_threshold: None,
            scores,
            evidence,
        });
    }

    Ok(ImportBatch {
        software: Software {
            name: PIPELINE.to_string(),
            version: None,
        },
        biosamples,
        ms_runs,
        spectra,
        proteins: build_proteins(&accessions, &descriptions, fasta),
        psms,
    })
}

/// `<run>.<start>.<end>.<charge>` -> (run, scan key). Values that do not follow
/// the convention keep the whole value as a native key under `fallback_run`.
fn split_spectrum(spectrum: &str, fallback_run: &str) -> (String, SpectrumKey) {
    let mut parts = spectrum.rsplitn(4, '.');
    let (_charge, _end, start, run) = (parts.next(), parts.next(), parts.next(), parts.next());
    match (run, start.and_then(|start| start.parse::<u32>().ok())) {
        (Some(run), Some(scan)) if !run.is_empty() => (run.to_string(), SpectrumKey::Scan(scan)),
        _ => (
            fallback_run.to_string(),
            SpectrumKey::Native(spectrum.to_string()),
        ),
    }
}

/// Parse `5M(15.9949), N-term(42.0106), C-term(-0.9840)`; `None` when any item
/// does not follow the convention (never silently dropped).
fn parse_assigned_modifications(text: &str) -> Option<Vec<Modification>> {
    text.split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(|item| {
            let (site, rest) = item.split_once('(')?;
            let mass_delta = rest
                .strip_suffix(')')?
                .parse::<f64>()
                .ok()
                .filter(|value| value.is_finite())?;
            let site = match site.to_ascii_lowercase().as_str() {
                "n-term" => ModSite::NTerm,
                "c-term" => ModSite::CTerm,
                _ => {
                    let digits = site.trim_end_matches(|ch: char| ch.is_ascii_alphabetic());
                    if digits.is_empty() || digits.len() == site.len() {
                        return None;
                    }
                    ModSite::Residue(digits.parse().ok()?)
                }
            };
            Some(Modification {
                site,
                mass_delta,
                name: None,
                accession: None,
            })
        })
        .collect()
}

fn build_proteins(
    accessions: &[String],
    descriptions: &BTreeMap<String, String>,
    fasta: Option<&[FastaRecord]>,
) -> Vec<Protein> {
    let index: HashMap<String, usize> = fasta.map(build_identifier_index).unwrap_or_default();
    accessions
        .iter()
        .map(|accession| {
            let record = fasta.and_then(|records| index.get(accession).map(|&at| &records[at]));
            Protein {
                accession: accession.clone(),
                description: descriptions.get(accession).cloned(),
                sequence: record.map(|record| record.sequence.clone()),
                is_decoy: record.map_or_else(|| is_decoy_accession(accession), |r| r.is_decoy),
                decoy_explicit: false,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spectrum_values_split_into_run_and_scan() {
        assert_eq!(
            split_spectrum("control-1.1001.1001.2", "x"),
            ("control-1".into(), SpectrumKey::Scan(1001))
        );
        assert_eq!(
            split_spectrum("20200101_HeLa.rep1.00042.00042.3", "x"),
            ("20200101_HeLa.rep1".into(), SpectrumKey::Scan(42))
        );
        assert_eq!(
            split_spectrum("odd-title", "exp"),
            ("exp".into(), SpectrumKey::Native("odd-title".into()))
        );
    }

    #[test]
    fn assigned_modifications_parse_sites_and_deltas() {
        let mods =
            parse_assigned_modifications("5M(15.9949), N-term(42.0106), C-term(-0.9840)").unwrap();
        let sites: Vec<(ModSite, f64)> = mods.iter().map(|m| (m.site, m.mass_delta)).collect();
        assert_eq!(
            sites,
            [
                (ModSite::Residue(5), 15.9949),
                (ModSite::NTerm, 42.0106),
                (ModSite::CTerm, -0.984)
            ]
        );
        assert_eq!(parse_assigned_modifications("").unwrap(), []);
    }

    #[test]
    fn malformed_modification_items_are_rejected() {
        for text in ["M(15.99)", "5M15.99", "5M(abc)", "5M(NaN)", "xx(1.0)"] {
            assert_eq!(parse_assigned_modifications(text), None, "{text}");
        }
    }
}
