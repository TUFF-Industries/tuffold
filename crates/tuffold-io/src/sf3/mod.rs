//! Scaffold `.sf3` projects → one [`ImportBatch`] per search engine.
//!
//! The layout was worked out from sample files alone (see
//! `docs/sf3-format.md`); nothing here comes from Scaffold's code. Reading
//! only: Tuffold never writes `.sf3`.
//!
//! Mapping conventions:
//! - Biosamples, categories, and MS runs come from `Experiment`; a run's
//!   spectra are those its `SpectrumCollation` lists. Scaffold's
//!   "Uncategorized Sample" means no category.
//! - Each engine's candidates for a spectrum are ranked by that engine's
//!   own score (`psms::engine_scores`), never by Scaffold's probability,
//!   which is not imported: Tuffold validates the engines, not Scaffold.
//!   Identical candidates (Scaffold keeps stale copies) count once.
//!   Engines whose score fields are unknown are not imported, and say so.
//! - Spectrum keys: a scan number where the title states one plainly
//!   (`name.850.850.3.pkl`, `... Scan 28842`, `scan=N`), else the title.
//! - Evidence: every protein Scaffold maps the peptide to, plus the one the
//!   engine reported. That one is stored cut at 22 characters (a decoy's
//!   `-DECOY` falls off, leaving its target's name), so a name that long is
//!   completed from Scaffold's map, as a decoy when it ends in part of
//!   `-DECOY`, else by its target's stored sequence (the peptide read
//!   forwards or backwards in it); otherwise the PSM is left out, counted. A stored protein sequence is kept
//!   only when it contains every peptide mapped to it.
//!
//! Checks: every block carries its id, every run holds as many spectra as
//! the project says, every modification sits on a residue its type allows,
//! and a PSM whose mass misses its spectrum's by more than
//! [`MAX_MASS_ERROR`] Da (after isotope errors) is reported, not imported.
//! A file failing any structural check, or with more than
//! [`MAX_MASS_FAILURES`] of its PSMs failing the mass check, is refused.

mod blockbased;
mod java;
mod peaks;
mod project;
mod psms;
mod records;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::ops::Bound;
use std::path::Path;

use tuffold_core::chemistry::PROTON;
use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, MsRun, Protein, Psm, Software, Spectrum, SpectrumKey,
};

pub use peaks::{Sf3Peaks, UnreadablePeaks};
pub(crate) use psms::expectation;
use psms::{Candidate, PEPTIDE_BLOCK, Peptides, Resolved, analysis_of, engine_scores, mass_error};
use records::{Container, Records, read_f32, read_f64, read_i32};

/// Daltons; far beyond any search tolerance, near any decoding error.
pub const MAX_MASS_ERROR: f64 = 1.0;
/// Share of PSMs allowed to fail the mass check before the file is refused.
pub const MAX_MASS_FAILURES: f64 = 0.05;

const SPECTRUM_BLOCK: usize = 48;
const ANALYSIS_BLOCK: usize = 180;
const PROTEIN_BLOCK: usize = 103;
const INT_LIST_BLOCK: usize = 256;
const ACCESSION: std::ops::Range<usize> = 8..31;
const FLANK_BEFORE: usize = 39;
const FLANK_AFTER: usize = 40;

#[derive(Debug, thiserror::Error)]
pub enum Sf3Error {
    #[error("not a Scaffold .sf3 project: {0}")]
    NotSf3(String),
    #[error("cannot read the .sf3 container: {0}")]
    Container(String),
    #[error("the .sf3 file has no {0}, which Tuffold needs")]
    Missing(String),
    #[error("{file}: {error}")]
    Java {
        file: String,
        error: java::JavaError,
    },
    #[error(
        "{file}: {problem}. This file's layout differs from the Scaffold versions Tuffold \
         has been checked on (4.0–5.3.4), so it was not imported rather than guessed at"
    )]
    Layout { file: String, problem: String },
}

pub(crate) fn layout(file: &str, problem: &str) -> Sf3Error {
    Sf3Error::Layout {
        file: file.to_string(),
        problem: problem.to_string(),
    }
}

/// What an `.sf3` import produced: one batch per search engine, notes on
/// everything left out or worth knowing, and the identified spectra's peaks
/// when the project stores them.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Sf3Import {
    pub batches: Vec<ImportBatch>,
    pub notes: Vec<String>,
    pub peaks: Option<Sf3Peaks>,
}

pub fn import_sf3(path: &Path) -> Result<Sf3Import, Sf3Error> {
    let container = Container::open(path)?;
    let runs = project::runs(&container)?;
    let programs = project::programs(&container)?;
    let types = project::modification_types(&container)?;
    let protein_map = project::protein_map(&container)?;
    let stored = project::proteins(&container)?;
    let names = Completions::new(&stored, &protein_map);
    let spectra = container.records("Spectrum", SPECTRUM_BLOCK)?;
    let titles = container.strings("SpectrumName")?;
    let analyses = container.records("SpectrumAnalysis", ANALYSIS_BLOCK)?;
    let protein_records = container.records("ProteinIdentification", PROTEIN_BLOCK)?;
    let records = container.records("PeptideIdentification", PEPTIDE_BLOCK)?;
    let tails = container.strings("PeptideSequence")?;
    let type_lists = container.records("ModificationTypeIds", INT_LIST_BLOCK)?;
    let position_lists = container.records("ModificationLocations", INT_LIST_BLOCK)?;
    let peptides = Peptides {
        tails: &tails,
        type_lists: &type_lists,
        position_lists: &position_lists,
        types: &types,
    };

    let options = peptides.residue_options(&records);
    let spectrum_run = spectrum_runs(&runs, &analyses, &spectra)?;
    let (keys, clashes) = spectrum_keys(&spectrum_run, &spectra, &titles);
    let mut notes = Notes {
        key_clashes: clashes,
        ..Notes::default()
    };

    // (engine, spectrum) → candidates, in record order.
    let mut groups: BTreeMap<(String, i32), Vec<Flanked<Candidate>>> = BTreeMap::new();
    let mut versions: HashMap<String, Option<String>> = HashMap::new();
    for id in records.ids() {
        let block = records
            .get(id)
            .ok_or_else(|| layout("PeptideIdentification", &format!("record {id} is missing")))?;
        let analysis = analysis_of(block);
        let Some(analysis_block) = analyses.get(analysis) else {
            return Err(layout(
                "PeptideIdentification",
                &format!("record {id} names no spectrum analysis"),
            ));
        };
        let spectrum = read_i32(analysis_block, 4).unwrap_or(0);
        let program = analysis_block[8];
        if !spectrum_run.contains_key(&spectrum) {
            notes.without_run.insert(spectrum);
            continue;
        }
        let (name, version) = programs
            .get(&program)
            .cloned()
            .unwrap_or_else(|| (format!("program {program}"), None));
        let Some(engine) = engine_scores(&name) else {
            *notes.unknown_engines.entry(name).or_default() += 1;
            continue;
        };
        if engine.inferred {
            notes.inferred.insert(name.clone());
        }
        let candidate = peptides.candidate(block, &engine).map_err(|problem| {
            layout(
                "PeptideIdentification",
                &format!("record {id}: {}", problem.0),
            )
        })?;
        versions.entry(name.clone()).or_insert(version);
        let flanks = [block[FLANK_BEFORE], block[FLANK_AFTER]];
        groups
            .entry((name, spectrum))
            .or_default()
            .push((candidate, flanks));
    }

    let mut psms: BTreeMap<String, Vec<(i32, Psm)>> = BTreeMap::new();
    let mut checked = 0_usize;
    for ((engine, spectrum), candidates) in groups {
        let block = spectra
            .get(spectrum)
            .ok_or_else(|| layout("Spectrum", &format!("spectrum {spectrum} is missing")))?;
        let run = &runs[spectrum_run[&spectrum]];
        let charge = read_i32(block, 20)
            .and_then(|charge| i16::try_from(charge).ok())
            .filter(|&charge| charge > 0);
        let mh = read_f64(block, 12).filter(|mh| mh.is_finite() && *mh > 0.0);
        let (Some(charge), Some(mh)) = (charge, mh) else {
            notes.without_charge += candidates.len();
            continue;
        };
        let mut resolved = Vec::with_capacity(candidates.len());
        for (candidate, flanks) in candidates {
            let Candidate {
                decoded,
                protein_record,
                primary,
                scores,
            } = candidate;
            match decoded.resolve(mh, &options, MAX_MASS_ERROR) {
                Some((peptide, mass)) => resolved.push((
                    Resolved {
                        peptide,
                        mass,
                        protein_record,
                        primary,
                        scores,
                    },
                    flanks,
                )),
                None => notes.unresolved += 1,
            }
        }
        for (rank, candidate, flanks) in rank(resolved, &mut notes) {
            if let Some(mass) = candidate.mass {
                checked += 1;
                let error = mass_error(mass, mh);
                if error.is_nan() || error.abs() > MAX_MASS_ERROR {
                    notes
                        .mass_failures
                        .push((candidate.peptide.proforma(), error));
                    continue;
                }
            }
            let Some(evidence) =
                evidence(&candidate, flanks, &protein_map, &protein_records, &names)
            else {
                notes.cut_accessions += 1;
                continue;
            };
            let z = f64::from(charge);
            psms.entry(engine.clone()).or_default().push((
                spectrum,
                Psm {
                    ms_run: run.name.clone(),
                    spectrum: keys[&spectrum].clone(),
                    rank,
                    charge,
                    calculated_mz: candidate.mass.map(|mass| (mass + z * PROTON) / z),
                    experimental_mz: Some((mh - PROTON) / z + PROTON),
                    pass_threshold: None,
                    scores: candidate.scores,
                    evidence,
                    peptide: candidate.peptide,
                },
            ));
        }
    }
    if checked > 0 && notes.mass_failures.len() as f64 > MAX_MASS_FAILURES * checked as f64 {
        return Err(layout(
            "PeptideIdentification",
            &format!(
                "{} of {checked} PSMs miss their spectrum's mass by more than {MAX_MASS_ERROR} Da",
                notes.mass_failures.len()
            ),
        ));
    }

    let proteins = proteins(&psms, &stored, &mut notes);
    let peak_refs = peak_refs(&psms, &spectra, &titles);
    let batches = psms
        .into_iter()
        .map(|(engine, psms)| {
            batch(
                &engine,
                versions.remove(&engine).flatten(),
                psms,
                &runs,
                &spectra,
                &titles,
                &keys,
                &proteins,
            )
        })
        .collect();
    let mut lines = notes.lines();
    // Peaks are evidence to look at, not what is validated: a project whose
    // peaks can't be read still imports, and says so.
    let peaks = match peaks::load(&container, peak_refs) {
        Ok(peaks) => peaks,
        Err(error) => {
            lines.push(format!(
                "The spectra's peaks were not read ({error}); attach the runs' spectra files to \
                 see annotated spectra."
            ));
            None
        }
    };
    Ok(Sf3Import {
        batches,
        notes: lines,
        peaks,
    })
}

/// Where each identified spectrum's peaks are (spectra naming no peak list
/// have none).
fn peak_refs(
    psms: &BTreeMap<String, Vec<(i32, Psm)>>,
    spectra: &Records,
    titles: &Records,
) -> Vec<peaks::PeakRef> {
    let mut seen = HashSet::new();
    let mut refs = Vec::new();
    for (spectrum, psm) in psms.values().flatten() {
        if !seen.insert(*spectrum) {
            continue;
        }
        let Some(block) = spectra.get(*spectrum) else {
            continue; // checked when runs were read
        };
        let peak_list = read_i32(block, 4).unwrap_or(0);
        if peak_list <= 0 {
            continue;
        }
        refs.push(peaks::PeakRef {
            run: psm.ms_run.clone(),
            key: psm.spectrum.clone(),
            native_id: titles
                .string(read_i32(block, 8).unwrap_or(0))
                .unwrap_or_else(|| format!("sf3 spectrum {spectrum}")),
            peak_list,
            rt_seconds: read_f32(block, 32)
                .filter(|rt| rt.is_finite() && *rt > 0.0)
                .map(f64::from),
            precursor_mz: psm.experimental_mz,
            precursor_charge: Some(psm.charge),
        });
    }
    refs
}

/// A candidate with the residues before and after it (ASCII bytes).
type Flanked<T> = (T, [u8; 2]);

/// Rank an engine's candidates for one spectrum by its primary score
/// (ties share a rank), dropping repeated copies of the same peptide form.
fn rank(candidates: Vec<(Resolved, [u8; 2])>, notes: &mut Notes) -> Vec<(u32, Resolved, [u8; 2])> {
    let mut seen = HashSet::new();
    let mut unique = Vec::with_capacity(candidates.len());
    for (candidate, flanks) in candidates {
        if !candidate.primary.is_finite() {
            notes.unscored += 1;
        } else if seen.insert(candidate.peptide.proforma()) {
            unique.push((candidate, flanks));
        } else {
            notes.duplicates += 1;
        }
    }
    let scores: Vec<f32> = unique
        .iter()
        .map(|(candidate, _)| candidate.primary)
        .collect();
    unique
        .into_iter()
        .map(|(candidate, flanks)| {
            let better = scores
                .iter()
                .filter(|&&score| score > candidate.primary)
                .count();
            (
                u32::try_from(better + 1).expect("few candidates"),
                candidate,
                flanks,
            )
        })
        .collect()
}

/// The accession a `ProteinIdentification` stores inline holds at most this
/// many characters before its `0x0a`; longer ones are cut there, so a
/// Scaffold 5 decoy `tr|A0QSM5|A0QSM5_MYCS2-DECOY` reads as its target.
const INLINE_ACCESSION_MAX: usize = 22;

/// What the engine's reported protein (stored inline) can be taken as.
#[derive(Debug, PartialEq, Eq)]
enum Reported {
    /// None stored.
    Missing,
    /// Shorter than the inline maximum, so whole.
    Whole(String),
    /// Possibly cut, and the one accession Scaffold maps the peptide to
    /// that starts with it.
    Completed(String),
    /// Possibly cut, and Scaffold's map lists the peptide, but not one
    /// accession starting with it: the map's proteins stand alone.
    MapOnly,
    /// Possibly cut, and the map doesn't list the peptide: whether it names
    /// a target or a decoy needs its sequence (`by_sequence`).
    Unknown(String),
}

fn reported(inline: &[u8], listed: Option<&[String]>) -> Reported {
    let end = inline
        .iter()
        .position(|&byte| byte == 0x0a)
        .unwrap_or(inline.len());
    let accession = String::from_utf8_lossy(&inline[..end]).into_owned();
    if accession.is_empty() {
        return Reported::Missing;
    }
    if end < INLINE_ACCESSION_MAX {
        return Reported::Whole(accession);
    }
    let Some(listed) = listed else {
        return match decoy_completion(&accession) {
            Some(decoy) => Reported::Completed(decoy),
            None => Reported::Unknown(accession),
        };
    };
    let mut completions = listed.iter().filter(|full| full.starts_with(&accession));
    match (completions.next(), completions.next()) {
        (Some(full), None) => Reported::Completed(full.clone()),
        _ => Reported::MapOnly,
    }
}

/// Scaffold's decoys are their target's accession plus `-DECOY`. A cut
/// accession ending in the start of that suffix is such a decoy, whose
/// target is the part before it: completing it keeps the cut from falling
/// on decoys only (`sp|A0QQ48|ECCE3_MYCS2`, 21 characters, is whole; its
/// decoy reads `sp|A0QQ48|ECCE3_MYCS2-`). Where a target's name fills the
/// field, target and decoy read the same and both are left out.
fn decoy_completion(cut: &str) -> Option<String> {
    const SUFFIX: &str = "-DECOY";
    (1..=SUFFIX.len())
        .find(|&kept| cut.len() > kept && cut.ends_with(&SUFFIX[..kept]))
        .map(|kept| format!("{}{SUFFIX}", &cut[..cut.len() - kept]))
}

/// Names a cut accession can complete to, sorted for prefix lookups.
pub(crate) struct Completions<'a> {
    stored: &'a HashMap<String, project::ProteinEntry>,
    /// Stored names with a sequence.
    sequenced: BTreeSet<String>,
    /// Decoy names in Scaffold's protein map.
    decoys: BTreeSet<String>,
}

impl<'a> Completions<'a> {
    fn new(
        stored: &'a HashMap<String, project::ProteinEntry>,
        protein_map: &HashMap<String, Vec<String>>,
    ) -> Self {
        Self {
            stored,
            sequenced: stored
                .iter()
                .filter(|(_, entry)| entry.sequence.is_some())
                .map(|(name, _)| name.clone())
                .collect(),
            decoys: protein_map
                .values()
                .flatten()
                .filter(|name| is_decoy_accession(name))
                .cloned()
                .collect(),
        }
    }

    fn sequence(&self, name: &str) -> Option<&'a str> {
        self.stored
            .get(name)
            .and_then(|entry| entry.sequence.as_deref())
    }

    /// The one name in `names` starting with `cut`.
    fn only<'s>(names: &'s BTreeSet<String>, cut: &str) -> Option<&'s str> {
        let mut found = names
            .range::<str, _>((Bound::Included(cut), Bound::Unbounded))
            .map(String::as_str)
            .take_while(|name| name.starts_with(cut));
        match (found.next(), found.next()) {
            (Some(name), None) => Some(name),
            _ => None,
        }
    }
}

/// A cut accession the map doesn't settle, settled by its target's stored
/// sequence (found under the cut name, the UniProt accession inside it, or
/// the one stored name starting with it). Scaffold's decoys are reversed
/// proteins: a peptide read forwards in the target's sequence (I and L
/// alike) is the target's, read backwards (whole, or but its last residue,
/// for decoys that keep the cleavage site) is its decoy's. Both, neither,
/// or no sequence: `None`. In the Scaffold 5 Sequest sample this settles
/// 121,324 of 122,732 such matches, 121,052 of them as decoys.
fn by_sequence(cut: &str, peptide: &str, names: &Completions<'_>) -> Option<String> {
    let target = if names.sequence(cut).is_some() {
        cut
    } else if let Some(uniprot) =
        tuffold_core::uniprot::uniprot_accession(cut).filter(|id| names.sequence(id).is_some())
    {
        uniprot
    } else {
        Completions::only(&names.sequenced, cut)?
    };
    let leucine = |text: &str| text.replace('I', "L");
    let sequence = leucine(names.sequence(target)?);
    let peptide = leucine(peptide);
    let reversed: String = peptide.chars().rev().collect();
    let reversed_but_last: String = peptide.chars().rev().skip(1).collect();
    let forwards = sequence.contains(&peptide);
    let backwards = sequence.contains(&reversed)
        || (reversed_but_last.len() >= 4 && sequence.contains(&reversed_but_last));
    match (forwards, backwards) {
        (true, false) => Some(target.to_string()),
        // The decoy's name as the map gives it, when it names it.
        (false, true) => Some(
            Completions::only(&names.decoys, cut)
                .map_or_else(|| format!("{cut}-DECOY"), str::to_string),
        ),
        _ => None,
    }
}

/// The PSM's proteins: the one the engine reported, and every protein
/// Scaffold maps its peptide to. `None` when the reported one may be cut
/// and nothing says which protein it is (the PSM is left out).
fn evidence(
    candidate: &Resolved,
    flanks: [u8; 2],
    protein_map: &HashMap<String, Vec<String>>,
    protein_records: &Records,
    names: &Completions<'_>,
) -> Option<Vec<Evidence>> {
    let flank = |byte: u8| (byte.is_ascii_alphabetic() || byte == b'-').then_some(byte as char);
    let peptide = &candidate.peptide.sequence;
    let listed = protein_map.get(peptide);
    let mut evidence: Vec<Evidence> = Vec::new();
    if let Some(block) = protein_records.get(candidate.protein_record) {
        let settled = match reported(&block[ACCESSION], listed.map(Vec::as_slice)) {
            Reported::Unknown(cut) => Reported::Completed(by_sequence(&cut, peptide, names)?),
            other => other,
        };
        match settled {
            Reported::Whole(accession) | Reported::Completed(accession) => {
                evidence.push(Evidence {
                    accession,
                    start: None,
                    end: None,
                    pre: flank(flanks[0]),
                    post: flank(flanks[1]),
                })
            }
            Reported::Missing | Reported::MapOnly | Reported::Unknown(_) => {}
        }
    }
    for accession in listed.into_iter().flatten() {
        if !evidence.iter().any(|known| &known.accession == accession) {
            evidence.push(Evidence {
                accession: accession.clone(),
                start: None,
                end: None,
                pre: None,
                post: None,
            });
        }
    }
    Some(evidence)
}

/// Every protein the PSMs map to; a stored sequence only when it holds all
/// of the protein's peptides (I and L alike).
fn proteins(
    psms: &BTreeMap<String, Vec<(i32, Psm)>>,
    stored: &HashMap<String, project::ProteinEntry>,
    notes: &mut Notes,
) -> HashMap<String, Protein> {
    let mut peptides: HashMap<&str, HashSet<&str>> = HashMap::new();
    for (_, psm) in psms.values().flatten() {
        for evidence in &psm.evidence {
            peptides
                .entry(evidence.accession.as_str())
                .or_default()
                .insert(psm.peptide.sequence.as_str());
        }
    }
    let leucine = |text: &str| text.replace('I', "L");
    peptides
        .into_iter()
        .map(|(accession, sequences)| {
            let entry = stored.get(accession);
            let sequence = entry
                .and_then(|entry| entry.sequence.clone())
                .filter(|sequence| {
                    let protein = leucine(sequence);
                    let holds = sequences
                        .iter()
                        .all(|peptide| protein.contains(&leucine(peptide)));
                    if !holds {
                        notes.dropped_sequences += 1;
                    }
                    holds
                });
            let protein = Protein {
                accession: accession.to_string(),
                description: entry.and_then(|entry| entry.description.clone()),
                sequence,
                is_decoy: is_decoy_accession(accession),
                decoy_explicit: false,
            };
            (accession.to_string(), protein)
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn batch(
    engine: &str,
    version: Option<String>,
    psms: Vec<(i32, Psm)>,
    runs: &[project::Run],
    spectra: &Records,
    titles: &Records,
    keys: &HashMap<i32, SpectrumKey>,
    proteins: &HashMap<String, Protein>,
) -> ImportBatch {
    let mut biosamples: Vec<BioSample> = Vec::new();
    for run in runs {
        if !biosamples.iter().any(|known| known.name == run.biosample) {
            biosamples.push(BioSample {
                name: run.biosample.clone(),
                category: run.category.clone(),
            });
        }
    }
    let mut spectrum_rows = Vec::new();
    let mut seen_spectra = HashSet::new();
    let mut accessions = BTreeMap::new();
    for (spectrum, psm) in &psms {
        for evidence in &psm.evidence {
            accessions.insert(evidence.accession.clone(), ());
        }
        if seen_spectra.insert(*spectrum) {
            let Some(block) = spectra.get(*spectrum) else {
                continue; // checked when runs were read
            };
            let rt = read_f32(block, 32).filter(|rt| rt.is_finite() && *rt > 0.0);
            spectrum_rows.push(Spectrum {
                ms_run: psm.ms_run.clone(),
                key: keys[spectrum].clone(),
                source_id: titles.string(read_i32(block, 8).unwrap_or(0)),
                precursor_mz: psm.experimental_mz,
                precursor_charge: Some(psm.charge),
                rt_seconds: rt.map(f64::from),
            });
        }
    }
    ImportBatch {
        software: Software {
            name: engine.to_string(),
            version,
        },
        biosamples,
        ms_runs: runs
            .iter()
            .map(|run| MsRun {
                name: run.name.clone(),
                biosample: run.biosample.clone(),
                spectra_file: None,
            })
            .collect(),
        spectra: spectrum_rows,
        proteins: accessions
            .keys()
            .filter_map(|accession| proteins.get(accession).cloned())
            .collect(),
        psms: psms.into_iter().map(|(_, psm)| psm).collect(),
    }
}

/// Spectrum id → index of its run, from each run's listed analyses.
fn spectrum_runs(
    runs: &[project::Run],
    analyses: &Records,
    spectra: &Records,
) -> Result<HashMap<i32, usize>, Sf3Error> {
    let mut owner = HashMap::new();
    for (index, run) in runs.iter().enumerate() {
        for &analysis in &run.analyses {
            let block = analyses.get(analysis).ok_or_else(|| {
                layout(
                    "SpectrumCollation",
                    &format!("run {} lists missing analysis {analysis}", run.name),
                )
            })?;
            let spectrum = read_i32(block, 4)
                .filter(|&spectrum| spectra.get(spectrum).is_some())
                .ok_or_else(|| {
                    layout(
                        "SpectrumAnalysis",
                        &format!("analysis {analysis} names no spectrum"),
                    )
                })?;
            if owner
                .insert(spectrum, index)
                .is_some_and(|other| other != index)
            {
                return Err(layout(
                    "SpectrumCollation",
                    &format!("spectrum {spectrum} is in two runs"),
                ));
            }
        }
    }
    Ok(owner)
}

/// Keys for every spectrum in a run; titles that would give two spectra of
/// one run the same key keep a unique form instead.
fn spectrum_keys(
    spectrum_run: &HashMap<i32, usize>,
    spectra: &Records,
    titles: &Records,
) -> (HashMap<i32, SpectrumKey>, usize) {
    let mut keys: HashMap<i32, SpectrumKey> = HashMap::new();
    let mut titles_of: HashMap<i32, Option<String>> = HashMap::new();
    for &spectrum in spectrum_run.keys() {
        let title = spectra
            .get(spectrum)
            .and_then(|block| titles.string(read_i32(block, 8).unwrap_or(0)));
        keys.insert(spectrum, key_from_title(title.as_deref(), spectrum));
        titles_of.insert(spectrum, title);
    }
    let mut uses: HashMap<(usize, &SpectrumKey), usize> = HashMap::new();
    for (spectrum, key) in &keys {
        *uses.entry((spectrum_run[spectrum], key)).or_default() += 1;
    }
    let clashing: Vec<i32> = keys
        .iter()
        .filter(|(spectrum, key)| uses[&(spectrum_run[*spectrum], *key)] > 1)
        .map(|(spectrum, _)| *spectrum)
        .collect();
    let clashes = clashing.len();
    for spectrum in clashing {
        let title = titles_of[&spectrum].as_deref().unwrap_or("spectrum");
        keys.insert(
            spectrum,
            SpectrumKey::Native(format!("{title} [sf3 spectrum {spectrum}]")),
        );
    }
    (keys, clashes)
}

pub(crate) fn key_from_title(title: Option<&str>, spectrum: i32) -> SpectrumKey {
    let Some(title) = title.map(str::trim).filter(|title| !title.is_empty()) else {
        return SpectrumKey::Native(format!("sf3 spectrum {spectrum}"));
    };
    if let SpectrumKey::Scan(scan) = SpectrumKey::from_native_id(title) {
        return SpectrumKey::Scan(scan);
    }
    // DTA-style: name.START.END.CHARGE[.dta|.pkl]
    let mut parts: Vec<&str> = title.split('.').collect();
    if parts
        .last()
        .is_some_and(|ext| ["dta", "pkl", "mgf"].contains(&ext.to_ascii_lowercase().as_str()))
    {
        parts.pop();
    }
    if let [.., start, end, charge] = parts.as_slice()
        && parts.len() >= 4
        && charge.parse::<u8>().is_ok()
        && let (Ok(start), Ok(end)) = (start.parse::<u32>(), end.parse::<u32>())
        && start == end
    {
        return SpectrumKey::Scan(start);
    }
    // "14262: Scan 28842"
    if let Some((_, scan)) = title.rsplit_once("Scan ")
        && let Ok(scan) = scan.parse::<u32>()
    {
        return SpectrumKey::Scan(scan);
    }
    SpectrumKey::Native(title.to_string())
}

#[derive(Default)]
struct Notes {
    without_run: HashSet<i32>,
    without_charge: usize,
    unknown_engines: BTreeMap<String, usize>,
    inferred: std::collections::BTreeSet<String>,
    unscored: usize,
    unresolved: usize,
    duplicates: usize,
    mass_failures: Vec<(String, f64)>,
    dropped_sequences: usize,
    key_clashes: usize,
    cut_accessions: usize,
}

impl Notes {
    fn lines(&self) -> Vec<String> {
        let mut lines = vec![
            "Scaffold's own probabilities were not imported: Tuffold validates the search \
             engines' scores itself."
                .to_string(),
        ];
        if !self.without_run.is_empty() {
            lines.push(format!(
                "{} spectra belong to no MS run in the project and were left out.",
                self.without_run.len()
            ));
        }
        if self.without_charge > 0 {
            lines.push(format!(
                "{} matches on spectra without a precursor charge or mass were left out.",
                self.without_charge
            ));
        }
        for (engine, count) in &self.unknown_engines {
            lines.push(format!(
                "{engine}: {count} matches not imported; Tuffold doesn't yet know where this \
                 engine's scores are stored."
            ));
        }
        for engine in &self.inferred {
            lines.push(format!(
                "{engine}: which stored values are its scores was worked out from the values \
                 alone; their names are marked (inferred) or give the offset."
            ));
        }
        if self.unscored > 0 {
            lines.push(format!(
                "{} matches without their engine's score were left out.",
                self.unscored
            ));
        }
        if self.unresolved > 0 {
            lines.push(format!(
                "{} matches with residues marked modified (lower case) but no single \
                 modification that fits the precursor mass were left out.",
                self.unresolved
            ));
        }
        if self.duplicates > 0 {
            lines.push(format!(
                "{} repeated copies of matches were counted once.",
                self.duplicates
            ));
        }
        if !self.mass_failures.is_empty() {
            let examples: Vec<String> = self
                .mass_failures
                .iter()
                .take(3)
                .map(|(peptide, error)| format!("{peptide} ({error:+.3} Da)"))
                .collect();
            lines.push(format!(
                "{} matches whose peptide mass misses the spectrum's by more than {MAX_MASS_ERROR} Da \
                 were left out, e.g. {}.",
                self.mass_failures.len(),
                examples.join(", ")
            ));
        }
        if self.key_clashes > 0 {
            lines.push(format!(
                "{} spectra share a scan number with another in their run, so they are keyed \
                 by title; peaks from a spectra file won't attach to them.",
                self.key_clashes
            ));
        }
        if self.cut_accessions > 0 {
            lines.push(format!(
                "{} matches were left out: Scaffold stores the engine's protein accession cut \
                 at {INLINE_ACCESSION_MAX} characters, and neither its protein map nor a stored \
                 sequence tells whether it names the target or its decoy.",
                self.cut_accessions
            ));
        }
        if self.dropped_sequences > 0 {
            lines.push(format!(
                "{} stored protein sequences don't contain their peptides and were dropped \
                 (coverage needs a FASTA for those proteins).",
                self.dropped_sequences
            ));
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_reported_accessions_are_completed_from_the_protein_map() {
        let inline = |text: &str| {
            let mut bytes = text.as_bytes().to_vec();
            bytes.push(0x0a);
            bytes.resize(23, 0x0a);
            bytes
        };
        let map = |accessions: &[&str]| -> Vec<String> {
            accessions.iter().map(|a| a.to_string()).collect()
        };
        // 21 characters: whole, whatever the map says.
        assert_eq!(
            reported(&inline("sp|P02768|ALBU_HUMAN1"), None),
            Reported::Whole("sp|P02768|ALBU_HUMAN1".into())
        );
        // 22 characters may be a decoy's name cut short (Scaffold 5).
        let cut = inline("tr|A0QSM5|A0QSM5_MYCS2");
        let decoy = map(&["tr|A0QSM5|A0QSM5_MYCS2-DECOY"]);
        assert_eq!(
            reported(&cut, Some(&decoy)),
            Reported::Completed("tr|A0QSM5|A0QSM5_MYCS2-DECOY".into())
        );
        let target = map(&["tr|A0QSM5|A0QSM5_MYCS2"]);
        assert_eq!(
            reported(&cut, Some(&target)),
            Reported::Completed("tr|A0QSM5|A0QSM5_MYCS2".into())
        );
        // Both, or neither: the map's proteins stand alone.
        let both = map(&["tr|A0QSM5|A0QSM5_MYCS2", "tr|A0QSM5|A0QSM5_MYCS2-DECOY"]);
        assert_eq!(reported(&cut, Some(&both)), Reported::MapOnly);
        assert_eq!(reported(&cut, Some(&map(&["P35527"]))), Reported::MapOnly);
        // Not in the map: the sequence decides (below)...
        assert_eq!(
            reported(&cut, None),
            Reported::Unknown("tr|A0QSM5|A0QSM5_MYCS2".into())
        );
        // ...unless it ends in part of the decoy suffix.
        assert_eq!(
            reported(&inline("sp|A0QQ48|ECCE3_MYCS2-"), None),
            Reported::Completed("sp|A0QQ48|ECCE3_MYCS2-DECOY".into())
        );
        assert_eq!(
            reported(&inline("sp|P0A|NAME_SPECI-DECO"), None),
            Reported::Completed("sp|P0A|NAME_SPECI-DECOY".into())
        );
        assert_eq!(
            reported(&inline("sp|P0A|NAMES_SP-DECOY"), None),
            Reported::Whole("sp|P0A|NAMES_SP-DECOY".into())
        );
        assert_eq!(
            reported(&inline("sp|P0A|NAMES_SPE-DECOY"), None),
            Reported::Completed("sp|P0A|NAMES_SPE-DECOY".into())
        );
        // No terminator at all is as cut as one at 22.
        assert_eq!(
            reported(&[b'A'; 23], None),
            Reported::Unknown("A".repeat(23))
        );
        assert_eq!(reported(&inline(""), None), Reported::Missing);
    }

    #[test]
    fn cut_accessions_outside_the_map_are_settled_by_the_targets_sequence() {
        let stored = HashMap::from([(
            "A0QRV5".to_string(),
            project::ProteinEntry {
                description: None,
                sequence: Some("MSTKLLPEPTIDEKGGAVVSIR".into()),
            },
        )]);
        let map = HashMap::new();
        let names = Completions::new(&stored, &map);
        let cut = "tr|A0QRV5|A0QRV5_MYCS2";
        let settle = |peptide: &str| by_sequence(cut, peptide, &names);
        // Forwards: the target, by the name its sequence is stored under.
        assert_eq!(settle("LLPEPTLDEK").as_deref(), Some("A0QRV5"));
        // Backwards (reversed protein): its decoy.
        assert_eq!(
            settle("EDITPEPLLK").as_deref(),
            Some("tr|A0QRV5|A0QRV5_MYCS2-DECOY")
        );
        // Backwards but the cleavage site kept.
        assert_eq!(
            settle("EDITPEPLLR").as_deref(),
            Some("tr|A0QRV5|A0QRV5_MYCS2-DECOY")
        );
        // Neither, or no sequence: unsettled.
        assert_eq!(settle("WWWWWWK"), None);
        assert_eq!(
            by_sequence("tr|B00000|B00000_MYCS2", "LLPEPTIDEK", &names),
            None
        );
        // A header longer than the field: the target is the one stored name
        // starting with the cut, and the decoy takes the map's name for it.
        let long = "gi|0000000B|emb|AAA64252.1|Ig";
        let stored = HashMap::from([(
            long.to_string(),
            project::ProteinEntry {
                description: None,
                sequence: Some("MSTKLLPEPTIDEKGG".into()),
            },
        )]);
        let map = HashMap::from([("OTHERK".to_string(), vec![format!("{long}-DECOY")])]);
        let names = Completions::new(&stored, &map);
        let cut = "gi|0000000B|emb|AAA642";
        assert_eq!(
            by_sequence(cut, "LLPEPTIDEK", &names).as_deref(),
            Some(long)
        );
        assert_eq!(
            by_sequence(cut, "EDITPEPLLK", &names),
            Some(format!("{long}-DECOY"))
        );
    }

    #[test]
    fn scan_numbers_come_only_from_plain_titles() {
        let key = |title: &str| key_from_title(Some(title), 7);
        assert_eq!(key("Control1.0850.0850.3.pkl"), SpectrumKey::Scan(850));
        assert_eq!(key("run.12.12.2.dta"), SpectrumKey::Scan(12));
        assert_eq!(key("14262: Scan 28842"), SpectrumKey::Scan(28842));
        assert_eq!(key("controllerType=0 scan=55"), SpectrumKey::Scan(55));
        // Summed spectra and ranges have no single scan.
        for title in [
            "run.12.14.2.dta",
            "8144: Sum of 3 scans in range 17366 to 17451 (rt=3862.31)",
            "Discoverer_Tutorial_SILAC-2-47_47",
        ] {
            assert_eq!(key(title), SpectrumKey::Native(title.into()), "{title}");
        }
        assert_eq!(
            key_from_title(None, 7),
            SpectrumKey::Native("sf3 spectrum 7".into())
        );
    }
}
