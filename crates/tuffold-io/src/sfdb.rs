//! Scaffold `.sfdb` files (Scaffold 5 and Scaffold DDA) → one
//! [`ImportBatch`] per search engine, read-only.
//!
//! An `.sfdb` is a plain SQLite database in an mzIdentML-like layout; see
//! `docs/sfdb-format.md` for the tables and what was checked. The rules are
//! the `.sf3` importer's:
//! - Only the engines' own scores come in (`ScoreIndex`, `PEPTIDE_SCORE_ARRAY`);
//!   Scaffold's probabilities and post-processor verdicts (Percolator) are
//!   left out and named. Ranks are by the engine's main score, over every
//!   PSM of a spectrum, before anything is left out; a spectrum whose best
//!   match can't be read is left out whole, so a lower match never stands
//!   in as rank 1.
//! - A peptide's modifications are the ones its row lists (residue
//!   positions 1-based, 0 for the N-terminus); the N-terminus each PSM
//!   carries is decided by its precursor (`peptides.rs`). Every PSM's mass
//!   must fit its precursor, or it is left out, and the file is refused
//!   when more than 5% miss.
//! - Decoys are Scaffold's `IS_DECOY` flags, else the accession convention.
//! - Runs come from the spectrum titles where they name the raw file
//!   (Scaffold 5), else from the MS sample, which is also the biosample.
//! - Peaks are read too (Scaffold keeps at most 300 in Scaffold 5), attached
//!   as partial peaks that a spectra file replaces.

mod peaks;
mod peptides;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;

use rusqlite::{Connection, OpenFlags};
use tuffold_core::chemistry::PROTON;
use tuffold_core::decoy::is_decoy_accession;
use tuffold_core::{
    BioSample, Evidence, ImportBatch, MsRun, Protein, Psm, Score, Software, Spectrum, SpectrumKey,
};

pub use peaks::SfdbPeaks;
use peptides::{Readings, Tolerance, peptide_neutral};

use crate::sf3::{expectation, key_from_title};
use crate::xml::run_name;

/// More PSMs than this missing their precursor refuses the file.
const MAX_MASS_MISSES: f64 = 0.05;
/// Programs whose output is a verdict on the engines, not engine output:
/// their scores are left out, and they are never an engine.
const VERDICTS: [&str; 3] = ["scaffold", "percolator", "peptideprophet"];

#[derive(Debug, thiserror::Error)]
pub enum SfdbError {
    #[error("not a Scaffold .sfdb file: {0}")]
    NotSfdb(String),
    #[error("cannot read the .sfdb file: {0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error(
        "{0}. The file's layout differs from the Scaffold .sfdb versions Tuffold knows \
         (Scaffold 5.1, Scaffold DDA 6.5), so it was not imported rather than guessed at"
    )]
    Layout(String),
}

/// One batch per search engine, notes on what was left out, and the stored
/// peaks of the identified spectra.
#[derive(Debug, Clone, PartialEq)]
pub struct SfdbImport {
    pub batches: Vec<ImportBatch>,
    pub notes: Vec<String>,
    pub peaks: Option<SfdbPeaks>,
}

const REQUIRED: [&str; 9] = [
    "PeptideSpectrumMatch",
    "Peptide",
    "PeptideEvidence",
    "Protein",
    "Spectrum",
    "MSSample",
    "ScoreIndex",
    "Modification",
    "AnalysisSoftware",
];

pub fn import_sfdb(path: &Path) -> Result<SfdbImport, SfdbError> {
    let conn = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(|error| SfdbError::NotSfdb(error.to_string()))?;
    let tables = table_names(&conn).map_err(|error| SfdbError::NotSfdb(error.to_string()))?;
    if let Some(missing) = REQUIRED.iter().find(|table| !tables.contains(**table)) {
        return Err(SfdbError::NotSfdb(format!("no {missing} table")));
    }
    let version: Option<String> = conn
        .query_row("SELECT version_name FROM Version LIMIT 1", [], |row| {
            row.get(0)
        })
        .ok();
    let mut notes = Notes::default();
    let engines = engines(&conn, &tables)?;
    let scores = scores(&conn, &mut notes)?;
    let modifications = peptides::modifications(&conn)?;
    let samples = samples(&conn)?;
    let spectra = spectra(&conn, &samples, &mut notes)?;
    let peptides = peptides::peptides(&conn, &modifications, &tables)?;
    let evidence = evidence(&conn)?;
    let tolerance = if tables.contains("IdentificationProtocol") {
        Tolerance::read(&conn)
    } else {
        Tolerance::default()
    };

    let raws = rank(read_psms(&conn, &engines, &scores, &mut notes)?, &scores);
    let context = Context {
        spectra: &spectra,
        peptides: &peptides,
        evidence: &evidence,
        scores: &scores,
        tolerance,
    };
    let mut drafts: Vec<Draft> = Vec::new();
    let mut lost: HashSet<(String, i64, String)> = HashSet::new();
    for raw in &raws {
        match context.draft(raw, &mut notes) {
            Ok(draft) => drafts.push(draft),
            Err(reason) => {
                if raw.rank == 1 {
                    lost.insert((raw.engine.clone(), raw.sample, raw.spectrum_id.clone()));
                }
                notes.left_out(reason);
            }
        }
    }
    if notes.mass_misses as f64 > MAX_MASS_MISSES * notes.checked.max(1) as f64 {
        return Err(SfdbError::Layout(format!(
            "{} of {} PSMs miss their precursor mass (e.g. {}), so the modifications were \
             probably misread",
            notes.mass_misses,
            notes.checked,
            notes.mass_examples.as_deref().unwrap_or("")
        )));
    }
    let before = drafts.len();
    drafts.retain(|draft| {
        !lost.contains(&(
            draft.engine.clone(),
            draft.sample,
            draft.spectrum_id.clone(),
        ))
    });
    notes.lost_spectra = lost.len();
    notes.lost_psms = before - drafts.len();
    if drafts.is_empty() {
        return Err(SfdbError::Layout("no PSMs could be read".into()));
    }
    let proteins = proteins(&conn, &drafts, &mut notes)?;
    let peak_refs = peak_refs(&drafts, &spectra);
    let layout = if tables.contains("SpectrumData") && stored_blobs(&conn) {
        peaks::Layout::Blob
    } else {
        peaks::Layout::Text
    };
    let versions: HashMap<String, Option<String>> = engines
        .values()
        .map(|engine| (engine.name.clone(), engine.version.clone()))
        .collect();
    notes.no_decoys = proteins.values().all(|protein| !protein.is_decoy);
    let batches = batches(
        drafts,
        &spectra,
        &samples,
        &proteins,
        &versions,
        version.as_deref(),
    )?;
    Ok(SfdbImport {
        batches,
        notes: notes.lines(&scores),
        peaks: Some(SfdbPeaks::new(path, layout, peak_refs)),
    })
}

/// One `PeptideSpectrumMatch` row with a main score, ranked.
struct Raw {
    id: i64,
    engine: String,
    sample: i64,
    spectrum_id: String,
    peptide_ref: String,
    charge: i64,
    values: Vec<f64>,
    main: f64,
    rank: u32,
}

struct Draft {
    engine: String,
    sample: i64,
    spectrum_id: String,
    run: String,
    key: SpectrumKey,
    psm: Psm,
}

/// What every PSM row is read against.
struct Context<'a> {
    spectra: &'a HashMap<SpectrumId, SpectrumInfo>,
    peptides: &'a HashMap<String, Readings>,
    evidence: &'a HashMap<String, Vec<Evidence>>,
    scores: &'a Scores,
    tolerance: Tolerance,
}

impl Context<'_> {
    /// The PSM a row is, or why it is left out (a mass miss is counted here).
    fn draft(&self, raw: &Raw, notes: &mut Notes) -> Result<Draft, String> {
        let info = self
            .spectra
            .get(&(raw.sample, raw.spectrum_id.clone()))
            .ok_or("its spectrum isn't stored")?;
        let readings = match self.peptides.get(&raw.peptide_ref) {
            Some(Ok(readings)) => readings,
            Some(Err(reason)) => return Err(reason.clone()),
            None => return Err("its peptide isn't stored".into()),
        };
        let charge = i16::try_from(raw.charge)
            .ok()
            .filter(|charge| *charge >= 1)
            .ok_or("it has no charge")?;
        let precursor = info
            .precursor_mz
            .map(|mz| (mz - PROTON) * f64::from(charge));
        let masses: Vec<Option<f64>> = readings
            .iter()
            .map(|reading| peptide_neutral(&reading.form))
            .collect();
        let chosen = match (precursor, masses.first().copied().flatten()) {
            (Some(precursor), Some(_)) => {
                notes.checked += 1;
                let fitting = masses
                    .iter()
                    .position(|mass| mass.is_some_and(|mass| self.tolerance.fits(mass, precursor)));
                let Some(index) = fitting else {
                    notes.mass_misses += 1;
                    notes.mass_examples.get_or_insert_with(|| {
                        format!(
                            "{} {charge}+ at m/z {:.4}",
                            raw.peptide_ref,
                            info.precursor_mz.unwrap_or_default()
                        )
                    });
                    return Err("it misses its precursor mass".into());
                };
                index
            }
            _ => {
                notes.unchecked += 1;
                0
            }
        };
        let reading = &readings[chosen];
        if !reading.named {
            notes.other_terminus += 1;
        }
        if reading.met_lost {
            notes.met_lost += 1;
        }
        let evidence = self
            .evidence
            .get(&raw.peptide_ref)
            .cloned()
            .unwrap_or_default();
        let calculated_mz =
            masses[chosen].map(|mass| (mass + f64::from(charge) * PROTON) / f64::from(charge));
        Ok(Draft {
            engine: raw.engine.clone(),
            sample: raw.sample,
            spectrum_id: raw.spectrum_id.clone(),
            run: info.run.clone(),
            key: info.key.clone(),
            psm: Psm {
                ms_run: info.run.clone(),
                spectrum: info.key.clone(),
                rank: raw.rank,
                peptide: reading.form.clone(),
                charge,
                calculated_mz,
                experimental_mz: info.precursor_mz,
                pass_threshold: None,
                scores: engine_scores(self.scores, &raw.values),
                evidence,
            },
        })
    }
}

fn table_names(conn: &Connection) -> rusqlite::Result<HashSet<String>> {
    let mut statement = conn.prepare("SELECT name FROM sqlite_master WHERE type = 'table'")?;
    let names = statement
        .query_map([], |row| row.get(0))?
        .collect::<rusqlite::Result<_>>()?;
    Ok(names)
}

fn is_verdict(name: &str) -> bool {
    let lowered = name.to_ascii_lowercase();
    VERDICTS.iter().any(|prefix| lowered.starts_with(prefix))
}

#[derive(Debug, Clone)]
struct Engine {
    name: String,
    version: Option<String>,
}

/// Identification protocol → search engine. A protocol credited to
/// Scaffold itself (Scaffold 5 does that) belongs to the file's only other
/// program; with several, the file can't say which ran, and is refused.
fn engines(conn: &Connection, tables: &HashSet<String>) -> Result<HashMap<i64, Engine>, SfdbError> {
    let mut statement =
        conn.prepare("SELECT ID, SOFTWARE_NAME, SOFTWARE_VERSION FROM AnalysisSoftware")?;
    let rows: Vec<(i64, String, Option<String>)> = statement
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))?
        .collect::<rusqlite::Result<_>>()?;
    let versions: HashMap<String, Option<String>> = rows
        .iter()
        .map(|(_, name, version)| (name.clone(), version.clone()))
        .collect();
    let software: HashMap<i64, String> = rows.into_iter().map(|(id, name, _)| (id, name)).collect();
    let others: BTreeSet<&String> = software.values().filter(|name| !is_verdict(name)).collect();
    let sole = (others.len() == 1)
        .then(|| others.iter().next().map(|name| (*name).clone()))
        .flatten();
    let mut statement = conn.prepare("SELECT DISTINCT PROTOCOL_ID FROM PeptideSpectrumMatch")?;
    let used: Vec<i64> = statement
        .query_map([], |row| row.get(0))?
        .collect::<rusqlite::Result<_>>()?;
    let mut protocols: HashMap<i64, i64> = HashMap::new();
    if tables.contains("IdentificationProtocol") {
        let mut statement =
            conn.prepare("SELECT ID, ANALYSIS_SOFTWARE_ID FROM IdentificationProtocol")?;
        protocols = statement
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
            .collect::<rusqlite::Result<_>>()?;
    }
    let mut engines = HashMap::new();
    for protocol in used {
        let named = protocols
            .get(&protocol)
            .and_then(|software_id| software.get(software_id))
            .filter(|name| !is_verdict(name))
            .cloned();
        let Some(engine) = named.or_else(|| sole.clone()) else {
            return Err(SfdbError::Layout(format!(
                "identification protocol {protocol} isn't tied to one search engine (the file \
                 lists {})",
                others
                    .iter()
                    .map(|name| name.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )));
        };
        let version = versions.get(&engine).cloned().flatten();
        engines.insert(
            protocol,
            Engine {
                name: engine,
                version,
            },
        );
    }
    Ok(engines)
}

/// The engine scores in `PEPTIDE_SCORE_ARRAY` order.
#[derive(Default)]
struct Scores {
    /// (index into the array, name, PSI-MS accession)
    kept: Vec<(usize, String, Option<String>)>,
    left_out: Vec<String>,
    /// Index of the score ranks are by.
    main: usize,
    lower_is_better: bool,
    /// Indices of Mascot's ion score and identity threshold, for its
    /// expectation value.
    mascot: Option<(usize, usize)>,
}

fn scores(conn: &Connection, notes: &mut Notes) -> Result<Scores, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT SCORE_INDEX, SCORE_NAME, SCORE_ACCESSION, SCORE_DECREASING FROM ScoreIndex \
         WHERE SCORE_TYPE = 'Peptide Score' ORDER BY SCORE_INDEX",
    )?;
    let rows: Vec<(i64, String, Option<String>, i64)> = statement
        .query_map([], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?))
        })?
        .collect::<rusqlite::Result<_>>()?;
    let mut scores = Scores::default();
    let mut main = None;
    let mut seen = HashSet::new();
    for (index, name, accession, decreasing) in rows {
        let Some(index) = usize::try_from(index).ok().filter(|index| *index < 64) else {
            return Err(SfdbError::Layout(format!("score {name} has index {index}")));
        };
        if !seen.insert(index) {
            return Err(SfdbError::Layout(format!(
                "two peptide scores share index {index}"
            )));
        }
        if is_verdict(&name) {
            scores.left_out.push(name);
            continue;
        }
        let lowered = name.to_ascii_lowercase();
        let ranks_on = !lowered.contains("threshold") && !lowered.contains("rank");
        if ranks_on && main.is_none() {
            main = Some(index);
            scores.lower_is_better = decreasing == 1;
        }
        scores.kept.push((index, name, accession));
    }
    let find = |wanted: &str| {
        scores
            .kept
            .iter()
            .find(|(_, name, _)| name.eq_ignore_ascii_case(wanted))
            .map(|(index, _, _)| *index)
    };
    scores.mascot = find("Mascot:score").zip(find("Mascot:identity threshold"));
    let Some(main) = main else {
        return Err(SfdbError::Layout(format!(
            "no engine score to rank by among {}",
            scores
                .kept
                .iter()
                .map(|(_, name, _)| name.as_str())
                .chain(scores.left_out.iter().map(String::as_str))
                .collect::<Vec<_>>()
                .join(", ")
        )));
    };
    scores.main = main;
    notes.scores_read = scores
        .kept
        .iter()
        .map(|(_, name, _)| name.clone())
        .collect();
    Ok(scores)
}

fn parse_scores(array: &str) -> Vec<f64> {
    array
        .split(',')
        .map(|value| value.trim().parse::<f64>().unwrap_or(f64::NAN))
        .collect()
}

fn value_at(values: &[f64], index: usize) -> Option<f64> {
    values.get(index).copied().filter(|value| value.is_finite())
}

fn engine_scores(scores: &Scores, values: &[f64]) -> Vec<Score> {
    let mut out: Vec<Score> = scores
        .kept
        .iter()
        .filter_map(|(index, name, accession)| {
            value_at(values, *index).map(|value| Score {
                name: name.clone(),
                accession: accession.clone(),
                value,
            })
        })
        .collect();
    if let Some((score, identity)) = scores.mascot
        && let (Some(score), Some(identity)) = (value_at(values, score), value_at(values, identity))
    {
        let value = expectation(score, identity);
        if value.is_finite() {
            out.push(Score {
                name: "Mascot:expectation value".into(),
                accession: None,
                value,
            });
        }
    }
    out
}

/// Every PSM row. A spectrum with a match lacking the main score can't be
/// ranked, so all its matches are left out, counted.
fn read_psms(
    conn: &Connection,
    engines: &HashMap<i64, Engine>,
    scores: &Scores,
    notes: &mut Notes,
) -> Result<Vec<Raw>, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT ID, MSSAMPLE_ID, SPECTRUM_ID, PEPTIDE_REF, ASSIGNED_CHARGE, \
         PEPTIDE_SCORE_ARRAY, PROTOCOL_ID FROM PeptideSpectrumMatch",
    )?;
    let mut rows = statement.query([])?;
    let mut raws = Vec::new();
    let mut unranked: HashSet<(String, i64, String)> = HashSet::new();
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        let protocol: i64 = row.get(6)?;
        let Some(engine) = engines.get(&protocol) else {
            return Err(SfdbError::Layout(format!(
                "PSM {id} names identification protocol {protocol}, which the file doesn't \
                 tie to a search engine"
            )));
        };
        let values = parse_scores(&row.get::<_, String>(5)?);
        let sample: i64 = row.get(1)?;
        let spectrum_id: String = row.get(2)?;
        let main = value_at(&values, scores.main);
        if main.is_none() {
            unranked.insert((engine.name.clone(), sample, spectrum_id.clone()));
        }
        raws.push(Raw {
            id,
            engine: engine.name.clone(),
            sample,
            spectrum_id,
            peptide_ref: row.get(3)?,
            charge: row.get(4)?,
            values,
            main: main.unwrap_or(f64::NAN),
            rank: 0,
        });
    }
    let before = raws.len();
    raws.retain(|raw| {
        !unranked.contains(&(raw.engine.clone(), raw.sample, raw.spectrum_id.clone()))
    });
    if before > raws.len() {
        notes.left_out_n(
            "a match of their spectrum lacks the engine's main score, so ranks are unknown",
            before - raws.len(),
        );
    }
    Ok(raws)
}

/// Rank each spectrum's PSMs (per engine) by the main score; ties keep the
/// file's order (by PSM id).
fn rank(mut raws: Vec<Raw>, scores: &Scores) -> Vec<Raw> {
    raws.sort_by(|a, b| {
        (&a.engine, a.sample, &a.spectrum_id)
            .cmp(&(&b.engine, b.sample, &b.spectrum_id))
            .then_with(|| {
                let order = a.main.total_cmp(&b.main);
                if scores.lower_is_better {
                    order
                } else {
                    order.reverse()
                }
            })
            .then(a.id.cmp(&b.id))
    });
    for index in 0..raws.len() {
        let rank = match index.checked_sub(1).map(|previous| &raws[previous]) {
            Some(previous)
                if (&previous.engine, previous.sample, &previous.spectrum_id)
                    == (
                        &raws[index].engine,
                        raws[index].sample,
                        &raws[index].spectrum_id,
                    ) =>
            {
                previous.rank + 1
            }
            _ => 1,
        };
        raws[index].rank = rank;
    }
    raws
}

/// MS sample id → name (Scaffold's biosample).
fn samples(conn: &Connection) -> Result<HashMap<i64, String>, SfdbError> {
    let mut statement = conn.prepare("SELECT ID, SAMPLE_NAME FROM MSSample")?;
    let samples: HashMap<i64, String> = statement
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .collect::<rusqlite::Result<_>>()?;
    if samples.is_empty() {
        return Err(SfdbError::Layout("no MS samples".into()));
    }
    Ok(samples)
}

/// What one stored spectrum is.
#[derive(Debug, Clone)]
struct SpectrumInfo {
    run: String,
    /// The raw file the title names, verbatim.
    file: Option<String>,
    key: SpectrumKey,
    title: Option<String>,
    precursor_mz: Option<f64>,
    precursor_charge: Option<i16>,
    rt_seconds: Option<f64>,
}

type SpectrumId = (i64, String);

fn spectra(
    conn: &Connection,
    samples: &HashMap<i64, String>,
    notes: &mut Notes,
) -> Result<HashMap<SpectrumId, SpectrumInfo>, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT MSSAMPLE_ID, SPECTRUM_ID, SPECTRUM_TITLE, ACQUISITION_MASS, ACQUISITION_CHARGE, \
         SCAN_START_TIME_IN_SEC FROM Spectrum",
    )?;
    let mut rows = statement.query([])?;
    let mut out: HashMap<SpectrumId, SpectrumInfo> = HashMap::new();
    while let Some(row) = rows.next()? {
        let sample: i64 = row.get(0)?;
        let spectrum_id: String = row.get(1)?;
        let title: Option<String> = row.get(2)?;
        let Some(sample_name) = samples.get(&sample) else {
            return Err(SfdbError::Layout(format!(
                "spectrum {spectrum_id} is in MS sample {sample}, which isn't listed"
            )));
        };
        let (run, file, key) = locate(title.as_deref(), sample_name, &spectrum_id);
        out.insert(
            (sample, spectrum_id),
            SpectrumInfo {
                run,
                file,
                key,
                title,
                precursor_mz: row
                    .get::<_, Option<f64>>(3)?
                    .filter(|mz| mz.is_finite() && *mz > 0.0),
                precursor_charge: row
                    .get::<_, Option<i64>>(4)?
                    .and_then(|charge| i16::try_from(charge).ok())
                    .filter(|charge| *charge > 0),
                rt_seconds: row.get::<_, Option<f64>>(5)?.filter(|rt| rt.is_finite()),
            },
        );
    }
    // Two spectra of one run with one key keep a unique form instead.
    let mut uses: HashMap<(String, SpectrumKey), usize> = HashMap::new();
    for info in out.values() {
        *uses
            .entry((info.run.clone(), info.key.clone()))
            .or_default() += 1;
    }
    for ((sample, spectrum_id), info) in &mut out {
        if uses[&(info.run.clone(), info.key.clone())] > 1 {
            notes.key_clashes += 1;
            info.key = SpectrumKey::Native(format!(
                "{} [sfdb spectrum {sample}:{spectrum_id}]",
                info.title.as_deref().unwrap_or("spectrum")
            ));
        }
    }
    Ok(out)
}

/// The run, raw file, and key of a spectrum from its title: Scaffold 5's
/// `File: "…\X.raw"; SpectrumID: "1445"; scans: "5791"` names both;
/// Scaffold DDA's `1-7_7` (id, first and last scan) names the scan and the
/// MS sample is the run; anything else as `.sf3` titles are read.
fn locate(
    title: Option<&str>,
    sample: &str,
    spectrum_id: &str,
) -> (String, Option<String>, SpectrumKey) {
    let title = title.map(str::trim).filter(|title| !title.is_empty());
    if let Some(title) = title
        && let Some(file) = quoted(title, "File:")
    {
        let key = quoted(title, "scans:")
            .and_then(|scans| scans.parse::<u32>().ok())
            .map_or_else(|| SpectrumKey::Native(title.to_string()), SpectrumKey::Scan);
        return (run_name(&file), Some(file), key);
    }
    let key = match title {
        Some(title) => dda_scan(title).map_or_else(
            || match key_from_title(Some(title), 0) {
                SpectrumKey::Scan(scan) => SpectrumKey::Scan(scan),
                _ => SpectrumKey::Native(title.to_string()),
            },
            SpectrumKey::Scan,
        ),
        None => SpectrumKey::Native(format!("sfdb spectrum {spectrum_id}")),
    };
    (sample.to_string(), None, key)
}

/// The value after `label` in `label "value"`.
fn quoted(title: &str, label: &str) -> Option<String> {
    let rest = &title[title.find(label)? + label.len()..];
    let rest = rest.trim_start().strip_prefix('"')?;
    Some(rest[..rest.find('"')?].to_string())
}

/// `ID-FIRST_LAST` with one scan.
fn dda_scan(title: &str) -> Option<u32> {
    let (id, scans) = title.split_once('-')?;
    let (first, last) = scans.split_once('_')?;
    id.parse::<u32>().ok()?;
    let (first, last) = (first.parse::<u32>().ok()?, last.parse::<u32>().ok()?);
    (first == last).then_some(first)
}

/// Peptide ref → where it sits on each protein, once per protein (a second
/// place on the same protein, a repeat, is not kept). A start below 1
/// (Scaffold DDA 7.1 is reported to write −1) is rebuilt from the stop and
/// the peptide's length.
fn evidence(conn: &Connection) -> Result<HashMap<String, Vec<Evidence>>, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT e.PEPTIDE_REF, e.PROTEIN_ACCESSION, e.PROTEIN_START_INDEX, \
         e.PROTEIN_STOP_INDEX, e.PRE_AA, e.POST_AA, length(p.PEPTIDE_SEQUENCE) \
         FROM PeptideEvidence e LEFT JOIN Peptide p ON p.PEPTIDE_REF = e.PEPTIDE_REF \
         ORDER BY e.PEPTIDE_REF, e.PROTEIN_ACCESSION, e.PROTEIN_START_INDEX",
    )?;
    let mut rows = statement.query([])?;
    let mut out: HashMap<String, Vec<Evidence>> = HashMap::new();
    while let Some(row) = rows.next()? {
        let reference: String = row.get(0)?;
        let accession: String = row.get(1)?;
        let position = |value: Option<i64>| {
            value
                .and_then(|value| u32::try_from(value).ok())
                .filter(|value| *value > 0)
        };
        let end = position(row.get(3)?);
        let length = position(row.get(6)?);
        let start = position(row.get(2)?).or_else(|| {
            end.zip(length)
                .filter(|(end, length)| end >= length)
                .map(|(end, length)| end - length + 1)
        });
        let residue = |value: Option<String>| {
            value
                .and_then(|value| value.chars().next())
                .filter(char::is_ascii_alphabetic)
        };
        let pre = residue(row.get(4)?);
        let post = residue(row.get(5)?);
        let list = out.entry(reference).or_default();
        if list.iter().any(|e| e.accession == accession) {
            continue;
        }
        list.push(Evidence {
            accession,
            start,
            end,
            pre,
            post,
        });
    }
    Ok(out)
}

/// The proteins PSMs name. A protein Scaffold flags `IS_DECOY` 1 is a decoy
/// (explicitly); one flagged 0 follows the accession convention, since the
/// flag is only Scaffold's own decoy pattern, not the engine's statement. A
/// stored sequence is kept only if it holds every peptide placed on it
/// (I = L).
fn proteins(
    conn: &Connection,
    drafts: &[Draft],
    notes: &mut Notes,
) -> Result<HashMap<String, Protein>, SfdbError> {
    let mut peptides_of: HashMap<&str, HashSet<&str>> = HashMap::new();
    for draft in drafts {
        for evidence in &draft.psm.evidence {
            peptides_of
                .entry(evidence.accession.as_str())
                .or_default()
                .insert(draft.psm.peptide.sequence.as_str());
        }
    }
    let mut statement = conn.prepare(
        "SELECT ACCESSION_NUMBER, ORIGINAL_PROTEIN_NAME, PROTEIN_SEQUENCE, IS_DECOY FROM Protein",
    )?;
    let mut rows = statement.query([])?;
    let mut out = HashMap::new();
    while let Some(row) = rows.next()? {
        let accession: String = row.get(0)?;
        let Some(peptides) = peptides_of.get(accession.as_str()) else {
            continue;
        };
        let description: Option<String> = row.get(1)?;
        let sequence: Option<String> = row.get(2)?;
        let flagged: Option<i64> = row.get(3)?;
        let sequence = sequence
            .map(|sequence| sequence.trim().to_ascii_uppercase())
            .filter(|sequence| !sequence.is_empty());
        let sequence = sequence.filter(|sequence| {
            let leucine = sequence.replace('I', "L");
            let holds = peptides
                .iter()
                .all(|peptide| leucine.contains(&peptide.replace('I', "L")));
            if !holds {
                notes.dropped_sequences += 1;
            }
            holds
        });
        let explicit = flagged == Some(1);
        out.insert(
            accession.clone(),
            Protein {
                description: description.filter(|text| !text.is_empty() && *text != accession),
                is_decoy: explicit || is_decoy_accession(&accession),
                decoy_explicit: explicit,
                accession,
                sequence,
            },
        );
    }
    notes.unlisted_proteins = peptides_of
        .keys()
        .filter(|accession| !out.contains_key(**accession))
        .count();
    Ok(out)
}

fn stored_blobs(conn: &Connection) -> bool {
    conn.query_row("SELECT 1 FROM SpectrumData LIMIT 1", [], |row| {
        row.get::<_, i64>(0)
    })
    .is_ok()
}

/// Where each identified spectrum's peaks are, once per spectrum.
fn peak_refs(drafts: &[Draft], spectra: &HashMap<SpectrumId, SpectrumInfo>) -> Vec<peaks::PeakRef> {
    let mut seen = HashSet::new();
    let mut ids: Vec<&SpectrumId> = spectra.keys().collect();
    ids.sort();
    let identified: HashSet<(&str, &SpectrumKey)> = drafts
        .iter()
        .map(|draft| (draft.run.as_str(), &draft.key))
        .collect();
    ids.into_iter()
        .filter_map(|id| {
            let info = &spectra[id];
            let wanted = identified.contains(&(info.run.as_str(), &info.key))
                && seen.insert((info.run.clone(), info.key.clone()));
            wanted.then(|| peaks::PeakRef {
                run: info.run.clone(),
                key: info.key.clone(),
                sample: id.0,
                spectrum_id: id.1.clone(),
                native_id: info
                    .title
                    .clone()
                    .unwrap_or_else(|| format!("sfdb spectrum {}", id.1)),
                rt_seconds: info.rt_seconds,
                precursor_mz: info.precursor_mz,
                precursor_charge: info.precursor_charge,
            })
        })
        .collect()
}

fn batches(
    drafts: Vec<Draft>,
    spectra: &HashMap<SpectrumId, SpectrumInfo>,
    samples: &HashMap<i64, String>,
    proteins: &HashMap<String, Protein>,
    engine_versions: &HashMap<String, Option<String>>,
    scaffold: Option<&str>,
) -> Result<Vec<ImportBatch>, SfdbError> {
    // Each run belongs to one MS sample; its file is the first title's, in
    // spectrum order.
    let mut ids: Vec<&SpectrumId> = spectra.keys().collect();
    ids.sort();
    let mut run_sample: BTreeMap<String, (String, Option<String>)> = BTreeMap::new();
    for id in ids {
        let info = &spectra[id];
        let name = &samples[&id.0];
        match run_sample.get(&info.run) {
            Some((other, _)) if other != name => {
                return Err(SfdbError::Layout(format!(
                    "run {} is in two MS samples ({other}, {name})",
                    info.run
                )));
            }
            Some(_) => {}
            None => {
                run_sample.insert(info.run.clone(), (name.clone(), info.file.clone()));
            }
        }
    }
    let infos: HashMap<(&str, &SpectrumKey), &SpectrumInfo> = spectra
        .values()
        .map(|info| ((info.run.as_str(), &info.key), info))
        .collect();
    let scaffold = scaffold.map(|version| version.trim_start_matches("Scaffold_").to_string());
    let mut by_engine: BTreeMap<String, Vec<Draft>> = BTreeMap::new();
    for draft in drafts {
        by_engine
            .entry(draft.engine.clone())
            .or_default()
            .push(draft);
    }
    let mut batches = Vec::new();
    for (engine, drafts) in by_engine {
        let mut runs: BTreeSet<String> = BTreeSet::new();
        let mut seen: HashSet<(String, SpectrumKey)> = HashSet::new();
        let mut batch_spectra = Vec::new();
        let mut used: BTreeMap<String, Protein> = BTreeMap::new();
        let mut psms = Vec::with_capacity(drafts.len());
        for draft in drafts {
            runs.insert(draft.run.clone());
            if seen.insert((draft.run.clone(), draft.key.clone())) {
                let info = infos[&(draft.run.as_str(), &draft.key)];
                batch_spectra.push(Spectrum {
                    ms_run: draft.run.clone(),
                    key: draft.key.clone(),
                    source_id: info.title.clone(),
                    precursor_mz: info.precursor_mz,
                    precursor_charge: info.precursor_charge,
                    rt_seconds: info.rt_seconds,
                });
            }
            for evidence in &draft.psm.evidence {
                used.entry(evidence.accession.clone()).or_insert_with(|| {
                    proteins
                        .get(&evidence.accession)
                        .cloned()
                        .unwrap_or_else(|| Protein {
                            accession: evidence.accession.clone(),
                            description: None,
                            sequence: None,
                            is_decoy: is_decoy_accession(&evidence.accession),
                            decoy_explicit: false,
                        })
                });
            }
            psms.push(draft.psm);
        }
        let biosamples: BTreeSet<&String> = runs.iter().map(|run| &run_sample[run].0).collect();
        let version = engine_versions.get(&engine).cloned().flatten();
        let version = match (version, &scaffold) {
            (Some(version), Some(scaffold)) => Some(format!("{version} (via Scaffold {scaffold})")),
            (None, Some(scaffold)) => Some(format!("via Scaffold {scaffold}")),
            (version, None) => version,
        };
        batches.push(ImportBatch {
            software: Software {
                name: engine,
                version,
            },
            biosamples: biosamples
                .into_iter()
                .map(|name| BioSample {
                    name: name.clone(),
                    category: None,
                })
                .collect(),
            ms_runs: runs
                .iter()
                .map(|run| MsRun {
                    name: run.clone(),
                    biosample: run_sample[run].0.clone(),
                    spectra_file: run_sample[run].1.clone(),
                })
                .collect(),
            spectra: batch_spectra,
            proteins: used.into_values().collect(),
            psms,
        });
    }
    Ok(batches)
}

#[derive(Default)]
struct Notes {
    left_out: BTreeMap<String, usize>,
    checked: usize,
    unchecked: usize,
    mass_misses: usize,
    mass_examples: Option<String>,
    other_terminus: usize,
    met_lost: usize,
    lost_spectra: usize,
    lost_psms: usize,
    key_clashes: usize,
    dropped_sequences: usize,
    unlisted_proteins: usize,
    scores_read: Vec<String>,
    no_decoys: bool,
}

impl Notes {
    fn left_out(&mut self, reason: String) {
        self.left_out_n(&reason, 1);
    }

    fn left_out_n(&mut self, reason: &str, count: usize) {
        *self.left_out.entry(reason.to_string()).or_default() += count;
    }

    fn lines(&self, scores: &Scores) -> Vec<String> {
        let mut lines = vec![format!(
            "Read the engine's own scores ({}); Scaffold's probabilities and other \
             verdicts{} were not imported: Tuffold validates the search engines' scores itself.",
            self.scores_read.join(", "),
            if scores.left_out.is_empty() {
                String::new()
            } else {
                format!(" ({})", scores.left_out.join(", "))
            }
        )];
        if self.no_decoys {
            lines.push(
                "The file holds no decoy proteins (no decoy search, or Scaffold kept none), so \
                 Tuffold can't estimate an FDR for it."
                    .into(),
            );
        }
        for (reason, count) in &self.left_out {
            lines.push(format!("{count} PSMs were left out: {reason}."));
        }
        if self.mass_misses > 0 {
            lines.push(format!(
                "{} of {} PSMs checked miss their precursor mass under every reading of their \
                 peptide (e.g. {}).",
                self.mass_misses,
                self.checked,
                self.mass_examples.as_deref().unwrap_or("")
            ));
        }
        if self.lost_spectra > 0 {
            lines.push(format!(
                "{} spectra whose best match was left out were left out whole ({} more PSMs), \
                 so a lower match never counts as rank 1.",
                self.lost_spectra, self.lost_psms
            ));
        }
        if self.unchecked > 0 {
            lines.push(format!(
                "{} PSMs couldn't be checked against their precursor mass (no precursor m/z, \
                 or a residue without a mass).",
                self.unchecked
            ));
        }
        if self.other_terminus > 0 {
            lines.push(format!(
                "{} PSMs fit their precursor only with N-terminal modifications other than \
                 their Scaffold peptide name shows ({} of them with the initial methionine \
                 lost); Scaffold files such PSMs under one peptide.",
                self.other_terminus, self.met_lost
            ));
        }
        if self.key_clashes > 0 {
            lines.push(format!(
                "{} spectra share a scan number in their run and keep their full title as key; \
                 peaks from a spectra file won't attach to them.",
                self.key_clashes
            ));
        }
        if self.dropped_sequences > 0 {
            lines.push(format!(
                "{} stored protein sequences don't contain their peptides and were not kept.",
                self.dropped_sequences
            ));
        }
        if self.unlisted_proteins > 0 {
            lines.push(format!(
                "{} proteins PSMs name aren't in the file's protein table; they are kept by \
                 accession alone (decoys by the accession convention).",
                self.unlisted_proteins
            ));
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titles_name_run_and_scan() {
        let (run, file, key) = locate(
            Some(r#"File: "Z:\CORE\RAW\P706_WT_7T.raw"; SpectrumID: "1445"; scans: "5791""#),
            "mascot_daemon_merge (F015729)",
            "index=2",
        );
        assert_eq!(run, "P706_WT_7T");
        assert_eq!(file.as_deref(), Some(r"Z:\CORE\RAW\P706_WT_7T.raw"));
        assert_eq!(key, SpectrumKey::Scan(5791));
        assert_eq!(
            locate(Some("1-7_7"), "Elute-1-1", "1"),
            ("Elute-1-1".to_string(), None, SpectrumKey::Scan(7))
        );
        // Summed scans have no single scan.
        assert_eq!(
            locate(Some("4-110_112"), "Elute-1-1", "4").2,
            SpectrumKey::Native("4-110_112".into())
        );
        assert_eq!(
            locate(None, "s", "9").2,
            SpectrumKey::Native("sfdb spectrum 9".into())
        );
    }

    fn raw(id: i64, spectrum: &str, main: f64) -> Raw {
        Raw {
            id,
            engine: "E".into(),
            sample: 1,
            spectrum_id: spectrum.into(),
            peptide_ref: String::new(),
            charge: 2,
            values: Vec::new(),
            main,
            rank: 0,
        }
    }

    #[test]
    fn ranks_follow_the_main_score_per_spectrum() {
        let higher = Scores::default();
        let ranked = rank(
            vec![
                raw(1, "a", 10.0),
                raw(2, "a", 30.0),
                raw(3, "b", 5.0),
                raw(4, "a", 30.0),
            ],
            &higher,
        );
        let ranks: Vec<(i64, u32)> = ranked.iter().map(|r| (r.id, r.rank)).collect();
        // Ties keep the file's order.
        assert_eq!(ranks, [(2, 1), (4, 2), (1, 3), (3, 1)]);
        let lower = Scores {
            lower_is_better: true,
            ..Scores::default()
        };
        let ranked = rank(vec![raw(1, "a", 0.9), raw(2, "a", 0.01)], &lower);
        assert_eq!(ranked[0].id, 2);
    }
}
