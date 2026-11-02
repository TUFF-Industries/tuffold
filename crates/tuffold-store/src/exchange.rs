//! Everything a results export (mzIdentML) carries: every validated
//! spectrum, target and decoy, with its peptide, the proteins the peptide
//! maps to, and every protein group, each with its q-value, PEP, and whether
//! it passes the saved thresholds.

use std::collections::HashMap;

use rusqlite::{Connection, params};
use tuffold_core::{ModSite, Modification, PeptideForm, SpectrumKey};
use tuffold_validate::{Level, Thresholds};

use crate::{StoreError, validation};

#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeRun {
    pub id: i64,
    pub name: String,
    pub spectra_file: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeProtein {
    pub id: i64,
    pub accession: String,
    pub description: Option<String>,
    pub sequence: Option<String>,
    pub is_decoy: bool,
}

/// Where a peptide occurs in a protein (`start` 1-based, when known).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ExchangeEvidence {
    pub peptide_id: i64,
    pub protein_id: i64,
    pub start: Option<i64>,
    pub end: Option<i64>,
    pub pre: Option<String>,
    pub post: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeSpectrum {
    pub id: i64,
    pub run_id: i64,
    /// `scan=N`, `index=N`, or a native id.
    pub key: String,
    /// The spectra file's own id for it (mzML's `controllerType=0
    /// controllerNumber=1 scan=N`), when the import recorded one that names
    /// this same spectrum.
    pub native_id: Option<String>,
    pub rt_seconds: Option<f64>,
    pub peptide_id: i64,
    pub charge: Option<i64>,
    pub experimental_mz: Option<f64>,
    pub calculated_mz: Option<f64>,
    pub is_decoy: bool,
    pub score: f64,
    pub q_value: f64,
    pub pep: f64,
    pub passes: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExchangePeptide {
    pub id: i64,
    pub form: PeptideForm,
    /// Peptide-level q-value and PEP, when the peptide was ranked.
    pub q_value: Option<f64>,
    pub pep: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeGroup {
    pub id: i64,
    pub is_decoy: bool,
    pub q_value: Option<f64>,
    pub pep: Option<f64>,
    pub passes: bool,
    /// (protein id, `leader` / `equivalent` / `subset`), leader first.
    pub members: Vec<(i64, String)>,
    /// (peptide id, `unique` / `razor` / `shared`).
    pub peptides: Vec<(i64, String)>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Exchange {
    pub thresholds: Thresholds,
    pub runs: Vec<ExchangeRun>,
    pub proteins: Vec<ExchangeProtein>,
    pub peptides: Vec<ExchangePeptide>,
    pub evidence: Vec<ExchangeEvidence>,
    pub spectra: Vec<ExchangeSpectrum>,
    pub groups: Vec<ExchangeGroup>,
}

pub(crate) fn exchange(conn: &Connection) -> Result<Exchange, StoreError> {
    let thresholds = validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default();
    let measure = validation::measure_column(thresholds.mode);
    let runs = conn
        .prepare("SELECT id, name, spectra_file FROM ms_runs ORDER BY id")?
        .query_map([], |row| {
            Ok(ExchangeRun {
                id: row.get(0)?,
                name: row.get(1)?,
                spectra_file: row.get(2)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    let psm_bound = thresholds.bound(Level::Psm);
    let spectra = conn
        .prepare(&format!(
            "SELECT s.id, s.ms_run_id, s.spectrum_key, s.rt_seconds, sr.peptide_id,
                    (SELECT MIN(p.charge) FROM psms p WHERE p.spectrum_id = s.id
                       AND p.peptide_id = sr.peptide_id AND p.rank = 1),
                    (SELECT MAX(p.experimental_mz) FROM psms p WHERE p.spectrum_id = s.id
                       AND p.peptide_id = sr.peptide_id AND p.rank = 1),
                    (SELECT MAX(p.calculated_mz) FROM psms p WHERE p.spectrum_id = s.id
                       AND p.peptide_id = sr.peptide_id AND p.rank = 1),
                    sr.is_decoy, sr.score, sr.q_value, sr.pep, sr.{measure} <= ?1,
                    s.source_id
             FROM spectrum_results sr JOIN spectra s ON s.id = sr.spectrum_id
             ORDER BY s.ms_run_id, s.id"
        ))?
        .query_map(params![psm_bound], |row| {
            let key: String = row.get(2)?;
            let native_id = row
                .get::<_, Option<String>>(13)?
                .filter(|id| SpectrumKey::from_native_id(id).canonical() == key);
            Ok(ExchangeSpectrum {
                id: row.get(0)?,
                run_id: row.get(1)?,
                key,
                native_id,
                rt_seconds: row.get(3)?,
                peptide_id: row.get(4)?,
                charge: row.get(5)?,
                experimental_mz: row.get(6)?,
                calculated_mz: row.get(7)?,
                is_decoy: row.get(8)?,
                score: row.get(9)?,
                q_value: row.get(10)?,
                pep: row.get(11)?,
                passes: row.get(12)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    let peptides = peptides(conn)?;
    let evidence = conn
        .prepare(
            "SELECT DISTINCT p.peptide_id, pp.protein_id, pp.start_pos, pp.end_pos, pp.pre,
                    pp.post
             FROM psm_proteins pp JOIN psms p ON p.id = pp.psm_id
             WHERE p.peptide_id IN (SELECT peptide_id FROM spectrum_results)
             ORDER BY 1, 2, 3",
        )?
        .query_map([], |row| {
            Ok(ExchangeEvidence {
                peptide_id: row.get(0)?,
                protein_id: row.get(1)?,
                start: row.get(2)?,
                end: row.get(3)?,
                pre: row.get(4)?,
                post: row.get(5)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    let proteins = conn
        .prepare(
            "SELECT id, accession, description, sequence, is_decoy FROM proteins
             WHERE id IN (SELECT protein_id FROM psm_proteins pp JOIN psms p ON p.id = pp.psm_id
                          WHERE p.peptide_id IN (SELECT peptide_id FROM spectrum_results))
                OR id IN (SELECT protein_id FROM protein_group_members)
             ORDER BY id",
        )?
        .query_map([], |row| {
            Ok(ExchangeProtein {
                id: row.get(0)?,
                accession: row.get(1)?,
                description: row.get(2)?,
                sequence: row.get(3)?,
                is_decoy: row.get(4)?,
            })
        })?
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Exchange {
        groups: groups(conn, &thresholds)?,
        thresholds,
        runs,
        proteins,
        peptides,
        evidence,
        spectra,
    })
}

/// Every peptide with a validated spectrum, with its modifications (loaded
/// in one pass rather than peptide by peptide).
fn peptides(conn: &Connection) -> Result<Vec<ExchangePeptide>, StoreError> {
    let mut modifications: HashMap<i64, Vec<Modification>> = HashMap::new();
    let mut statement = conn.prepare(
        "SELECT peptide_id, site_kind, position, mass_delta, name, accession
         FROM modifications
         WHERE peptide_id IN (SELECT peptide_id FROM spectrum_results)
         ORDER BY peptide_id, ordinal",
    )?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        let site = match row.get::<_, String>(1)?.as_str() {
            "n_term" => ModSite::NTerm,
            "c_term" => ModSite::CTerm,
            _ => ModSite::Residue(row.get::<_, Option<u32>>(2)?.unwrap_or(0)),
        };
        modifications
            .entry(row.get(0)?)
            .or_default()
            .push(Modification {
                site,
                mass_delta: row.get(3)?,
                name: row.get(4)?,
                accession: row.get(5)?,
            });
    }
    let mut statement = conn.prepare(
        "SELECT p.id, p.sequence, pr.q_value, pr.pep FROM peptides p
         LEFT JOIN peptide_results pr ON pr.peptide_id = p.id
         WHERE p.id IN (SELECT peptide_id FROM spectrum_results) ORDER BY p.id",
    )?;
    let mut rows = statement.query([])?;
    let mut peptides = Vec::new();
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        let sequence: String = row.get(1)?;
        let form = PeptideForm::new(&sequence, modifications.remove(&id).unwrap_or_default())
            .map_err(|error| StoreError::Damaged(format!("peptide {id}: {error}")))?;
        peptides.push(ExchangePeptide {
            id,
            form,
            q_value: row.get(2)?,
            pep: row.get(3)?,
        });
    }
    Ok(peptides)
}

fn groups(conn: &Connection, thresholds: &Thresholds) -> Result<Vec<ExchangeGroup>, StoreError> {
    let measure = validation::measure_column(thresholds.mode);
    let mut groups: Vec<ExchangeGroup> = conn
        .prepare(&format!(
            "SELECT id, is_decoy, q_value, pep,
                    COALESCE({measure} <= ?1, 0) AND peptides >= ?2
             FROM protein_groups ORDER BY id"
        ))?
        .query_map(
            params![thresholds.bound(Level::Protein), thresholds.min_peptides],
            |row| {
                Ok(ExchangeGroup {
                    id: row.get(0)?,
                    is_decoy: row.get(1)?,
                    q_value: row.get(2)?,
                    pep: row.get(3)?,
                    passes: row.get(4)?,
                    members: Vec::new(),
                    peptides: Vec::new(),
                })
            },
        )?
        .collect::<Result<Vec<_>, _>>()?;
    let index: HashMap<i64, usize> = groups.iter().enumerate().map(|(i, g)| (g.id, i)).collect();
    let mut statement = conn.prepare(
        "SELECT group_id, protein_id, role FROM protein_group_members
         ORDER BY group_id, role != 'leader', protein_id",
    )?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        if let Some(&at) = index.get(&row.get::<_, i64>(0)?) {
            groups[at].members.push((row.get(1)?, row.get(2)?));
        }
    }
    let mut statement = conn.prepare(
        "SELECT group_id, peptide_id, role FROM protein_group_peptides ORDER BY group_id, peptide_id",
    )?;
    let mut rows = statement.query([])?;
    while let Some(row) = rows.next()? {
        if let Some(&at) = index.get(&row.get::<_, i64>(0)?) {
            groups[at].peptides.push((row.get(1)?, row.get(2)?));
        }
    }
    Ok(groups)
}
