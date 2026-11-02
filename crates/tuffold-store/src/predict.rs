//! What deep-learned rescoring reads: the rank-1 PSMs' peptides and spectra.

use std::collections::HashMap;

use rusqlite::{Connection, OptionalExtension};
use tuffold_core::{ModSite, Modification, PeptideForm};

use crate::StoreError;
use crate::peaks::{self, Peaks};

/// A rank-1 PSM with its spectrum's peaks.
#[derive(Debug, Clone, PartialEq)]
pub struct PsmPeaks {
    pub psm: i64,
    pub peptide: i64,
    pub charge: i32,
    pub peaks: Peaks,
}

/// Every peptide a rank-1 PSM names, with its modifications.
pub(crate) fn rank1_peptides(conn: &Connection) -> Result<HashMap<i64, PeptideForm>, StoreError> {
    let mut modifications: HashMap<i64, Vec<Modification>> = HashMap::new();
    let mut statement = conn.prepare(
        "SELECT peptide_id, site_kind, position, mass_delta, name, accession
         FROM modifications
         WHERE peptide_id IN (SELECT peptide_id FROM psms WHERE rank = 1)
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
        "SELECT id, sequence FROM peptides
         WHERE id IN (SELECT peptide_id FROM psms WHERE rank = 1)",
    )?;
    let mut rows = statement.query([])?;
    let mut forms = HashMap::new();
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        let sequence: String = row.get(1)?;
        let form = PeptideForm::new(&sequence, modifications.remove(&id).unwrap_or_default())
            .map_err(|error| StoreError::Damaged(format!("peptide {id}: {error}")))?;
        forms.insert(id, form);
    }
    Ok(forms)
}

/// Whether any rank-1 PSM's spectrum has peaks (what the deep-learning
/// features are computed from).
pub(crate) fn has_rank1_peaks(conn: &Connection) -> Result<bool, StoreError> {
    Ok(conn.query_row(
        "SELECT EXISTS (SELECT 1 FROM psms p JOIN spectra s ON s.id = p.spectrum_id
         WHERE p.rank = 1 AND s.peaks IS NOT NULL)",
        [],
        |row| row.get(0),
    )?)
}

/// The rank-1 PSMs whose spectra have peaks, with their spectra.
pub(crate) fn rank1_with_peaks(conn: &Connection) -> Result<Vec<(i64, i64)>, StoreError> {
    let mut statement = conn.prepare(
        "SELECT p.id, p.spectrum_id FROM psms p JOIN spectra s ON s.id = p.spectrum_id
         WHERE p.rank = 1 AND s.peaks IS NOT NULL ORDER BY p.id",
    )?;
    let ids = statement
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .collect::<Result<_, _>>()?;
    Ok(ids)
}

/// Rank-1 PSMs whose spectrum has peaks, `chunk` at a time (a project's
/// peaks may not fit in memory at once).
pub(crate) fn for_each_rank1_spectrum<E: From<StoreError>>(
    conn: &Connection,
    chunk: usize,
    mut visit: impl FnMut(&[PsmPeaks]) -> Result<(), E>,
) -> Result<(), E> {
    let mut statement = conn
        .prepare(
            "SELECT p.id, p.peptide_id, p.charge, s.peaks FROM psms p
             JOIN spectra s ON s.id = p.spectrum_id
             WHERE p.rank = 1 AND s.peaks IS NOT NULL ORDER BY p.id",
        )
        .map_err(StoreError::from)?;
    let mut rows = statement.query([]).map_err(StoreError::from)?;
    let mut batch = Vec::with_capacity(chunk);
    while let Some(row) = rows.next().map_err(StoreError::from)? {
        let psm: i64 = row.get(0).map_err(StoreError::from)?;
        let blob: Vec<u8> = row.get(3).map_err(StoreError::from)?;
        let peaks = peaks::decode(&blob)
            .ok_or_else(|| StoreError::Damaged(format!("the peaks of PSM {psm}'s spectrum")))?;
        batch.push(PsmPeaks {
            psm,
            peptide: row.get(1).map_err(StoreError::from)?,
            charge: row.get(2).map_err(StoreError::from)?,
            peaks,
        });
        if batch.len() == chunk {
            visit(&batch)?;
            batch.clear();
        }
    }
    if !batch.is_empty() {
        visit(&batch)?;
    }
    Ok(())
}

/// Visit the stored peaks (if any) of the requested spectra one at a time.
pub(crate) fn for_each_spectrum_peaks(
    conn: &Connection,
    spectrum_ids: &[i64],
    mut visit: impl FnMut(i64, Option<&Peaks>),
) -> Result<(), StoreError> {
    if spectrum_ids.is_empty() {
        return Ok(());
    }
    let mut statement = conn.prepare("SELECT peaks FROM spectra WHERE id = ?1")?;
    for &spectrum_id in spectrum_ids {
        let blob: Option<Vec<u8>> = statement
            .query_row([spectrum_id], |row| row.get(0))
            .optional()?
            .flatten();
        match blob {
            None => visit(spectrum_id, None),
            Some(blob) => {
                let peaks = peaks::decode(&blob).ok_or_else(|| {
                    StoreError::Damaged(format!("the peaks of spectrum {spectrum_id}"))
                })?;
                visit(spectrum_id, Some(&peaks));
            }
        }
    }
    Ok(())
}
