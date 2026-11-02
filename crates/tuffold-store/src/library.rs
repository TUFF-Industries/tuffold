//! What a spectral library is built from: for every peptide ion (peptide
//! form at one charge) with a passing target spectrum, its best-scoring
//! passing spectrum that has peaks attached, and how many passing spectra it
//! had in all.

use std::collections::HashMap;

use rusqlite::{Connection, params};
use tuffold_core::PeptideForm;
use tuffold_validate::Level;

use crate::browse::peptide_form;
use crate::peaks::{self, Peaks};
use crate::{StoreError, validation};

#[derive(Debug, Clone, PartialEq)]
pub struct LibraryEntry {
    pub form: PeptideForm,
    pub charge: i64,
    pub run: String,
    /// `scan=N`, `index=N`, or a native id.
    pub key: String,
    pub rt_seconds: Option<f64>,
    pub q_value: f64,
    pub score: f64,
    /// Passing target spectra of this peptide ion, the chosen one included.
    pub replicates: usize,
    /// The leader of the passing group holding it as unique or razor; for
    /// a peptide in no passing group, the match's first target protein.
    pub protein: Option<String>,
    /// The residues before and after the peptide (`-` at a protein end):
    /// where it sits in the leader's sequence (else another target protein
    /// of the match), else as the engine reported them; `None` when nothing
    /// says.
    pub flanks: Option<(char, char)>,
    pub peaks: Peaks,
}

/// Library entries, by best score; and how many peptide ions passed without
/// any spectrum with peaks (left out).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Library {
    pub entries: Vec<LibraryEntry>,
    pub without_peaks: usize,
}

struct Candidate {
    spectrum_id: i64,
    run: String,
    key: String,
    rt_seconds: Option<f64>,
    q_value: f64,
    score: f64,
    has_peaks: bool,
}

pub(crate) fn library(conn: &Connection) -> Result<Library, StoreError> {
    let thresholds = validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default();
    let measure = validation::measure_column(thresholds.mode);
    let mut statement = conn.prepare(&format!(
        "SELECT sr.peptide_id,
                (SELECT MIN(p.charge) FROM psms p
                  WHERE p.spectrum_id = s.id AND p.peptide_id = sr.peptide_id AND p.rank = 1),
                s.id, r.name, s.spectrum_key, s.rt_seconds, sr.q_value, sr.score,
                s.peaks IS NOT NULL
         FROM spectrum_results sr
         JOIN spectra s ON s.id = sr.spectrum_id
         JOIN ms_runs r ON r.id = s.ms_run_id
         WHERE sr.is_decoy = 0 AND sr.{measure} <= ?1"
    ))?;
    let mut ions: HashMap<(i64, i64), Vec<Candidate>> = HashMap::new();
    let mut rows = statement.query(params![thresholds.bound(Level::Psm)])?;
    while let Some(row) = rows.next()? {
        let Some(charge) = row.get::<_, Option<i64>>(1)? else {
            continue;
        };
        ions.entry((row.get(0)?, charge))
            .or_default()
            .push(Candidate {
                spectrum_id: row.get(2)?,
                run: row.get(3)?,
                key: row.get(4)?,
                rt_seconds: row.get(5)?,
                q_value: row.get(6)?,
                score: row.get(7)?,
                has_peaks: row.get(8)?,
            });
    }
    let proteins = leaders(conn, &thresholds)?;
    let mut library = Library::default();
    let mut keys: Vec<(i64, i64)> = ions.keys().copied().collect();
    keys.sort_unstable();
    for (peptide_id, charge) in keys {
        let candidates = &ions[&(peptide_id, charge)];
        let Some(best) = candidates.iter().filter(|c| c.has_peaks).max_by(|a, b| {
            a.score
                .total_cmp(&b.score)
                .then(b.spectrum_id.cmp(&a.spectrum_id))
        }) else {
            library.without_peaks += 1;
            continue;
        };
        let blob: Vec<u8> = conn.query_row(
            "SELECT peaks FROM spectra WHERE id = ?1",
            [best.spectrum_id],
            |row| row.get(0),
        )?;
        let Some(peaks) = peaks::decode(&blob) else {
            library.without_peaks += 1;
            continue;
        };
        let form = peptide_form(conn, peptide_id)?;
        let leader = proteins.get(&peptide_id);
        let flanks = match leader
            .and_then(|(_, sequence)| sequence.as_deref())
            .and_then(|sequence| flanks_in(sequence, &form.sequence))
        {
            Some(flanks) => Some(flanks),
            None => other_flanks(conn, best.spectrum_id, peptide_id, &form.sequence)?,
        };
        library.entries.push(LibraryEntry {
            form,
            charge,
            run: best.run.clone(),
            key: best.key.clone(),
            rt_seconds: best.rt_seconds,
            q_value: best.q_value,
            score: best.score,
            replicates: candidates.len(),
            protein: match leader {
                Some((accession, _)) => Some(accession.clone()),
                None => first_target(conn, best.spectrum_id, peptide_id)?,
            },
            flanks,
            peaks,
        });
    }
    library.entries.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| a.form.proforma().cmp(&b.form.proforma()))
            .then(a.charge.cmp(&b.charge))
    });
    Ok(library)
}

/// The residues around the first place `peptide` occurs in `protein` (I
/// and L alike), `-` at an end.
fn flanks_in(protein: &str, peptide: &str) -> Option<(char, char)> {
    let leucine = |text: &str| text.to_ascii_uppercase().replace('I', "L");
    let (haystack, needle) = (leucine(protein), leucine(peptide));
    let start = haystack.find(&needle)?;
    let end = start + needle.len();
    let before = start
        .checked_sub(1)
        .map_or('-', |i| char::from(protein.as_bytes()[i]));
    let after = protein.as_bytes().get(end).map_or('-', |&b| char::from(b));
    Some((before.to_ascii_uppercase(), after.to_ascii_uppercase()))
}

/// Flanks for a peptide whose leader doesn't say (no passing group, or no
/// sequence): from the first target protein of this spectrum's match that
/// has a sequence holding it, else as the engine reported them.
fn other_flanks(
    conn: &Connection,
    spectrum_id: i64,
    peptide_id: i64,
    peptide: &str,
) -> Result<Option<(char, char)>, StoreError> {
    let mut statement = conn.prepare(
        "SELECT pr.sequence, pp.pre, pp.post FROM psms p
         JOIN psm_proteins pp ON pp.psm_id = p.id
         JOIN proteins pr ON pr.id = pp.protein_id
         WHERE p.spectrum_id = ?1 AND p.peptide_id = ?2 AND pr.is_decoy = 0
         ORDER BY p.rank, pr.id",
    )?;
    let rows = statement.query_map(params![spectrum_id, peptide_id], |row| {
        Ok((
            row.get::<_, Option<String>>(0)?,
            row.get::<_, Option<String>>(1)?,
            row.get::<_, Option<String>>(2)?,
        ))
    })?;
    let mut reported = None;
    for row in rows {
        let (sequence, pre, post) = row?;
        if let Some(flanks) = sequence.as_deref().and_then(|s| flanks_in(s, peptide)) {
            return Ok(Some(flanks));
        }
        let first = |text: Option<String>| text.and_then(|t| t.chars().next());
        if reported.is_none()
            && let (Some(pre), Some(post)) = (first(pre), first(post))
        {
            reported = Some((pre, post));
        }
    }
    Ok(reported)
}

/// The first target protein this spectrum's match maps to (for a peptide
/// in no passing group).
fn first_target(
    conn: &Connection,
    spectrum_id: i64,
    peptide_id: i64,
) -> Result<Option<String>, StoreError> {
    let mut statement = conn.prepare(
        "SELECT pr.accession FROM psms p
         JOIN psm_proteins pp ON pp.psm_id = p.id
         JOIN proteins pr ON pr.id = pp.protein_id
         WHERE p.spectrum_id = ?1 AND p.peptide_id = ?2 AND pr.is_decoy = 0
         ORDER BY p.rank, pr.id LIMIT 1",
    )?;
    let mut rows = statement.query(params![spectrum_id, peptide_id])?;
    Ok(match rows.next()? {
        Some(row) => Some(row.get(0)?),
        None => None,
    })
}

/// Each peptide's passing unique or razor group's leader: accession and
/// sequence.
fn leaders(
    conn: &Connection,
    thresholds: &tuffold_validate::Thresholds,
) -> Result<HashMap<i64, (String, Option<String>)>, StoreError> {
    let measure = validation::measure_column(thresholds.mode);
    let mut statement = conn.prepare(&format!(
        "SELECT gp.peptide_id, p.accession, p.sequence
         FROM protein_group_peptides gp
         JOIN protein_groups g ON g.id = gp.group_id
         JOIN proteins p ON p.id = g.leader_protein_id
         WHERE gp.role IN ('unique', 'razor') AND g.is_decoy = 0
           AND g.{measure} <= ?1 AND g.peptides >= ?2"
    ))?;
    let leaders = statement
        .query_map(
            params![thresholds.bound(Level::Protein), thresholds.min_peptides],
            |row| Ok((row.get(0)?, (row.get(1)?, row.get(2)?))),
        )?
        .collect::<Result<HashMap<_, _>, _>>()?;
    Ok(leaders)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flanks_come_from_where_the_peptide_sits() {
        assert_eq!(flanks_in("MKPEPTIDEKSAM", "PEPTIDEK"), Some(('K', 'S')));
        // At the protein's ends: '-'.
        assert_eq!(flanks_in("PEPTIDEK", "PEPTIDEK"), Some(('-', '-')));
        // I and L alike; the protein's own residues are reported.
        assert_eq!(flanks_in("RLEAKR", "IEAK"), Some(('R', 'R')));
        assert_eq!(flanks_in("MKAAA", "PEPTIDEK"), None);
    }
}
