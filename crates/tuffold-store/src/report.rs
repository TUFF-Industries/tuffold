//! Whole-project listings for reports and exports: every peptide and every
//! spectrum passing the saved thresholds (the browsing queries page through
//! one protein or peptide at a time instead).

use std::collections::HashMap;

use rusqlite::{Connection, params};
use tuffold_validate::{Level, SiteScore, format_sites};

use crate::browse::{Biosample, biosamples};
use crate::{StoreError, validation};

#[derive(Debug, Clone, PartialEq)]
pub struct PeptideReportRow {
    pub peptide_id: i64,
    pub proforma: String,
    pub sequence: String,
    pub q_value: f64,
    /// 1 − PEP.
    pub probability: f64,
    pub score: f64,
    /// Passing target spectra assigned this peptide.
    pub spectra: usize,
    /// Parallel to [`PeptideReport::biosamples`].
    pub spectra_by_sample: Vec<usize>,
    /// Passing protein groups holding it: (leader accession, role), the
    /// unique or razor group first.
    pub groups: Vec<(String, String)>,
    /// The peptide's localization is that of its passing spectrum with the
    /// highest minimum site probability (the best spectrum's), not a combined
    /// estimate over its spectra.
    pub localized_peptide: Option<String>,
    pub site_probabilities: Option<String>,
    pub best_ascore: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumReportRow {
    pub run: String,
    pub biosample: String,
    /// `scan=N`, `index=N`, or a native id.
    pub key: String,
    pub charge: Option<i64>,
    pub precursor_mz: Option<f64>,
    pub rt_seconds: Option<f64>,
    pub proforma: String,
    pub q_value: f64,
    pub probability: f64,
    pub score: f64,
    pub engines: u32,
    pub localized_peptide: Option<String>,
    pub site_probabilities: Option<String>,
    pub best_ascore: Option<f64>,
}

/// Every passing target peptide and spectrum.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Report {
    pub biosamples: Vec<Biosample>,
    /// Best score first.
    pub peptides: Vec<PeptideReportRow>,
    /// By run, then spectrum.
    pub spectra: Vec<SpectrumReportRow>,
}

fn count(value: i64) -> usize {
    usize::try_from(value).unwrap_or_default()
}

/// What the report prints of a scored spectrum's localization.
#[derive(Clone)]
struct ReportLocalization {
    localized_peptide: String,
    site_probabilities: String,
    best_ascore: f64,
}

/// The localization of every passing target spectrum that was scored (and
/// has all of its summary values), in a fixed number of queries; the same
/// values `browse::spectrum_localization` reads one spectrum at a time.
fn spectrum_localizations(
    conn: &Connection,
    measure: &str,
    psm_bound: f64,
) -> Result<HashMap<i64, ReportLocalization>, StoreError> {
    let mut sites: HashMap<i64, Vec<SiteScore>> = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT t.spectrum_id, t.position, t.residue, t.mass_delta, t.modification_name,
                t.probability, t.ascore, t.determining_matched, t.determining_total, t.assigned
         FROM ptm_site_scores t
         JOIN ptm_spectrum_summary ps ON ps.spectrum_id = t.spectrum_id
         JOIN spectrum_results sr ON sr.spectrum_id = t.spectrum_id
         WHERE t.is_decoy = 0 AND ps.status = 'scored'
           AND sr.is_decoy = 0 AND sr.{measure} <= ?1
         ORDER BY t.spectrum_id, t.position, t.mass_delta"
    ))?;
    let mut rows = statement.query([psm_bound])?;
    while let Some(row) = rows.next()? {
        let residue: String = row.get(2)?;
        sites.entry(row.get(0)?).or_default().push(SiteScore {
            position: row.get(1)?,
            residue: residue.chars().next().unwrap_or('?'),
            mass_delta: row.get(3)?,
            modification_name: row.get(4)?,
            probability: row.get(5)?,
            ascore: row.get(6)?,
            determining_ions_matched: row.get(7)?,
            determining_ions_total: row.get(8)?,
            assigned: row.get::<_, i64>(9)? != 0,
            is_decoy: false,
        });
    }

    let mut localizations = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT ps.spectrum_id, ps.localized_peptide, ps.best_ascore, ps.min_probability
         FROM ptm_spectrum_summary ps
         JOIN spectrum_results sr ON sr.spectrum_id = ps.spectrum_id
         WHERE ps.status = 'scored' AND sr.is_decoy = 0 AND sr.{measure} <= ?1"
    ))?;
    let mut rows = statement.query([psm_bound])?;
    while let Some(row) = rows.next()? {
        let spectrum_id: i64 = row.get(0)?;
        let summary: (Option<String>, Option<f64>, Option<f64>) =
            (row.get(1)?, row.get(2)?, row.get(3)?);
        let (Some(localized_peptide), Some(best_ascore), Some(_)) = summary else {
            continue;
        };
        let spectrum_sites = sites.remove(&spectrum_id).unwrap_or_default();
        localizations.insert(
            spectrum_id,
            ReportLocalization {
                localized_peptide,
                site_probabilities: format_sites(&spectrum_sites),
                best_ascore,
            },
        );
    }
    Ok(localizations)
}

/// Each peptide's best scored passing spectrum: the highest minimum site
/// probability, then Ascore, then spectrum score, then the lower spectrum id
/// (the order `browse::best_peptide_localization` uses).
fn best_localized_spectra(
    conn: &Connection,
    measure: &str,
    psm_bound: f64,
) -> Result<HashMap<i64, i64>, StoreError> {
    let mut best = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT sr.peptide_id, ps.spectrum_id
         FROM spectrum_results sr
         JOIN ptm_spectrum_summary ps ON ps.spectrum_id = sr.spectrum_id
         WHERE sr.is_decoy = 0 AND sr.{measure} <= ?1 AND ps.status = 'scored'
         ORDER BY sr.peptide_id, ps.min_probability DESC, ps.best_ascore DESC,
                  sr.score DESC, sr.spectrum_id"
    ))?;
    let mut rows = statement.query([psm_bound])?;
    while let Some(row) = rows.next()? {
        best.entry(row.get::<_, i64>(0)?)
            .or_insert(row.get::<_, i64>(1)?);
    }
    Ok(best)
}

pub(crate) fn report(conn: &Connection) -> Result<Report, StoreError> {
    let thresholds = validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default();
    let measure = validation::measure_column(thresholds.mode);
    let psm_bound = thresholds.bound(Level::Psm);
    let biosamples = biosamples(conn)?;
    let column: HashMap<i64, usize> = biosamples
        .iter()
        .enumerate()
        .map(|(i, b)| (b.id, i))
        .collect();

    let mut by_sample: HashMap<i64, Vec<usize>> = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT sr.peptide_id, r.biosample_id, COUNT(*)
         FROM spectrum_results sr
         JOIN spectra s ON s.id = sr.spectrum_id
         JOIN ms_runs r ON r.id = s.ms_run_id
         WHERE sr.is_decoy = 0 AND sr.{measure} <= ?1
         GROUP BY sr.peptide_id, r.biosample_id"
    ))?;
    let mut rows = statement.query([psm_bound])?;
    while let Some(row) = rows.next()? {
        let counts = by_sample
            .entry(row.get(0)?)
            .or_insert_with(|| vec![0; biosamples.len()]);
        if let Some(&s) = column.get(&row.get::<_, i64>(1)?) {
            counts[s] = count(row.get(2)?);
        }
    }

    let mut groups: HashMap<i64, Vec<(String, String)>> = HashMap::new();
    let mut statement = conn.prepare(&format!(
        "SELECT gp.peptide_id, p.accession, gp.role
         FROM protein_group_peptides gp
         JOIN protein_groups g ON g.id = gp.group_id
         JOIN proteins p ON p.id = g.leader_protein_id
         WHERE g.is_decoy = 0 AND g.{measure} <= ?1 AND g.peptides >= ?2
         ORDER BY gp.peptide_id, CASE gp.role WHEN 'shared' THEN 1 ELSE 0 END, p.accession"
    ))?;
    let mut rows = statement.query(params![
        thresholds.bound(Level::Protein),
        thresholds.min_peptides
    ])?;
    while let Some(row) = rows.next()? {
        groups
            .entry(row.get(0)?)
            .or_default()
            .push((row.get(1)?, row.get(2)?));
    }

    let mut localizations = spectrum_localizations(conn, measure, psm_bound)?;
    let best_spectra = best_localized_spectra(conn, measure, psm_bound)?;

    let peptide_rows = conn
        .prepare(&format!(
            "SELECT pe.id, pe.proforma, pe.sequence, r.q_value, r.pep, r.score
             FROM peptide_results r JOIN peptides pe ON pe.id = r.peptide_id
             WHERE r.is_decoy = 0 AND r.{measure} <= ?1
             ORDER BY r.score DESC, pe.proforma"
        ))?
        .query_map([thresholds.bound(Level::Peptide)], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                row.get::<_, String>(1)?,
                row.get::<_, String>(2)?,
                row.get::<_, f64>(3)?,
                row.get::<_, f64>(4)?,
                row.get::<_, f64>(5)?,
            ))
        })?
        .collect::<Result<Vec<_>, _>>()?;

    let mut peptides = Vec::with_capacity(peptide_rows.len());
    for (peptide_id, proforma, sequence, q_value, pep, score) in peptide_rows {
        let spectra_by_sample = by_sample
            .remove(&peptide_id)
            .unwrap_or_else(|| vec![0; biosamples.len()]);
        let loc = best_spectra
            .get(&peptide_id)
            .and_then(|spectrum_id| localizations.get(spectrum_id));
        let (localized_peptide, site_probabilities, best_ascore) = match loc {
            Some(loc) => (
                Some(loc.localized_peptide.clone()),
                Some(loc.site_probabilities.clone()),
                Some(loc.best_ascore),
            ),
            None => (None, None, None),
        };
        peptides.push(PeptideReportRow {
            peptide_id,
            proforma,
            sequence,
            q_value,
            probability: 1.0 - pep,
            score,
            spectra: spectra_by_sample.iter().sum(),
            spectra_by_sample,
            groups: groups.remove(&peptide_id).unwrap_or_default(),
            localized_peptide,
            site_probabilities,
            best_ascore,
        });
    }

    let spectrum_rows = conn
        .prepare(&format!(
            "SELECT s.id, r.name, b.name, s.spectrum_key, s.precursor_charge, s.precursor_mz,
                    s.rt_seconds, pe.proforma, sr.q_value, sr.pep, sr.score, sr.engines
             FROM spectrum_results sr
             JOIN spectra s ON s.id = sr.spectrum_id
             JOIN ms_runs r ON r.id = s.ms_run_id
             JOIN biosamples b ON b.id = r.biosample_id
             JOIN peptides pe ON pe.id = sr.peptide_id
             WHERE sr.is_decoy = 0 AND sr.{measure} <= ?1
             ORDER BY r.name, s.id"
        ))?
        .query_map([psm_bound], |row| {
            Ok((
                row.get::<_, i64>(0)?,
                SpectrumReportRow {
                    run: row.get(1)?,
                    biosample: row.get(2)?,
                    key: row.get(3)?,
                    charge: row.get(4)?,
                    precursor_mz: row.get(5)?,
                    rt_seconds: row.get(6)?,
                    proforma: row.get(7)?,
                    q_value: row.get(8)?,
                    probability: 1.0 - row.get::<_, f64>(9)?,
                    score: row.get(10)?,
                    engines: row.get(11)?,
                    localized_peptide: None,
                    site_probabilities: None,
                    best_ascore: None,
                },
            ))
        })?
        .collect::<Result<Vec<_>, _>>()?;

    let mut spectra = Vec::with_capacity(spectrum_rows.len());
    for (spectrum_id, mut row) in spectrum_rows {
        if let Some(loc) = localizations.remove(&spectrum_id) {
            row.localized_peptide = Some(loc.localized_peptide);
            row.site_probabilities = Some(loc.site_probabilities);
            row.best_ascore = Some(loc.best_ascore);
        }
        spectra.push(row);
    }

    Ok(Report {
        biosamples,
        peptides,
        spectra,
    })
}
