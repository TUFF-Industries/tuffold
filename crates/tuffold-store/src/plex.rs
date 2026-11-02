//! Isobaric plexes (schema v9): which runs make up a TMT plex (one run, or
//! the fractions of one labelled pool) and which biosample each reporter
//! channel holds. Channel biosamples are ordinary biosamples (categories
//! apply as usual) that have no runs of their own.

use std::collections::HashMap;

use rusqlite::{Connection, OptionalExtension, params};
use tuffold_core::reporter::{Impurities, Plex};
use tuffold_quant::reporter::{correct, extract, extract_summed};
use tuffold_validate::Level;

use crate::browse::Biosample;
use crate::{StoreError, peaks, validation};

/// A plex as stored or to be stored.
#[derive(Debug, Clone, PartialEq)]
pub struct PlexDefinition {
    pub name: String,
    pub plex: Plex,
    pub runs: Vec<String>,
    /// (channel, biosample name), in the plex's channel order.
    pub channels: Vec<(String, String)>,
    /// The reagent lot's isotope impurities, when given; reporters are
    /// corrected with them.
    pub impurities: Option<Impurities>,
    /// Channels holding the pooled reference sample that bridges this plex
    /// to others (in channel order).
    pub references: Vec<String>,
}

/// One protein group's passing PSMs in a plex, with their reporter
/// intensities (parallel to [`ReporterList::channels`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ReporterListRow {
    pub group_id: i64,
    pub accession: String,
    pub description: Option<String>,
    pub psms: Vec<Vec<Option<f64>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReporterList {
    pub plex: PlexDefinition,
    /// Each assigned channel's biosample, in channel order.
    pub biosamples: Vec<Biosample>,
    pub rows: Vec<ReporterListRow>,
    /// Passing spectra of the plex's groups that have no peaks to read.
    pub without_peaks: usize,
    /// PSMs whose reporters were read from their SPS-MS3 scan rather than
    /// the MS/MS spectrum itself.
    pub from_ms3: usize,
    /// Passing spectra left out because their precursor purity was below
    /// the cut, and ones whose purity is unknown (no MS1 scans or isolation
    /// window), which are kept.
    pub impure: usize,
    pub purity_unknown: usize,
    /// For each reagent set the plex could have been labelled with instead
    /// ([`Plex::alternatives`]), how many of these PSMs would have every
    /// assigned channel at its m/z.
    pub alternatives: Vec<(Plex, usize)>,
}

pub(crate) fn set_plex(
    conn: &mut Connection,
    definition: &PlexDefinition,
) -> Result<(), StoreError> {
    if definition.name.trim().is_empty() {
        return Err(StoreError::UnnamedPlex);
    }
    if definition.channels.is_empty() {
        return Err(StoreError::NoChannels);
    }
    if let Some(unassigned) = definition
        .references
        .iter()
        .find(|r| !definition.channels.iter().any(|(c, _)| c == *r))
    {
        return Err(StoreError::UnassignedReference(unassigned.clone()));
    }
    if let Some(given) = definition
        .impurities
        .as_ref()
        .map(Impurities::plex)
        .filter(|&given| given != definition.plex)
    {
        return Err(StoreError::ImpuritiesForOtherPlex {
            given: given.name(),
            plex: definition.plex.name(),
        });
    }
    let known: Vec<String> = definition
        .plex
        .channels()
        .into_iter()
        .map(|c| c.name)
        .collect();
    let mut seen = Vec::new();
    for (channel, _) in &definition.channels {
        if !known.contains(channel) {
            return Err(StoreError::UnknownChannel {
                channel: channel.clone(),
                plex: definition.plex.name(),
                known,
            });
        }
        if seen.contains(channel) {
            return Err(StoreError::DuplicateChannel(channel.clone()));
        }
        seen.push(channel.clone());
    }
    let tx = conn.transaction()?;
    let mut run_ids = Vec::new();
    for run in &definition.runs {
        let id: Option<i64> = tx
            .query_row("SELECT id FROM ms_runs WHERE name = ?1", [run], |row| {
                row.get(0)
            })
            .optional()?;
        let Some(id) = id else {
            let known = tx
                .prepare("SELECT name FROM ms_runs ORDER BY name")?
                .query_map([], |row| row.get(0))?
                .collect::<Result<Vec<String>, _>>()?;
            return Err(StoreError::UnknownRun {
                run: run.clone(),
                known,
            });
        };
        let other: Option<String> = tx
            .query_row(
                "SELECT x.name FROM plex_runs r JOIN plexes x ON x.id = r.plex_id
                 WHERE r.ms_run_id = ?1 AND x.name != ?2",
                params![id, definition.name],
                |row| row.get(0),
            )
            .optional()?;
        if let Some(plex) = other {
            return Err(StoreError::RunInPlex {
                run: run.clone(),
                plex,
            });
        }
        run_ids.push(id);
    }
    tx.execute("DELETE FROM plexes WHERE name = ?1", [&definition.name])?;
    tx.execute(
        "INSERT INTO plexes (name, label, impurities) VALUES (?1, ?2, ?3)",
        params![
            definition.name,
            definition.plex.name(),
            definition.impurities.as_ref().map(Impurities::to_text)
        ],
    )?;
    let plex_id = tx.last_insert_rowid();
    for id in run_ids {
        tx.execute(
            "INSERT INTO plex_runs (ms_run_id, plex_id) VALUES (?1, ?2)",
            params![id, plex_id],
        )?;
    }
    for (channel, sample) in &definition.channels {
        tx.execute(
            "INSERT OR IGNORE INTO biosamples (name) VALUES (?1)",
            [sample],
        )?;
        tx.execute(
            "INSERT INTO plex_channels (plex_id, channel, biosample_id, is_reference)
             SELECT ?1, ?2, id, ?4 FROM biosamples WHERE name = ?3",
            params![
                plex_id,
                channel,
                sample,
                definition.references.contains(channel)
            ],
        )?;
    }
    tx.commit()?;
    Ok(())
}

/// Remove a plex: its run and channel rows go with it (`ON DELETE
/// CASCADE`); its channel biosamples stay, as ordinary biosamples.
pub(crate) fn remove_plex(conn: &Connection, name: &str) -> Result<(), StoreError> {
    if conn.execute("DELETE FROM plexes WHERE name = ?1", [name])? == 0 {
        let known = conn
            .prepare("SELECT name FROM plexes ORDER BY name")?
            .query_map([], |row| row.get(0))?
            .collect::<Result<Vec<String>, _>>()?;
        return Err(StoreError::UnknownPlex {
            name: name.to_string(),
            known,
        });
    }
    Ok(())
}

pub(crate) fn plexes(conn: &Connection) -> Result<Vec<PlexDefinition>, StoreError> {
    let stored: Vec<(i64, String, String, Option<String>)> = conn
        .prepare("SELECT id, name, label, impurities FROM plexes ORDER BY name")?
        .query_map([], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?))
        })?
        .collect::<Result<_, _>>()?;
    let mut out = Vec::new();
    for (id, name, label, impurities) in stored {
        let plex = Plex::parse(&label).map_err(StoreError::Damaged)?;
        let impurities = impurities
            .map(|text| Impurities::parse(&text, plex))
            .transpose()
            .map_err(|e| StoreError::Damaged(format!("plex {name}'s impurities: {e}")))?;
        let runs = conn
            .prepare(
                "SELECT m.name FROM plex_runs r JOIN ms_runs m ON m.id = r.ms_run_id
                 WHERE r.plex_id = ?1 ORDER BY m.name",
            )?
            .query_map([id], |row| row.get(0))?
            .collect::<Result<Vec<String>, _>>()?;
        let assigned: HashMap<String, (String, bool)> = conn
            .prepare(
                "SELECT c.channel, b.name, c.is_reference FROM plex_channels c
                 JOIN biosamples b ON b.id = c.biosample_id WHERE c.plex_id = ?1",
            )?
            .query_map([id], |row| Ok((row.get(0)?, (row.get(1)?, row.get(2)?))))?
            .collect::<Result<_, _>>()?;
        let order = plex.channels();
        let channels = order
            .iter()
            .filter_map(|c| Some((c.name.clone(), assigned.get(&c.name)?.0.clone())))
            .collect();
        let references = order
            .iter()
            .filter(|c| {
                assigned
                    .get(&c.name)
                    .is_some_and(|(_, reference)| *reference)
            })
            .map(|c| c.name.clone())
            .collect();
        out.push(PlexDefinition {
            name,
            plex,
            runs,
            channels,
            impurities,
            references,
        });
    }
    Ok(out)
}

/// Protein groups passing the saved thresholds, with their unique and razor
/// peptides' passing spectra in the plex's runs and those spectra's reporter
/// intensities (in accession order).
pub(crate) fn reporter_list(
    conn: &Connection,
    plex_name: &str,
    min_purity: Option<f64>,
) -> Result<ReporterList, StoreError> {
    let all = plexes(conn)?;
    let Some(plex) = all.iter().find(|p| p.name == plex_name).cloned() else {
        return Err(StoreError::UnknownPlex {
            name: plex_name.to_string(),
            known: all.into_iter().map(|p| p.name).collect(),
        });
    };
    let tolerance_ppm = plex.plex.tolerance_ppm();
    let read = if plex.plex.is_itraq() {
        extract_summed
    } else {
        extract
    };
    let plex_id: i64 = conn.query_row(
        "SELECT id FROM plexes WHERE name = ?1",
        [plex_name],
        |row| row.get(0),
    )?;
    let mz: HashMap<String, f64> = plex
        .plex
        .channels()
        .into_iter()
        .map(|c| (c.name, c.mz))
        .collect();
    let targets: Vec<f64> = plex.channels.iter().map(|(c, _)| mz[c]).collect();
    let unmixing = plex.impurities.as_ref().map(|impurities| {
        let all = plex.plex.channels();
        let assigned: Vec<usize> = plex
            .channels
            .iter()
            .filter_map(|(c, _)| all.iter().position(|a| &a.name == c))
            .collect();
        let mz: Vec<f64> = all.iter().map(|c| c.mz).collect();
        (impurities.matrix(), mz, assigned)
    });
    let alternatives: Vec<(Plex, Vec<f64>)> = plex
        .plex
        .alternatives()
        .into_iter()
        .filter_map(|other| {
            let mz: HashMap<String, f64> = other
                .channels()
                .into_iter()
                .map(|c| (c.name, c.mz))
                .collect();
            let targets = plex
                .channels
                .iter()
                .map(|(c, _)| mz.get(c).copied())
                .collect::<Option<Vec<f64>>>()?;
            Some((other, targets))
        })
        .collect();
    let mut complete_elsewhere = vec![0; alternatives.len()];
    let biosamples = plex
        .channels
        .iter()
        .map(|(_, sample)| {
            conn.query_row(
                "SELECT id, name, category FROM biosamples WHERE name = ?1",
                [sample],
                |row| {
                    Ok(Biosample {
                        id: row.get(0)?,
                        name: row.get(1)?,
                        category: row.get(2)?,
                    })
                },
            )
            .map_err(Into::into)
        })
        .collect::<Result<Vec<_>, StoreError>>()?;
    let thresholds = validation::settings(conn)?
        .map(|(settings, _)| settings.thresholds)
        .unwrap_or_default();
    let measure = validation::measure_column(thresholds.mode);
    let mut rows: Vec<ReporterListRow> = conn
        .prepare(&format!(
            "SELECT g.id, p.accession, p.description
             FROM protein_groups g JOIN proteins p ON p.id = g.leader_protein_id
             WHERE g.is_decoy = 0 AND g.{measure} <= ?1 AND g.peptides >= ?2"
        ))?
        .query_map(
            params![thresholds.bound(Level::Protein), thresholds.min_peptides],
            |row| {
                Ok(ReporterListRow {
                    group_id: row.get(0)?,
                    accession: row.get(1)?,
                    description: row.get(2)?,
                    psms: Vec::new(),
                })
            },
        )?
        .collect::<Result<_, _>>()?;
    let index: HashMap<i64, usize> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (row.group_id, i))
        .collect();
    let mut statement = conn.prepare(&format!(
        "SELECT gp.group_id, s.peaks, s.ms3_peaks, s.precursor_purity
         FROM protein_group_peptides gp
         JOIN spectrum_results sr ON sr.peptide_id = gp.peptide_id
         JOIN spectra s ON s.id = sr.spectrum_id
         JOIN plex_runs pr ON pr.ms_run_id = s.ms_run_id
         WHERE gp.role IN ('unique', 'razor') AND pr.plex_id = ?1
           AND sr.is_decoy = 0 AND sr.{measure} <= ?2"
    ))?;
    let mut found = statement.query(params![plex_id, thresholds.bound(Level::Psm)])?;
    let (mut without_peaks, mut from_ms3, mut impure, mut purity_unknown) = (0, 0, 0, 0);
    while let Some(row) = found.next()? {
        let Some(&i) = index.get(&row.get::<_, i64>(0)?) else {
            continue;
        };
        let purity: Option<f64> = row.get(3)?;
        // SPS-MS3: the reporters are in the MS3 scan acquired from this one.
        let ms3 = row
            .get::<_, Option<Vec<u8>>>(2)?
            .and_then(|blob| peaks::decode(&blob));
        if ms3.is_some() {
            from_ms3 += 1;
        }
        let Some(peaks) = ms3.or(row
            .get::<_, Option<Vec<u8>>>(1)?
            .and_then(|blob| peaks::decode(&blob)))
        else {
            without_peaks += 1;
            continue;
        };
        // Among spectra with peaks to read: some left out as impure would
        // also have missed a channel (counted only as impure).
        match (min_purity, purity) {
            (Some(min), Some(purity)) if purity < min => {
                impure += 1;
                continue;
            }
            (_, None) => purity_unknown += 1,
            _ => {}
        }
        for ((_, other), count) in alternatives.iter().zip(&mut complete_elsewhere) {
            let found = extract(&peaks.mz, &peaks.intensity, other, tolerance_ppm);
            if found.iter().all(Option::is_some) {
                *count += 1;
            }
        }
        let reporters = match &unmixing {
            None => read(&peaks.mz, &peaks.intensity, &targets, tolerance_ppm),
            Some((matrix, mz, assigned)) => corrected(
                &read(&peaks.mz, &peaks.intensity, mz, tolerance_ppm),
                matrix,
                assigned,
            ),
        };
        rows[i].psms.push(reporters);
    }
    rows.retain(|row| !row.psms.is_empty());
    rows.sort_by(|a, b| a.accession.cmp(&b.accession));
    Ok(ReporterList {
        plex,
        biosamples,
        rows,
        without_peaks,
        from_ms3,
        impure,
        purity_unknown,
        alternatives: alternatives
            .into_iter()
            .map(|(other, _)| other)
            .zip(complete_elsewhere)
            .collect(),
    })
}

/// A PSM's assigned channels after unmixing all of the reagent set's
/// (`all`, in channel order; an unassigned channel without a reporter holds
/// nothing). A PSM missing an assigned channel is left as read: it is left
/// out of quantitation anyway.
fn corrected(all: &[Option<f64>], matrix: &[Vec<f64>], assigned: &[usize]) -> Vec<Option<f64>> {
    if assigned.iter().any(|&i| all[i].is_none()) {
        return assigned.iter().map(|&i| all[i]).collect();
    }
    let observed: Vec<f64> = all.iter().map(|v| v.unwrap_or(0.0)).collect();
    let unmixed = correct(&observed, matrix);
    assigned.iter().map(|&i| Some(unmixed[i])).collect()
}
