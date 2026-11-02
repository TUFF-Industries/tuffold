//! Several plexes on one scale: internal reference scaling (IRS; Plubell et
//! al., Mol Cell Proteomics 2017).
//!
//! Reporter intensities compare well within a plex but not between plexes:
//! which PSMs each run picked for MS/MS, and how intense they were, differ, so
//! one protein's sum lands on a different scale in every plex. When each plex
//! holds the same reference sample (a pool of all samples) in one or more
//! channels, that reference measures the scale: for each protein, every plex
//! is shifted so its reference reads the mean of the plexes' references
//! (log2, i.e. their geometric mean), and the protein's other channels move
//! with it. The reference channels then say nothing more and are left out of
//! the result.
//!
//! A plex whose reference channels have no value for a protein can't be put
//! on the common scale, so its values for that protein are left out and
//! counted. A protein that only one plex measured keeps its values (its
//! reference is the mean), comparable within that plex.

use std::collections::{HashMap, HashSet};

use crate::intensity::{IntensityQuantities, IntensityRow};
use crate::reporter::ReporterQuantities;
use crate::{QuantError, Sample};

/// One plex's quantitation ([`super::quantify_reporters`]) and which of its
/// samples are reference channels.
#[derive(Debug, Clone, PartialEq)]
pub struct PlexQuantities {
    pub name: String,
    pub quant: ReporterQuantities,
    /// Parallel to the plex's samples.
    pub reference: Vec<bool>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct Bridged {
    /// Every plex's non-reference channels, plex by plex; normalization shifts
    /// are each plex's own.
    pub quantities: IntensityQuantities,
    /// Each sample's plex and channel, as "plex channel".
    pub channels: Vec<String>,
    /// Which plex (index into the input) each sample came from.
    pub plex_of: Vec<usize>,
    /// (plex, protein) values left out because the plex's reference channels
    /// had no value for the protein.
    pub unbridged: usize,
    pub psms_used: usize,
    pub psms_incomplete: usize,
}

/// Put several plexes' proteins on one scale through their reference
/// channels (see the module notes).
pub fn bridge(plexes: Vec<PlexQuantities>) -> Result<Bridged, QuantError> {
    for plex in &plexes {
        if !plex.reference.contains(&true) {
            return Err(QuantError::NoReference(plex.name.clone()));
        }
        if !plex.reference.contains(&false) {
            return Err(QuantError::OnlyReference(plex.name.clone()));
        }
    }
    let mut out = Bridged::default();
    // Where each plex's kept samples start in the output.
    let mut offsets = Vec::new();
    for (index, plex) in plexes.iter().enumerate() {
        offsets.push(out.quantities.samples.len());
        let q = &plex.quant.quantities;
        for (i, sample) in q.samples.iter().enumerate() {
            if plex.reference[i] {
                continue;
            }
            out.quantities.samples.push(Sample {
                name: sample.name.clone(),
                category: sample.category.clone(),
            });
            out.quantities.shifts.push(q.shifts[i]);
            out.quantities.normalized.push(q.normalized[i]);
            out.channels
                .push(format!("{} {}", plex.name, plex.quant.channels[i]));
            out.plex_of.push(index);
        }
        out.psms_used += plex.quant.psms_used;
        out.psms_incomplete += plex.quant.psms_incomplete;
        out.quantities.shared_peptides += q.shared_peptides;
    }
    let n = out.quantities.samples.len();

    // Each protein, in first-seen order: its row in every plex measuring it.
    let mut order: Vec<(i64, String)> = Vec::new();
    let mut seen = HashSet::new();
    for plex in &plexes {
        for row in &plex.quant.quantities.rows {
            if seen.insert(row.group_id) {
                order.push((row.group_id, row.accession.clone()));
            }
        }
    }
    let rows_of: Vec<HashMap<i64, &IntensityRow>> = plexes
        .iter()
        .map(|plex| {
            plex.quant
                .quantities
                .rows
                .iter()
                .map(|r| (r.group_id, r))
                .collect()
        })
        .collect();
    for (group_id, accession) in order {
        let found: Vec<(usize, &IntensityRow, Option<f64>)> = plexes
            .iter()
            .zip(&rows_of)
            .enumerate()
            .filter_map(|(index, (plex, rows))| {
                let row = *rows.get(&group_id)?;
                Some((index, row, reference_level(row, &plex.reference)))
            })
            .collect();
        let levels: Vec<f64> = found.iter().filter_map(|(_, _, level)| *level).collect();
        out.unbridged += found.len() - levels.len();
        if levels.is_empty() {
            continue;
        }
        let target = levels.iter().sum::<f64>() / levels.len() as f64;
        let mut log2 = vec![None; n];
        let mut peptides = vec![0; n];
        for (index, row, level) in found {
            let Some(level) = level else { continue };
            let reference = &plexes[index].reference;
            let kept = (0..row.log2.len()).filter(|&i| !reference[i]);
            for (slot, i) in (offsets[index]..).zip(kept) {
                log2[slot] = row.log2[i].map(|v| v + target - level);
                peptides[slot] = row.peptides[i];
            }
        }
        out.quantities.rows.push(IntensityRow {
            group_id,
            accession,
            intensities: log2.iter().map(|v| v.map(f64::exp2)).collect(),
            log2,
            peptides,
            top3: vec![None; n],
        });
    }
    Ok(out)
}

/// The mean log2 of a protein's reference channels, if all have a value.
fn reference_level(row: &IntensityRow, reference: &[bool]) -> Option<f64> {
    let values: Vec<f64> = row
        .log2
        .iter()
        .zip(reference)
        .filter(|(_, is_reference)| **is_reference)
        .map(|(v, _)| *v)
        .collect::<Option<_>>()?;
    Some(values.iter().sum::<f64>() / values.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(name: &str) -> Sample {
        Sample {
            name: name.into(),
            category: None,
        }
    }

    fn row(group_id: i64, log2: [Option<f64>; 3]) -> IntensityRow {
        IntensityRow {
            group_id,
            accession: format!("P{group_id}"),
            intensities: log2.iter().map(|v| v.map(f64::exp2)).collect(),
            log2: log2.to_vec(),
            peptides: vec![2; 3],
            top3: vec![None; 3],
        }
    }

    /// A plex of a reference channel (126) and two samples.
    fn plex(name: &str, rows: Vec<IntensityRow>) -> PlexQuantities {
        PlexQuantities {
            name: name.into(),
            quant: ReporterQuantities {
                quantities: IntensityQuantities {
                    samples: vec![
                        sample(&format!("{name}-pool")),
                        sample(&format!("{name}-a")),
                        sample(&format!("{name}-b")),
                    ],
                    shifts: vec![0.0, 0.1, -0.1],
                    normalized: vec![true; 3],
                    shared_peptides: 5,
                    rows,
                },
                channels: vec!["126".into(), "127".into(), "128".into()],
                psms_used: 10,
                psms_incomplete: 1,
            },
            reference: vec![true, false, false],
        }
    }

    #[test]
    fn plexes_meet_at_the_mean_of_their_references() {
        let one = plex(
            "one",
            vec![
                row(1, [Some(4.0), Some(5.0), Some(6.0)]),
                row(2, [Some(1.0), Some(2.0), Some(3.0)]),
                row(3, [None, Some(7.0), Some(7.0)]),
            ],
        );
        let two = plex(
            "two",
            vec![
                row(1, [Some(2.0), Some(2.0), Some(4.0)]),
                row(3, [Some(3.0), Some(4.0), Some(5.0)]),
            ],
        );
        let bridged = bridge(vec![one, two]).unwrap();
        let q = &bridged.quantities;
        let names: Vec<&str> = q.samples.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(names, ["one-a", "one-b", "two-a", "two-b"]);
        assert_eq!(
            bridged.channels,
            ["one 127", "one 128", "two 127", "two 128"]
        );
        assert_eq!(bridged.plex_of, [0, 0, 1, 1]);
        assert_eq!(q.shifts, [0.1, -0.1, 0.1, -0.1]);
        assert_eq!((bridged.psms_used, bridged.psms_incomplete), (20, 2));
        let log2 = |id: i64| {
            q.rows
                .iter()
                .find(|r| r.group_id == id)
                .map(|r| r.log2.clone())
                .unwrap()
        };
        // P1's references read 4 and 2: both move to 3, one by −1, two by +1.
        assert_eq!(log2(1), [Some(4.0), Some(5.0), Some(3.0), Some(5.0)]);
        // P2 only in one: kept as it is.
        assert_eq!(log2(2), [Some(2.0), Some(3.0), None, None]);
        // P3: one has no reference value, so only two's values remain.
        assert_eq!(log2(3), [None, None, Some(4.0), Some(5.0)]);
        assert_eq!(bridged.unbridged, 1);
        assert_eq!(
            q.rows.iter().find(|r| r.group_id == 1).unwrap().peptides,
            [2, 2, 2, 2]
        );
    }

    #[test]
    fn every_plex_needs_a_reference_and_something_else() {
        let mut none = plex("one", Vec::new());
        none.reference = vec![false; 3];
        assert!(matches!(
            bridge(vec![none]),
            Err(QuantError::NoReference(name)) if name == "one"
        ));
        let mut all = plex("one", Vec::new());
        all.reference = vec![true; 3];
        assert!(matches!(
            bridge(vec![all]),
            Err(QuantError::OnlyReference(name)) if name == "one"
        ));
    }
}
