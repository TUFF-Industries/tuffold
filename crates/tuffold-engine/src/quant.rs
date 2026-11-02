//! Quantitation of a validated project, and its TSV exports: spectral counts
//! ([`quantify_project`]) or MS1 precursor intensity
//! ([`quantify_intensity`]).
//!
//! Both come from the protein groups at the saved thresholds (each group's
//! unique and razor peptides per biosample), so changing a threshold changes
//! what is quantified, as it changes what is listed.

use std::collections::HashMap;
use std::io::{self, Write};

use tuffold_core::reporter::Plex;
use tuffold_quant::anova::AnovaComparison;
use tuffold_quant::intensity::{
    IntensityComparison, IntensityQuantities, IntensityRow, ProteinPeptides, quantify_intensities,
};
use tuffold_quant::qc::SampleQc;
use tuffold_quant::reporter::{
    PlexQuantities, ProteinReporters, ReporterQuantities, bridge,
    quantify_reporters as quantify_plex,
};
use tuffold_quant::{
    Comparison, ProteinCounts, QuantError, QuantRow, Quantities, Sample, quantify,
};
use tuffold_store::{Biosample, Ms1Run, PlexDefinition, Project, StoreError};

#[derive(Debug, thiserror::Error)]
pub enum QuantifyError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error(transparent)]
    Quant(#[from] QuantError),
    #[error("this project has no current validation results; validate it first")]
    NotValidated,
    #[error(
        "no run has MS1 data; add the runs' mzML files (with their MS1 scans) to measure \
         precursor intensity"
    )]
    NoMs1,
    #[error(
        "this project has no plex; define one with `tuffold plex` to quantify by reporter ions"
    )]
    NoPlex,
}

/// Reporter-ion (TMT) quantitation: one plex, or several on one scale.
#[derive(Debug, Clone, PartialEq)]
pub struct ReporterQuant {
    /// The plexes quantified: one, or every plex of the project bridged by
    /// their reference channels (internal reference scaling).
    pub plexes: Vec<PlexDefinition>,
    pub quantities: ReporterQuantities,
    /// Per plex, reagent sets that would have given more PSMs every channel
    /// than the plex's own, with that count: the plex was probably defined
    /// with the wrong one (the two TMT6 sets share channel names).
    pub better_reagents: Vec<(String, Plex, usize)>,
    /// Passing spectra without peaks to read reporters from.
    pub without_peaks: usize,
    /// Passing spectra whose reporters came from their SPS-MS3 scan.
    pub from_ms3: usize,
    /// The precursor purity cut applied, if any; passing spectra left out
    /// below it, and ones whose purity is unknown (kept).
    pub min_purity: Option<f64>,
    pub impure: usize,
    pub purity_unknown: usize,
    /// When bridging: (plex, protein group) values left out because the
    /// plex's reference channels had no value for the group.
    pub unbridged: usize,
}

/// Spectra whose precursor made up less than this share of their isolation
/// window are left out of reporter quantitation (see `tuffold_quant::purity`).
pub const REPORTER_MIN_PURITY: f64 = 0.5;

/// Normalized reporter intensities of every protein group passing the
/// thresholds: in the named plex, or, without a name, in the project's only
/// plex or across all of them, bridged by their reference channels. Spectra
/// below [`REPORTER_MIN_PURITY`] are left out.
pub fn quantify_reporters(
    project: &Project,
    plex: Option<&str>,
) -> Result<ReporterQuant, QuantifyError> {
    quantify_reporters_with(project, plex, Some(REPORTER_MIN_PURITY))
}

/// [`quantify_reporters`] with another purity cut (`None`: every spectrum).
pub fn quantify_reporters_with(
    project: &Project,
    plex: Option<&str>,
    min_purity: Option<f64>,
) -> Result<ReporterQuant, QuantifyError> {
    require_validated(project)?;
    let names: Vec<String> = match plex {
        Some(name) => vec![name.to_string()],
        None => project.plexes()?.into_iter().map(|p| p.name).collect(),
    };
    if names.is_empty() {
        return Err(QuantifyError::NoPlex);
    }
    let mut plexes = Vec::new();
    let mut quantities = Vec::new();
    let mut better_reagents = Vec::new();
    let (mut without_peaks, mut from_ms3, mut impure, mut purity_unknown) = (0, 0, 0, 0);
    for name in &names {
        let one = one_plex(project, name, min_purity)?;
        impure += one.impure;
        purity_unknown += one.purity_unknown;
        better_reagents.extend(one.better.into_iter().map(|(p, n)| (name.clone(), p, n)));
        without_peaks += one.without_peaks;
        from_ms3 += one.from_ms3;
        plexes.push(one.definition);
        quantities.push(one.quantities);
    }
    if let ([quantities], [_]) = (quantities.as_slice(), plexes.as_slice()) {
        return Ok(ReporterQuant {
            plexes,
            quantities: quantities.clone(),
            better_reagents,
            without_peaks,
            from_ms3,
            min_purity,
            impure,
            purity_unknown,
            unbridged: 0,
        });
    }
    let bridged = bridge(
        plexes
            .iter()
            .zip(quantities)
            .map(|(definition, quant)| PlexQuantities {
                name: definition.name.clone(),
                reference: definition
                    .channels
                    .iter()
                    .map(|(c, _)| definition.references.contains(c))
                    .collect(),
                quant,
            })
            .collect(),
    )?;
    Ok(ReporterQuant {
        plexes,
        quantities: ReporterQuantities {
            quantities: bridged.quantities,
            channels: bridged.channels,
            psms_used: bridged.psms_used,
            psms_incomplete: bridged.psms_incomplete,
        },
        better_reagents,
        without_peaks,
        from_ms3,
        min_purity,
        impure,
        purity_unknown,
        unbridged: bridged.unbridged,
    })
}

/// One plex quantified on its own.
struct OnePlex {
    definition: PlexDefinition,
    quantities: ReporterQuantities,
    /// Reagent sets completing more PSMs than the plex's own.
    better: Vec<(Plex, usize)>,
    /// Passing spectra without peaks.
    without_peaks: usize,
    /// Passing spectra read from their SPS-MS3 scan.
    from_ms3: usize,
    impure: usize,
    purity_unknown: usize,
}

fn one_plex(
    project: &Project,
    name: &str,
    min_purity: Option<f64>,
) -> Result<OnePlex, QuantifyError> {
    let list = project.reporter_list(name, min_purity)?;
    let samples = list
        .biosamples
        .iter()
        .map(|b| Sample {
            name: b.name.clone(),
            category: b.category.clone(),
        })
        .collect();
    let channels = list.plex.channels.iter().map(|(c, _)| c.clone()).collect();
    let proteins = list
        .rows
        .into_iter()
        .map(|row| ProteinReporters {
            group_id: row.group_id,
            accession: row.accession,
            psms: row.psms,
        })
        .collect();
    let quantities = quantify_plex(samples, channels, proteins)?;
    let better = list
        .alternatives
        .into_iter()
        .filter(|&(_, complete)| complete > quantities.psms_used)
        .collect();
    Ok(OnePlex {
        definition: list.plex,
        quantities,
        better,
        without_peaks: list.without_peaks,
        from_ms3: list.from_ms3,
        impure: list.impure,
        purity_unknown: list.purity_unknown,
    })
}

/// Tab-separated export of reporter quantitation: per protein group the PSMs
/// summed, then each channel's normalized intensity and its log2, then the
/// comparison's statistics when one is given (rows in its order). Unknown
/// values are empty cells.
pub fn write_reporter_tsv(
    out: &mut impl Write,
    quant: &ReporterQuant,
    comparison: Option<&IntensityComparison>,
) -> io::Result<()> {
    let quantities = &quant.quantities.quantities;
    let labels: Vec<String> = quantities
        .samples
        .iter()
        .zip(&quant.quantities.channels)
        .map(|(s, c)| format!("{} ({c})", clean(&s.name)))
        .collect();
    let mut header = vec!["accession".to_string(), "psms".to_string()];
    for prefix in ["intensity", "log2_normalized"] {
        header.extend(labels.iter().map(|l| format!("{prefix}:{l}")));
    }
    if let Some(c) = comparison {
        header.extend([
            format!("mean_log2:{}", clean(&c.a)),
            format!("mean_log2:{}", clean(&c.b)),
            format!("log2_fold:{}/{}", clean(&c.b), clean(&c.a)),
            "welch_p".into(),
            "welch_q".into(),
        ]);
    }
    writeln!(out, "{}", header.join("\t"))?;
    let plain = |v: Option<f64>| v.map_or_else(String::new, |v| v.to_string());
    let scientific = |v: Option<f64>| v.map_or_else(String::new, |v| format!("{v:e}"));
    let cells = |row: &IntensityRow| {
        let mut cells = vec![
            clean(&row.accession),
            row.peptides
                .first()
                .copied()
                .unwrap_or_default()
                .to_string(),
        ];
        cells.extend(row.intensities.iter().map(|&v| plain(v)));
        cells.extend(row.log2.iter().map(|&v| plain(v)));
        cells
    };
    match comparison {
        None => {
            for row in &quantities.rows {
                writeln!(out, "{}", cells(row).join("\t"))?;
            }
        }
        Some(c) => {
            let by_id: HashMap<i64, &IntensityRow> = quantities
                .rows
                .iter()
                .map(|row| (row.group_id, row))
                .collect();
            for stats in &c.rows {
                let Some(row) = by_id.get(&stats.group_id) else {
                    continue;
                };
                let mut cells = cells(row);
                cells.extend([
                    plain(stats.mean_a),
                    plain(stats.mean_b),
                    plain(stats.log2_fold),
                    scientific(stats.p),
                    scientific(stats.q),
                ]);
                writeln!(out, "{}", cells.join("\t"))?;
            }
        }
    }
    Ok(())
}

/// MS1 intensity quantitation, and how much MS1 data each run had.
#[derive(Debug, Clone, PartialEq)]
pub struct IntensityQuant {
    pub quantities: IntensityQuantities,
    pub runs: Vec<Ms1Run>,
    /// Biosamples left without intensities because some of their runs have
    /// no MS1 data: summing only the traced fractions would make them look
    /// lower.
    pub incomplete: Vec<String>,
    /// Run measurements from match-between-runs transfers.
    pub transferred: usize,
}

fn require_validated(project: &Project) -> Result<(), QuantifyError> {
    let current = project
        .validation_summary()?
        .is_some_and(|summary| summary.validated_at.is_some());
    if current {
        Ok(())
    } else {
        Err(QuantifyError::NotValidated)
    }
}

/// Biosamples with a run lacking MS1 data.
fn incomplete_biosamples(biosamples: &[Biosample], runs: &[Ms1Run]) -> Vec<String> {
    biosamples
        .iter()
        .filter(|b| {
            runs.iter()
                .any(|r| r.biosample == b.name && r.ms1_scans.is_none_or(|n| n == 0))
        })
        .map(|b| b.name.clone())
        .collect()
}

/// Normalized log2 MS1 intensities for every protein group passing the
/// thresholds.
pub fn quantify_intensity(project: &Project) -> Result<IntensityQuant, QuantifyError> {
    require_validated(project)?;
    let list = project.intensity_list()?;
    if !list.has_ms1() {
        return Err(QuantifyError::NoMs1);
    }
    let incomplete = incomplete_biosamples(&list.biosamples, &list.runs);
    let transferred = list.transferred;
    let left_out: Vec<bool> = list
        .biosamples
        .iter()
        .map(|b| incomplete.contains(&b.name))
        .collect();
    let samples = list
        .biosamples
        .into_iter()
        .map(|b| Sample {
            name: b.name,
            category: b.category,
        })
        .collect();
    let proteins = list
        .rows
        .into_iter()
        .map(|row| ProteinPeptides {
            group_id: row.group_id,
            accession: row.accession,
            transferred: row.peptides.iter().map(|p| p.transferred.clone()).collect(),
            peptides: row
                .peptides
                .into_iter()
                .map(|p| {
                    p.intensities
                        .into_iter()
                        .zip(&left_out)
                        .map(|(v, &out)| v.filter(|_| !out))
                        .collect()
                })
                .collect(),
        })
        .collect();
    Ok(IntensityQuant {
        quantities: quantify_intensities(samples, proteins)?,
        runs: list.runs,
        incomplete,
        transferred,
    })
}

/// Normalized counts and NSAF for every protein group passing the thresholds.
pub fn quantify_project(project: &Project) -> Result<Quantities, QuantifyError> {
    require_validated(project)?;
    let list = project.protein_list()?;
    let samples = list
        .biosamples
        .into_iter()
        .map(|b| Sample {
            name: b.name,
            category: b.category,
        })
        .collect();
    let proteins = list
        .rows
        .into_iter()
        .map(|row| ProteinCounts {
            group_id: row.group_id,
            accession: row.accession,
            length: row.length,
            counts: row.spectra_by_sample,
            exclusive_counts: row.exclusive_spectra_by_sample,
            exclusive_peptides: row.exclusive_peptides_by_sample,
            peptides: row.peptides_by_sample,
            observable: row.observable_peptides,
        })
        .collect();
    Ok(quantify(samples, proteins)?)
}

/// Tab-separated export: one row per protein group with counts, normalized
/// counts, and NSAF per biosample, then the comparison's statistics when one
/// is given (rows then follow its order, most significant first; p- and
/// q-values in scientific notation). Unknown values are empty cells, never
/// zeros.
pub fn write_quant_tsv(
    out: &mut impl Write,
    quantities: &Quantities,
    comparison: Option<&Comparison>,
) -> io::Result<()> {
    let mut header = vec!["accession".to_string()];
    for prefix in [
        "count",
        "normalized",
        "nsaf",
        "exclusive",
        "exclusive_peptides",
        "empai",
    ] {
        header.extend(
            quantities
                .samples
                .iter()
                .map(|s| format!("{prefix}:{}", clean(&s.name))),
        );
    }
    if let Some(c) = comparison {
        let pair = format!("{}/{}", clean(&c.b), clean(&c.a));
        header.extend([
            format!("mean_normalized:{}", clean(&c.a)),
            format!("mean_normalized:{}", clean(&c.b)),
            format!("log2_fold:{pair}"),
            "fisher_p".into(),
            "fisher_q".into(),
            "welch_p".into(),
            "welch_q".into(),
            "mann_whitney_p".into(),
            "mann_whitney_q".into(),
        ]);
    }
    writeln!(out, "{}", header.join("\t"))?;
    match comparison {
        None => {
            for row in &quantities.rows {
                writeln!(out, "{}", row_cells(row).join("\t"))?;
            }
        }
        Some(c) => {
            let by_id: HashMap<i64, &QuantRow> = quantities
                .rows
                .iter()
                .map(|row| (row.group_id, row))
                .collect();
            for stats in &c.rows {
                let row = by_id.get(&stats.group_id).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{} is compared but not quantified", stats.accession),
                    )
                })?;
                let optional = |v: Option<f64>| v.map_or_else(String::new, |v| format!("{v:e}"));
                let mut cells = row_cells(row);
                cells.extend([
                    stats.mean_a.to_string(),
                    stats.mean_b.to_string(),
                    stats.log2_fold.to_string(),
                    format!("{:e}", stats.fisher_p),
                    format!("{:e}", stats.fisher_q),
                    optional(stats.welch_p),
                    optional(stats.welch_q),
                    optional(stats.mann_whitney_p),
                    optional(stats.mann_whitney_q),
                ]);
                writeln!(out, "{}", cells.join("\t"))?;
            }
        }
    }
    Ok(())
}

/// Accession, counts, normalized counts, NSAF (empty without a length),
/// exclusive spectra and peptides, and emPAI (empty without a sequence).
fn row_cells(row: &QuantRow) -> Vec<String> {
    let n = row.counts.len();
    let optional = |values: &Option<Vec<f64>>| match values {
        Some(values) => values.iter().map(f64::to_string).collect(),
        None => vec![String::new(); n],
    };
    let counts = |values: &[usize]| -> Vec<String> {
        if values.len() == n {
            values.iter().map(usize::to_string).collect()
        } else {
            vec![String::new(); n]
        }
    };
    let mut cells = vec![clean(&row.accession)];
    cells.extend(row.counts.iter().map(usize::to_string));
    cells.extend(row.normalized.iter().map(f64::to_string));
    cells.extend(optional(&row.nsaf));
    cells.extend(counts(&row.exclusive_counts));
    cells.extend(counts(&row.exclusive_peptides));
    cells.extend(optional(&row.empai));
    cells
}

/// Tab-separated export of MS1 intensities: one row per protein group with
/// median-polish intensities, their log2, and the peptides measured per biosample, then the comparison's statistics when one is
/// given (rows then follow its order). Unknown values are empty cells.
pub fn write_intensity_tsv(
    out: &mut impl Write,
    quantities: &IntensityQuantities,
    comparison: Option<&IntensityComparison>,
) -> io::Result<()> {
    let mut header = vec!["accession".to_string()];
    for prefix in ["intensity", "log2_normalized", "peptides", "top3"] {
        header.extend(
            quantities
                .samples
                .iter()
                .map(|s| format!("{prefix}:{}", clean(&s.name))),
        );
    }
    if let Some(c) = comparison {
        header.extend([
            format!("mean_log2:{}", clean(&c.a)),
            format!("mean_log2:{}", clean(&c.b)),
            format!("measured:{}", clean(&c.a)),
            format!("measured:{}", clean(&c.b)),
            format!("log2_fold:{}/{}", clean(&c.b), clean(&c.a)),
            "welch_p".into(),
            "welch_q".into(),
        ]);
    }
    writeln!(out, "{}", header.join("\t"))?;
    let by_id: HashMap<i64, &IntensityRow> = quantities
        .rows
        .iter()
        .map(|row| (row.group_id, row))
        .collect();
    let plain = |v: Option<f64>| v.map_or_else(String::new, |v| v.to_string());
    let scientific = |v: Option<f64>| v.map_or_else(String::new, |v| format!("{v:e}"));
    let cells = |row: &IntensityRow| {
        let mut cells = vec![clean(&row.accession)];
        cells.extend(row.intensities.iter().map(|&v| plain(v)));
        cells.extend(row.log2.iter().map(|&v| plain(v)));
        cells.extend(row.peptides.iter().map(usize::to_string));
        cells.extend(row.top3.iter().map(|&v| plain(v)));
        cells
    };
    match comparison {
        None => {
            for row in &quantities.rows {
                writeln!(out, "{}", cells(row).join("\t"))?;
            }
        }
        Some(c) => {
            for stats in &c.rows {
                let row = by_id.get(&stats.group_id).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{} is compared but not quantified", stats.accession),
                    )
                })?;
                let mut cells = cells(row);
                cells.extend([
                    plain(stats.mean_a),
                    plain(stats.mean_b),
                    stats.measured_a.to_string(),
                    stats.measured_b.to_string(),
                    plain(stats.log2_fold),
                    scientific(stats.p),
                    scientific(stats.q),
                ]);
                writeln!(out, "{}", cells.join("\t"))?;
            }
        }
    }
    Ok(())
}

/// Tab-separated export of a comparison of every category: per category
/// the mean log2 value and how many samples measured it, then the range
/// (highest − lowest mean), F, p, and q. Unknown values are empty cells.
pub fn write_anova_tsv(out: &mut impl Write, anova: &AnovaComparison) -> io::Result<()> {
    let mut header = vec!["accession".to_string()];
    header.extend(
        anova
            .categories
            .iter()
            .map(|c| format!("mean_log2:{}", clean(c))),
    );
    header.extend(
        anova
            .categories
            .iter()
            .map(|c| format!("measured:{}", clean(c))),
    );
    header.extend(["log2_range", "anova_f", "anova_p", "anova_q"].map(String::from));
    writeln!(out, "{}", header.join("\t"))?;
    let plain = |v: Option<f64>| v.map_or_else(String::new, |v| v.to_string());
    let scientific = |v: Option<f64>| v.map_or_else(String::new, |v| format!("{v:e}"));
    for row in &anova.rows {
        let mut cells = vec![clean(&row.accession)];
        cells.extend(row.means.iter().map(|&m| plain(m)));
        cells.extend(row.measured.iter().map(usize::to_string));
        cells.extend([
            plain(row.range),
            plain(row.f),
            scientific(row.p),
            scientific(row.q),
        ]);
        writeln!(out, "{}", cells.join("\t"))?;
    }
    Ok(())
}

/// Sample QC, one row per sample: groups measured, median, PC1 and PC2 (their
/// share of the variance in the header; empty without a PCA), and Pearson's r
/// with every sample (empty where it couldn't be computed).
pub fn write_qc_tsv(out: &mut impl Write, qc: &SampleQc) -> io::Result<()> {
    let mut header: Vec<String> = ["sample", "category", "measured_groups", "median_log2"]
        .map(String::from)
        .to_vec();
    match &qc.pca {
        Some(pca) => {
            header.push(format!("PC1 ({:.1}%)", pca.explained[0] * 100.0));
            header.push(format!("PC2 ({:.1}%)", pca.explained[1] * 100.0));
        }
        None => header.extend(["PC1", "PC2"].map(String::from)),
    }
    header.extend(qc.samples.iter().map(|s| format!("r:{}", clean(&s.name))));
    writeln!(out, "{}", header.join("\t"))?;
    let number = |v: Option<f64>| v.map_or(String::new(), |v| format!("{v:.4}"));
    for (s, sample) in qc.samples.iter().enumerate() {
        let mut row = vec![
            clean(&sample.name),
            clean(sample.category.as_deref().unwrap_or_default()),
            sample.measured.to_string(),
            number(sample.median_log2),
            number(qc.pca.as_ref().map(|p| p.scores[s].0)),
            number(qc.pca.as_ref().map(|p| p.scores[s].1)),
        ];
        row.extend(qc.correlation[s].iter().map(|&r| number(r)));
        writeln!(out, "{}", row.join("\t"))?;
    }
    out.flush()
}

/// Tabs and line breaks inside a field would break the table.
fn clean(field: &str) -> String {
    field.replace(['\t', '\n', '\r'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_biosample_with_an_untraced_run_is_incomplete() {
        let biosample = |id, name: &str| Biosample {
            id,
            name: name.into(),
            category: None,
        };
        let run = |run: &str, biosample: &str, ms1_scans| Ms1Run {
            run: run.into(),
            biosample: biosample.into(),
            ms1_scans,
            traced: 0,
            measured: 0,
            untraced: 0,
        };
        // b has two fractions and only one was traced; c has none at all.
        let incomplete = incomplete_biosamples(
            &[biosample(1, "a"), biosample(2, "b"), biosample(3, "c")],
            &[
                run("a1", "a", Some(900)),
                run("b1", "b", Some(800)),
                run("b2", "b", None),
                run("c1", "c", Some(0)),
            ],
        );
        assert_eq!(incomplete, ["b", "c"]);
    }
}

#[cfg(test)]
mod reporter_tests {
    //! TMT reporter quantitation end to end: the demo's spectra get TMT6
    //! reporter ions (one protein tripled in three channels, channel 131
    //! loaded twice), then a plex is defined, quantified, and compared.

    use std::collections::{HashMap, HashSet};

    use tuffold_core::reporter::{Impurities, Plex};
    use tuffold_core::{PeakSpectrum, SourceFormat};
    use tuffold_store::{ImportSource, PlexDefinition, Project};
    use tuffold_validate::Settings;

    use crate::demo::generate;
    use crate::import::attach_scans;
    use crate::{quantify_reporters, set_categories, set_plex, validate_project};

    /// The demo project, validated, with TMT6 reporters on every MS/MS
    /// spectrum: `reporters(run, scan index, up)` gives each channel's
    /// intensity, `up` for the spectra of the changed protein (the one with
    /// the most rank-1 spectra). With `ms3` the reporters go into an MS3
    /// scan after each MS/MS one, naming it as parent (SPS-MS3), and the
    /// MS/MS spectra get none. Returns the project, its runs, and that
    /// protein's accession.
    fn tmt_project(
        path: &std::path::Path,
        ms3: bool,
        reporters: impl Fn(&str, usize, bool) -> Vec<f64>,
    ) -> (Project, Vec<String>, String) {
        labelled_project(path, Plex::Tmt6, 0.0, ms3, reporters)
    }

    /// [`tmt_project`] for any reagent set, with its reporters `ppm_off` their
    /// channels' m/z.
    fn labelled_project(
        path: &std::path::Path,
        plex: Plex,
        ppm_off: f64,
        ms3: bool,
        reporters: impl Fn(&str, usize, bool) -> Vec<f64>,
    ) -> (Project, Vec<String>, String) {
        let demo = generate(7);
        let mut spectra_of: HashMap<&str, HashSet<(String, tuffold_core::SpectrumKey)>> =
            HashMap::new();
        for psm in demo.batch.psms.iter().filter(|p| p.rank == 1) {
            for evidence in &psm.evidence {
                spectra_of
                    .entry(evidence.accession.as_str())
                    .or_default()
                    .insert((psm.ms_run.clone(), psm.spectrum.clone()));
            }
        }
        let (changed, changed_spectra) = spectra_of
            .iter()
            .filter(|(accession, _)| {
                !accession.starts_with("DECOY") && !accession.starts_with("rev")
            })
            .max_by_key(|(accession, spectra)| (spectra.len(), *accession))
            .map(|(a, s)| (a.to_string(), s.clone()))
            .unwrap();
        let mut project = Project::create(path).unwrap();
        let source = ImportSource {
            path: "TMT test (synthetic)".into(),
            format: SourceFormat::Synthetic,
        };
        project.import(&source, &demo.batch).unwrap();
        let channels = plex.channels();
        for (run, scans) in demo.peaks.clone() {
            let scans: Vec<PeakSpectrum> = scans
                .into_iter()
                .enumerate()
                .flat_map(|(i, mut scan)| {
                    if scan.ms_level < 2 {
                        return vec![scan];
                    }
                    let up = changed_spectra.contains(&(run.clone(), scan.key.clone()));
                    let mut with_reporters = if ms3 {
                        PeakSpectrum {
                            key: tuffold_core::SpectrumKey::Native(format!("ms3-{i}")),
                            native_id: format!("ms3-{i}"),
                            ms_level: 3,
                            mz: Vec::new(),
                            intensity: Vec::new(),
                            parent: Some(scan.native_id.clone()),
                            ..scan.clone()
                        }
                    } else {
                        scan.clone()
                    };
                    for (channel, value) in channels.iter().zip(reporters(&run, i, up)) {
                        with_reporters.mz.push(channel.mz * (1.0 + ppm_off * 1e-6));
                        with_reporters.intensity.push(value as f32);
                    }
                    if ms3 {
                        scan.parent = None;
                        vec![scan, with_reporters]
                    } else {
                        vec![with_reporters]
                    }
                })
                .collect();
            attach_scans(&mut project, &run, scans.into_iter().map(Ok)).unwrap();
        }
        validate_project(&mut project, &Settings::default()).unwrap();
        let runs = demo.batch.ms_runs.iter().map(|r| r.name.clone()).collect();
        (project, runs, changed)
    }

    /// A few percent of deterministic noise per scan and channel, so tests
    /// have variance.
    fn noise(scan: usize, channel: usize) -> f64 {
        1.0 + 0.03 * ((scan * 7 + channel * 13) as f64).sin()
    }

    /// An 8-plex end to end, reporters 50 ppm off (a TOF's, or a neighbour's
    /// isotope at this nominal mass): the protein that tripled in the last
    /// four channels comes back tripled once the certificate is applied, and
    /// the spill past 120's empty place goes nowhere.
    #[test]
    fn an_itraq_plex_recovers_its_ratios_through_its_certificate() {
        let dir = tempfile::tempdir().unwrap();
        let plex = Plex::Itraq8;
        let spill = 0.10;
        let (mut project, runs, changed) = labelled_project(
            &dir.path().join("itraq.tuff"),
            plex,
            50.0,
            false,
            |_, i, up| {
                let truth: Vec<f64> = (0..8)
                    .map(|c| {
                        let change = if up && c >= 4 { 3.0 } else { 1.0 };
                        10_000.0 * change * noise(i, c)
                    })
                    .collect();
                // 116 (index 3) loses 10% of its reporters to 117.
                let mut observed = truth.clone();
                observed[3] -= spill * truth[3];
                observed[4] += spill * truth[3];
                observed
            },
        );
        let samples = ["c1", "c2", "c3", "c4", "t1", "t2", "t3", "t4"];
        let certificate = "113 0/0/0/0\n114 0/0/0/0\n115 0/0/0/0\n116 0/0/10/0\n\
                           117 0/0/0/0\n118 0/0/0/0\n119 0/0/0/0\n121 0/0/0/0\n";
        let definition = PlexDefinition {
            name: "itraq".into(),
            plex,
            runs,
            channels: plex
                .channels()
                .iter()
                .zip(samples)
                .map(|(c, s)| (c.name.clone(), s.to_string()))
                .collect(),
            impurities: Some(Impurities::parse(certificate, plex).unwrap()),
            references: Vec::new(),
        };
        set_plex(&mut project, &definition).unwrap();
        assert_eq!(project.plexes().unwrap(), std::slice::from_ref(&definition));
        let groups: Vec<(&str, Option<&str>)> = samples
            .iter()
            .map(|s| {
                let side = if s.starts_with('c') {
                    "control"
                } else {
                    "treated"
                };
                (*s, Some(side))
            })
            .collect();
        set_categories(&mut project, &groups).unwrap();

        let quant = quantify_reporters(&project, None).unwrap();
        let r = &quant.quantities;
        assert_eq!(r.psms_incomplete, 0, "60 ppm reads what is 50 off");
        assert!(r.psms_used > 100, "{}", r.psms_used);
        let comparison = r.quantities.compare("control", "treated").unwrap();
        let up = comparison
            .rows
            .iter()
            .find(|row| row.accession.contains(&changed))
            .expect("the changed protein is quantified");
        assert!(
            (up.log2_fold.unwrap() - 3f64.log2()).abs() < 0.05,
            "{:?}",
            up.log2_fold
        );
    }

    #[test]
    fn a_plex_recovers_its_ratios_after_loading_correction() {
        let dir = tempfile::tempdir().unwrap();
        let channels = Plex::Tmt6.channels();
        let (mut project, runs, changed) =
            tmt_project(&dir.path().join("tmt.tuff"), false, |_, i, up| {
                let truth: Vec<f64> = (0..6)
                    .map(|c| {
                        let change = if up && c >= 3 { 3.0 } else { 1.0 };
                        // Channel 131 was loaded twice.
                        let loading = if c == 5 { 2.0 } else { 1.0 };
                        10_000.0 * change * loading * noise(i, c)
                    })
                    .collect();
                // The 129 reagent's impurity: 20% of its reporters two ¹³C
                // lighter, i.e. at 127 (both carry a ¹⁵N).
                let mut observed = truth.clone();
                observed[1] += 0.2 * truth[3];
                observed[3] *= 0.8;
                observed
            });
        let samples = ["c1", "c2", "c3", "t1", "t2", "t3"];
        let definition = PlexDefinition {
            name: "pool".into(),
            plex: Plex::Tmt6,
            runs,
            channels: channels
                .iter()
                .zip(samples)
                .map(|(c, s)| (c.name.clone(), s.to_string()))
                .collect(),
            impurities: None,
            references: Vec::new(),
        };
        set_plex(&mut project, &definition).unwrap();
        assert_eq!(project.plexes().unwrap(), std::slice::from_ref(&definition));
        let groups: Vec<(&str, Option<&str>)> = samples
            .iter()
            .map(|s| {
                (
                    *s,
                    Some(if s.starts_with('c') {
                        "control"
                    } else {
                        "treated"
                    }),
                )
            })
            .collect();
        set_categories(&mut project, &groups).unwrap();

        let quant = quantify_reporters(&project, None).unwrap();
        let r = &quant.quantities;
        assert_eq!(r.psms_incomplete, 0);
        assert!(r.psms_used > 100, "{}", r.psms_used);
        // 131's double loading is normalized away: its shift is log2(2) above the
        // others' (shifts are relative to each protein's mean over channels).
        let shifts = &r.quantities.shifts;
        assert!((shifts[5] - shifts[0] - 1.0).abs() < 0.05, "{shifts:?}");
        let comparison = r.quantities.compare("control", "treated").unwrap();
        let row = |accession: &str| {
            comparison
                .rows
                .iter()
                .find(|row| row.accession.contains(accession))
        };
        // 129's spill into 127 lifts the control side of the changed protein
        // only: 127 reads (1 + 0.2 · 3)/(1 + 0.2) = 4/3 of the others, so the
        // fold loses log2(4/3)/3 = 0.14.
        let up = row(&changed).expect("the changed protein is quantified");
        let compressed = 3f64.log2() - up.log2_fold.unwrap();
        assert!(
            (compressed - (4f64 / 3.0).log2() / 3.0).abs() < 0.05,
            "{compressed}"
        );

        // With the lot's certificate the spill is undone.
        let certificate = "126 0 0 0 0\n127 0 0 0 0\n128 0 0 0 0\n129 20 0 0 0\n130 0 0 0 0\n\
                           131 0 0 0 0\n";
        let corrected = PlexDefinition {
            impurities: Some(Impurities::parse(certificate, Plex::Tmt6).unwrap()),
            ..definition.clone()
        };
        set_plex(&mut project, &corrected).unwrap();
        let quant = quantify_reporters(&project, None).unwrap();
        assert!(quant.plexes[0].impurities.is_some());
        let comparison = quant
            .quantities
            .quantities
            .compare("control", "treated")
            .unwrap();
        let row = |accession: &str| {
            comparison
                .rows
                .iter()
                .find(|row| row.accession.contains(accession))
        };
        let up = row(&changed).expect("the changed protein is quantified");
        assert!(
            (up.log2_fold.unwrap() - 3f64.log2()).abs() < 0.02,
            "{changed}: {:?}",
            up.log2_fold
        );
        assert!(up.q.unwrap() < 0.01, "{:?}", up.q);

        // Spectra measured as impure are left out at the default cut and
        // kept without one; unmeasured ones are kept and counted.
        let all = quant.quantities.psms_used;
        rusqlite::Connection::open(project.path().unwrap())
            .unwrap()
            .execute(
                "UPDATE spectra SET precursor_purity = CASE WHEN id % 2 = 0 THEN 0.2 ELSE 0.9 END
                 WHERE id % 3 = 0",
                [],
            )
            .unwrap();
        let filtered = quantify_reporters(&project, None).unwrap();
        assert!(filtered.impure > 0);
        assert_eq!(filtered.quantities.psms_used + filtered.impure, all);
        assert!(filtered.purity_unknown > filtered.impure);
        let unfiltered = crate::quantify_reporters_with(&project, None, None).unwrap();
        assert_eq!(
            (unfiltered.quantities.psms_used, unfiltered.impure),
            (all, 0)
        );
        let unchanged = comparison
            .rows
            .iter()
            .filter(|row| !row.accession.contains(&changed))
            .filter_map(|row| row.log2_fold)
            .map(f64::abs)
            .fold(0.0, f64::max);
        assert!(unchanged < 0.1, "largest unchanged fold {unchanged}");
    }

    #[test]
    fn two_plexes_meet_through_their_reference_channel() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("two.tuff");
        // Runs 3 and 4 are a second plex whose spectra were each picked at
        // their own intensity (up to ±1.5 log2): within-plex normalization
        // can't undo that per protein, the reference channel can. 126 holds
        // the pool of the five samples (the changed protein at 2.2×).
        let second: Vec<String> = generate(7).batch.ms_runs[2..]
            .iter()
            .map(|r| r.name.clone())
            .collect();
        let (mut project, runs, changed) = tmt_project(&path, false, |run, i, up| {
            let scale = if second.iter().any(|r| r == run) {
                (1.5 * (i as f64 * 0.37).sin()).exp2()
            } else {
                1.0
            };
            (0..6)
                .map(|c| {
                    let change = match (up, c) {
                        (false, _) | (true, 1 | 2) => 1.0,
                        (true, 0) => 2.2,
                        (true, _) => 3.0,
                    };
                    10_000.0 * change * scale * noise(i, c)
                })
                .collect()
        });
        let roles = ["pool", "c1", "c2", "t1", "t2", "t3"];
        for (name, runs) in [("a", &runs[..2]), ("b", &runs[2..])] {
            let definition = PlexDefinition {
                name: name.into(),
                plex: Plex::Tmt6,
                runs: runs.to_vec(),
                channels: Plex::Tmt6
                    .channels()
                    .iter()
                    .zip(roles)
                    .map(|(c, role)| (c.name.clone(), format!("{name}-{role}")))
                    .collect(),
                impurities: None,
                references: vec!["126".into()],
            };
            set_plex(&mut project, &definition).unwrap();
        }
        let groups: Vec<(String, Option<&str>)> = ["a", "b"]
            .iter()
            .flat_map(|plex| {
                roles[1..].iter().map(move |role| {
                    let category = if role.starts_with('c') {
                        "control"
                    } else {
                        "treated"
                    };
                    (format!("{plex}-{role}"), Some(category))
                })
            })
            .collect();
        let groups: Vec<(&str, Option<&str>)> =
            groups.iter().map(|(n, c)| (n.as_str(), *c)).collect();
        set_categories(&mut project, &groups).unwrap();

        let quant = quantify_reporters(&project, None).unwrap();
        assert_eq!(quant.plexes.len(), 2);
        let q = &quant.quantities.quantities;
        let names: Vec<&str> = q.samples.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(
            names,
            [
                "a-c1", "a-c2", "a-t1", "a-t2", "a-t3", "b-c1", "b-c2", "b-t1", "b-t2", "b-t3"
            ]
        );
        assert_eq!(quant.quantities.channels[5], "b 127");
        // On one scale, a protein reads the same in both plexes' c1.
        let mut gaps: Vec<f64> = q
            .rows
            .iter()
            .filter_map(|r| Some((r.log2[0]? - r.log2[5]?).abs()))
            .collect();
        gaps.sort_by(f64::total_cmp);
        assert!(gaps.len() > 50, "{}", gaps.len());
        assert!(
            gaps[gaps.len() / 2] < 0.05,
            "median gap {}",
            gaps[gaps.len() / 2]
        );
        let comparison = q.compare("control", "treated").unwrap();
        let fold = |accession: &str| {
            comparison
                .rows
                .iter()
                .find(|r| r.accession.contains(accession))
                .and_then(|r| r.log2_fold)
        };
        let up = fold(&changed).unwrap();
        assert!((up - 3f64.log2()).abs() < 0.1, "{up}");
        let unchanged = comparison
            .rows
            .iter()
            .filter(|r| !r.accession.contains(&changed))
            .filter_map(|r| r.log2_fold)
            .map(f64::abs)
            .fold(0.0, f64::max);
        assert!(unchanged < 0.1, "largest unchanged fold {unchanged}");

        // Without a reference in b there is nothing to bridge with.
        let b = project
            .plexes()
            .unwrap()
            .into_iter()
            .find(|p| p.name == "b")
            .unwrap();
        set_plex(
            &mut project,
            &PlexDefinition {
                references: Vec::new(),
                ..b
            },
        )
        .unwrap();
        let error = quantify_reporters(&project, None).unwrap_err().to_string();
        assert!(
            error.starts_with("plex b has no reference channel"),
            "{error}"
        );
        // Each plex can still be quantified on its own.
        assert_eq!(
            quantify_reporters(&project, Some("b"))
                .unwrap()
                .quantities
                .channels
                .len(),
            6
        );
    }

    #[test]
    fn sps_ms3_reporters_are_read_from_the_ms3_scan() {
        let dir = tempfile::tempdir().unwrap();
        let (mut project, runs, changed) =
            tmt_project(&dir.path().join("ms3.tuff"), true, |_, i, up| {
                (0..6)
                    .map(|c| {
                        let change = if up && c >= 3 { 3.0 } else { 1.0 };
                        10_000.0 * change * noise(i, c)
                    })
                    .collect()
            });
        let samples = ["c1", "c2", "c3", "t1", "t2", "t3"];
        let definition = PlexDefinition {
            name: "pool".into(),
            plex: Plex::Tmt6,
            runs,
            channels: Plex::Tmt6
                .channels()
                .iter()
                .zip(samples)
                .map(|(c, s)| (c.name.clone(), s.to_string()))
                .collect(),
            impurities: None,
            references: Vec::new(),
        };
        set_plex(&mut project, &definition).unwrap();
        let groups: Vec<(&str, Option<&str>)> = samples
            .iter()
            .map(|s| {
                (
                    *s,
                    Some(if s.starts_with('c') {
                        "control"
                    } else {
                        "treated"
                    }),
                )
            })
            .collect();
        set_categories(&mut project, &groups).unwrap();

        let quant = quantify_reporters(&project, None).unwrap();
        let r = &quant.quantities;
        // Every passing spectrum has its MS3, and the MS/MS ones hold no
        // reporters, so all used PSMs came from MS3 scans.
        assert!(r.psms_used > 100, "{}", r.psms_used);
        assert_eq!(r.psms_incomplete, 0);
        assert_eq!(quant.from_ms3, r.psms_used);
        let comparison = r.quantities.compare("control", "treated").unwrap();
        let up = comparison
            .rows
            .iter()
            .find(|row| row.accession.contains(&changed))
            .and_then(|row| row.log2_fold)
            .unwrap();
        assert!((up - 3f64.log2()).abs() < 0.1, "{up}");
    }
}
