//! `tuffold samples`, `tuffold plex`, and `tuffold quant`.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use tuffold_core::reporter::{Impurities, Plex};
use tuffold_engine::{
    REPORTER_MIN_PURITY, kept_impurities, kept_references, quantify_intensity, quantify_project,
    quantify_reporters_with, write_anova_tsv, write_intensity_tsv, write_quant_tsv,
    write_reporter_tsv,
};
use tuffold_store::{PlexDefinition, Project};

/// List biosamples, after applying `NAME=CATEGORY` assignments (an empty
/// category clears it).
pub fn samples(path: &Path, assignments: &[String]) -> Result<(), String> {
    let mut project = Project::open(path).map_err(|e| e.to_string())?;
    let parsed = assignments
        .iter()
        .map(|assignment| {
            let (name, category) = assignment.split_once('=').ok_or_else(|| {
                format!("'{assignment}': write NAME=CATEGORY (or NAME= to clear)")
            })?;
            Ok((name.trim(), Some(category.trim()).filter(|c| !c.is_empty())))
        })
        .collect::<Result<Vec<_>, String>>()?;
    if !parsed.is_empty() {
        tuffold_engine::set_categories(&mut project, &parsed).map_err(|e| e.to_string())?;
    }
    let list = project.protein_list().map_err(|e| e.to_string())?;
    for biosample in list.biosamples {
        let category = biosample.category.as_deref().unwrap_or("(no category)");
        out!("{}\t{category}", biosample.name);
    }
    Ok(())
}

/// Quantify, optionally compare two categories, and optionally export TSV.
pub fn quant(path: &Path, compare: Option<&[String]>, export: Option<&Path>) -> Result<(), String> {
    let project = Project::open(path).map_err(|e| e.to_string())?;
    let quantities = quantify_project(&project).map_err(|e| e.to_string())?;
    out!(
        "{} protein groups in {} biosamples ({} spectra counted)",
        quantities.rows.len(),
        quantities.samples.len(),
        quantities.totals.iter().sum::<usize>()
    );
    let comparison = match compare {
        Some([a, b]) => Some(quantities.compare(a, b).map_err(|e| e.to_string())?),
        Some(_) => return Err("--compare takes two categories: A B".into()),
        None => None,
    };
    if let Some(c) = &comparison {
        let significant = c.rows.iter().filter(|r| r.q(c.primary) <= 0.05).count();
        out!(
            "{} ({}) vs {} ({}): {significant} protein groups differ at 5% FDR by {}",
            c.b,
            c.samples_b.join(", "),
            c.a,
            c.samples_a.join(", "),
            c.primary.describe()
        );
        for row in c.rows.iter().take(10) {
            out!(
                "  {:<28} log2 fold {:>6.2}  q {:.2e}",
                row.accession,
                row.log2_fold,
                row.q(c.primary)
            );
        }
    }
    if let Some(export) = export {
        let file = File::create(export).map_err(|e| format!("{}: {e}", export.display()))?;
        write_quant_tsv(&mut BufWriter::new(file), &quantities, comparison.as_ref())
            .map_err(|e| format!("{}: {e}", export.display()))?;
        out!("wrote {}", export.display());
    }
    Ok(())
}

/// MS1 intensity quantitation, optionally comparing two categories and
/// exporting TSV.
pub fn intensity(
    path: &Path,
    compare: Option<&[String]>,
    export: Option<&Path>,
) -> Result<(), String> {
    let project = Project::open(path).map_err(|e| e.to_string())?;
    let quant = quantify_intensity(&project).map_err(|e| e.to_string())?;
    for run in &quant.runs {
        match run.ms1_scans {
            Some(scans) if scans > 0 => out!(
                "run {}: {} MS1 scans, {} of {} traced precursors with a peak{}",
                run.run,
                scans,
                run.measured,
                run.traced,
                if run.untraced > 0 {
                    format!(
                        "; {} identified since its spectra file was added (add it again to \
                         trace them)",
                        run.untraced
                    )
                } else {
                    String::new()
                }
            ),
            _ => out!("run {}: no MS1 data; add its mzML to measure it", run.run),
        }
    }
    if quant.transferred > 0 {
        out!(
            "{} run measurements from match-between-runs transfers (tuffold mbr)",
            quant.transferred
        );
    }
    if !quant.incomplete.is_empty() {
        out!(
            "left out of intensities (some of their runs have no MS1 data): {}",
            quant.incomplete.join(", ")
        );
    }
    let q = &quant.quantities;
    let measured = q
        .rows
        .iter()
        .filter(|r| r.log2.iter().any(Option::is_some))
        .count();
    out!(
        "{measured} of {} protein groups have an intensity in {} biosamples (median polish); {}",
        q.rows.len(),
        q.samples.len(),
        if q.normalized.iter().any(|&n| n) {
            format!(
                "normalized on {} peptides measured in two or more biosamples",
                q.shared_peptides
            )
        } else {
            "not normalized (no two biosamples share enough peptides)".to_string()
        }
    );
    let comparison = match compare {
        Some([a, b]) => Some(q.compare(a, b).map_err(|e| e.to_string())?),
        Some(_) => return Err("--compare takes two categories: A B".into()),
        None => None,
    };
    if let Some(c) = &comparison {
        let tested = c.rows.iter().filter(|r| r.p.is_some()).count();
        let significant = c
            .rows
            .iter()
            .filter(|r| r.q.is_some_and(|q| q <= 0.05))
            .count();
        out!(
            "{} ({}) vs {} ({}): {significant} of {tested} tested protein groups differ at 5% \
             FDR by Welch's t-test (groups need two values in each category to be tested)",
            c.b,
            c.samples_b.join(", "),
            c.a,
            c.samples_a.join(", "),
        );
        for row in c.rows.iter().filter(|r| r.q.is_some()).take(10) {
            out!(
                "  {:<28} log2 fold {:>6.2}  q {:.2e}",
                row.accession,
                row.log2_fold.unwrap_or(f64::NAN),
                row.q.unwrap_or(1.0)
            );
        }
    }
    if let Some(export) = export {
        let file = File::create(export).map_err(|e| format!("{}: {e}", export.display()))?;
        write_intensity_tsv(&mut BufWriter::new(file), q, comparison.as_ref())
            .map_err(|e| format!("{}: {e}", export.display()))?;
        out!("wrote {}", export.display());
    }
    Ok(())
}

/// Compare every category at once (one-way ANOVA), by spectral counts, MS1
/// intensity, or reporter ions, optionally exporting the table.
pub fn anova(
    path: &Path,
    by: &str,
    plex: Option<&str>,
    min_purity: Option<f64>,
    export: Option<&Path>,
) -> Result<(), String> {
    let project = Project::open(path).map_err(|e| e.to_string())?;
    let (anova, values) = match by {
        "intensity" => {
            let quant = quantify_intensity(&project).map_err(|e| e.to_string())?;
            (quant.quantities.anova(), "log2 median-polish intensities")
        }
        "reporters" => {
            let quant =
                quantify_reporters_with(&project, plex, min_purity).map_err(|e| e.to_string())?;
            (
                quant.quantities.quantities.anova(),
                "log2 reporter intensities",
            )
        }
        _ => {
            let quantities = quantify_project(&project).map_err(|e| e.to_string())?;
            (quantities.anova(), "log2(normalized counts + 1)")
        }
    };
    let anova = anova.map_err(|e| e.to_string())?;
    let tested = anova.rows.iter().filter(|r| r.p.is_some()).count();
    let significant = anova
        .rows
        .iter()
        .filter(|r| r.q.is_some_and(|q| q <= 0.05))
        .count();
    let groups: Vec<String> = anova
        .categories
        .iter()
        .zip(&anova.members)
        .map(|(c, m)| format!("{c} ({})", m.join(", ")))
        .collect();
    out!("{}", groups.join(" vs "));
    out!(
        "{significant} of {tested} tested protein groups differ across categories at 5% FDR \
         (one-way ANOVA on {values}; {} groups untestable)",
        anova.rows.len() - tested
    );
    for row in anova.rows.iter().filter(|r| r.q.is_some()).take(10) {
        out!(
            "  {:<28} range {:>5.2} log2  F {:>7.2}  q {:.2e}",
            row.accession,
            row.range.unwrap_or(f64::NAN),
            row.f.unwrap_or(f64::NAN),
            row.q.unwrap_or(1.0)
        );
    }
    if let Some(export) = export {
        let file = File::create(export).map_err(|e| format!("{}: {e}", export.display()))?;
        write_anova_tsv(&mut BufWriter::new(file), &anova)
            .map_err(|e| format!("{}: {e}", export.display()))?;
        out!("wrote {}", export.display());
    }
    Ok(())
}

/// Sample QC of one measure: per-sample numbers, PCA, and the correlations
/// least like the rest.
pub fn qc(
    path: &Path,
    by: &str,
    plex: Option<&str>,
    min_purity: Option<f64>,
    export: Option<&Path>,
) -> Result<(), String> {
    let project = Project::open(path).map_err(|e| e.to_string())?;
    let (qc, values) = match by {
        "intensity" => (
            quantify_intensity(&project)
                .map_err(|e| e.to_string())?
                .quantities
                .sample_qc(),
            "log2 median-polish intensities",
        ),
        "reporters" => (
            quantify_reporters_with(&project, plex, min_purity)
                .map_err(|e| e.to_string())?
                .quantities
                .quantities
                .sample_qc(),
            "log2 reporter intensities",
        ),
        _ => (
            quantify_project(&project)
                .map_err(|e| e.to_string())?
                .sample_qc(),
            "log2(normalized counts + 1)",
        ),
    };
    out!(
        "{} samples, {} protein groups ({values})",
        qc.samples.len(),
        qc.groups
    );
    match &qc.pca {
        Some(pca) => out!(
            "PCA on the {} groups measured in every sample: PC1 {:.1}%, PC2 {:.1}%",
            pca.groups,
            pca.explained[0] * 100.0,
            pca.explained[1] * 100.0
        ),
        None => out!("No PCA: it needs three samples and three groups measured in all of them"),
    }
    for (s, sample) in qc.samples.iter().enumerate() {
        let others: Vec<f64> = qc.correlation[s]
            .iter()
            .enumerate()
            .filter(|(other, _)| *other != s)
            .filter_map(|(_, r)| *r)
            .collect();
        let score = qc.pca.as_ref().map(|p| p.scores[s]);
        out!(
            "  {:<24} {:<14} {:>6} groups  median {:>6}  mean r {:>5}  PC1 {:>8}  PC2 {:>8}",
            sample.name,
            sample.category.as_deref().unwrap_or("-"),
            sample.measured,
            sample.median_log2.map_or("-".into(), |v| format!("{v:.2}")),
            if others.is_empty() {
                "-".to_string()
            } else {
                format!("{:.3}", others.iter().sum::<f64>() / others.len() as f64)
            },
            score.map_or("-".into(), |p| format!("{:.2}", p.0)),
            score.map_or("-".into(), |p| format!("{:.2}", p.1)),
        );
    }
    if let Some(export) = export {
        let file = File::create(export).map_err(|e| format!("{}: {e}", export.display()))?;
        tuffold_engine::write_qc_tsv(&mut BufWriter::new(file), &qc)
            .map_err(|e| format!("{}: {e}", export.display()))?;
        out!("wrote {}", export.display());
    }
    Ok(())
}

/// Define a plex (with `--label`), or remove one (with `--remove`), then list
/// the project's plexes.
#[allow(clippy::too_many_arguments)]
pub fn plex(
    path: &Path,
    name: Option<&str>,
    label: Option<&str>,
    remove: bool,
    runs: &[String],
    impurities: Option<&str>,
    references: &[String],
    channels: &[String],
) -> Result<(), String> {
    let mut project = Project::open(path).map_err(|e| e.to_string())?;
    if remove {
        let name = name.ok_or("name the plex to remove: tuffold plex PROJECT NAME --remove")?;
        tuffold_engine::remove_plex(&mut project, name).map_err(|e| e.to_string())?;
    }
    if label.is_none() && (impurities.is_some() || !references.is_empty()) {
        return Err("give the plex's --label with --impurities or --reference".into());
    }
    if let Some(label) = label {
        let name = name.ok_or("name the plex: tuffold plex PROJECT NAME --label ...")?;
        let plex = Plex::parse(label)?;
        let mut assigned = Vec::new();
        for assignment in channels {
            let (channel, sample) = assignment
                .split_once('=')
                .filter(|(c, s)| !c.trim().is_empty() && !s.trim().is_empty())
                .ok_or_else(|| {
                    format!("'{assignment}': write CHANNEL=SAMPLE, e.g. 126=control-1")
                })?;
            assigned.push((channel.trim().to_string(), sample.trim().to_string()));
        }
        if assigned.is_empty() {
            return Err("assign at least one channel: CHANNEL=SAMPLE ...".into());
        }
        let impurities = match impurities {
            None => kept_impurities(&project, name, plex).map_err(|e| e.to_string())?,
            Some("none") => None,
            Some(file) => {
                let text = std::fs::read_to_string(file).map_err(|e| format!("{file}: {e}"))?;
                Some(Impurities::parse(&text, plex).map_err(|e| format!("{file}: {e}"))?)
            }
        };
        let references = match references {
            [] => kept_references(&project, name, plex).map_err(|e| e.to_string())?,
            [none] if none == "none" => Vec::new(),
            given => given.to_vec(),
        };
        let definition = PlexDefinition {
            name: name.to_string(),
            plex,
            runs: runs.to_vec(),
            channels: assigned,
            impurities,
            references,
        };
        tuffold_engine::set_plex(&mut project, &definition).map_err(|e| e.to_string())?;
    }
    let plexes = project.plexes().map_err(|e| e.to_string())?;
    if plexes.is_empty() {
        out!("no plexes; define one with --label");
    }
    for p in plexes {
        out!(
            "{} ({}): runs {}; {}",
            p.name,
            p.plex,
            p.runs.join(", "),
            if p.impurities.is_some() {
                "corrected for the lot's isotope impurities"
            } else {
                "no impurity certificate (reporters as read)"
            }
        );
        for (channel, sample) in &p.channels {
            let reference = if p.references.contains(channel) {
                "\treference"
            } else {
                ""
            };
            out!("  {channel}\t{sample}{reference}");
        }
    }
    Ok(())
}

/// Quantify a plex by reporter ions, optionally compare two categories, and
/// optionally export TSV.
pub fn reporters(
    path: &Path,
    plex: Option<&str>,
    min_purity: Option<f64>,
    compare: Option<&[String]>,
    export: Option<&Path>,
) -> Result<(), String> {
    let project = Project::open(path).map_err(|e| e.to_string())?;
    let quant = quantify_reporters_with(&project, plex, min_purity).map_err(|e| e.to_string())?;
    let r = &quant.quantities;
    let q = &r.quantities;
    let plexes: Vec<String> = quant
        .plexes
        .iter()
        .map(|p| {
            let corrected = if p.impurities.is_some() {
                ", impurity-corrected"
            } else {
                ""
            };
            format!("{} ({}{corrected})", p.name, p.plex)
        })
        .collect();
    out!(
        "{} {}: {} protein groups from {} PSMs with every channel ({} missing a reporter left out{})",
        if plexes.len() == 1 { "plex" } else { "plexes" },
        plexes.join(", "),
        q.rows.len(),
        r.psms_used,
        r.psms_incomplete,
        if quant.without_peaks > 0 {
            format!("; {} passing spectra have no peaks", quant.without_peaks)
        } else {
            String::new()
        }
    );
    match quant.min_purity {
        Some(min) => out!(
            "precursor purity ≥ {min}: {} spectra below it left out, {} without a measured \
             purity kept",
            quant.impure,
            quant.purity_unknown
        ),
        None => out!("no precursor purity cut"),
    }
    if quant.from_ms3 > 0 {
        out!(
            "reporters of {} PSMs read from their SPS-MS3 scans",
            quant.from_ms3
        );
    }
    if quant.plexes.len() > 1 {
        out!(
            "bridged through each plex's reference channels (internal reference scaling), \
             which are left out of the samples; {} protein values left out where a plex's \
             reference had none",
            quant.unbridged
        );
    }
    for (plex, other, complete) in &quant.better_reagents {
        out!(
            "note: in plex {plex}, at the {other} reagents' reporter m/z, {complete} PSMs have \
             every channel; it was probably labelled with those: redefine it with --label {other}"
        );
    }
    for ((sample, channel), shift) in q.samples.iter().zip(&r.channels).zip(&q.shifts) {
        out!(
            "  {channel}\t{}\t{}\tnormalized by {:+.3} log2",
            sample.name,
            sample.category.as_deref().unwrap_or("(no category)"),
            -shift
        );
    }
    let comparison = match compare {
        Some([a, b]) => Some(q.compare(a, b).map_err(|e| e.to_string())?),
        Some(_) => return Err("--compare takes two categories: A B".into()),
        None => None,
    };
    if let Some(c) = &comparison {
        let tested = c.rows.iter().filter(|r| r.p.is_some()).count();
        let significant = c
            .rows
            .iter()
            .filter(|r| r.q.is_some_and(|q| q <= 0.05))
            .count();
        out!(
            "{} vs {}: {significant} of {tested} tested protein groups differ at 5% FDR (Welch's t-test)",
            c.b,
            c.a
        );
        for row in c.rows.iter().filter(|r| r.q.is_some()).take(10) {
            out!(
                "  {:<28} log2 fold {:>6.2}  q {:.2e}",
                row.accession,
                row.log2_fold.unwrap_or(f64::NAN),
                row.q.unwrap_or(1.0)
            );
        }
    }
    if let Some(export) = export {
        let file = File::create(export).map_err(|e| format!("{}: {e}", export.display()))?;
        write_reporter_tsv(&mut BufWriter::new(file), &quant, comparison.as_ref())
            .map_err(|e| format!("{}: {e}", export.display()))?;
        out!("wrote {}", export.display());
    }
    Ok(())
}

/// The purity cut from `--min-purity`: the default when not given, none for
/// `none`, else a fraction from 0 to 1.
pub fn min_purity(given: Option<&str>) -> Result<Option<f64>, String> {
    match given {
        None => Ok(Some(REPORTER_MIN_PURITY)),
        Some("none") => Ok(None),
        Some(text) => text
            .parse::<f64>()
            .ok()
            .filter(|v| (0.0..=1.0).contains(v))
            .map(Some)
            .ok_or_else(|| format!("--min-purity {text}: give a fraction from 0 to 1, or none")),
    }
}
