//! What each MCP tool does, as plain functions over the engine crates that
//! return JSON. No MCP types here, so every tool is tested directly.
//!
//! Read tools open the project read-only in effect (opening never writes);
//! write tools are the same workflows the CLI and app run. Lists are paged
//! and always say how many rows exist in all, so nothing is dropped
//! silently.

use std::path::{Path, PathBuf};

use serde_json::{Value, json};
use tuffold_core::coverage::covered_fraction;
use tuffold_engine::{
    OutcomeKind, SearchSettings, apply_thresholds, create_demo_project, duplicate_note,
    import_inputs, quantify_project, search_files, validate_project,
};
use tuffold_quant::PrimaryTest;
use tuffold_store::{LevelCounts, Project, ValidationSummary};
use tuffold_validate::{Settings, ThresholdRequest, Thresholds};

/// The most rows one call returns; ask again with an offset for more.
pub const MAX_ROWS: usize = 500;

type Result<T> = std::result::Result<T, String>;

/// Project paths must be absolute: a relative one would resolve against the
/// server's working directory, which the client does not control.
fn absolute(path: &Path) -> Result<&Path> {
    if path.is_absolute() {
        Ok(path)
    } else {
        Err(format!(
            "{}: give the project as an absolute path",
            path.display()
        ))
    }
}

fn open(path: &Path) -> Result<Project> {
    Project::open(absolute(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

/// Run a write workflow on a project, creating it when missing. `work`
/// returns its report and whether anything was imported; when nothing was,
/// the call is an error, and a project this call created is removed again
/// rather than left behind empty.
fn importing(
    path: &Path,
    work: impl FnOnce(&mut Project) -> Result<(Value, bool)>,
) -> Result<Value> {
    let path = absolute(path)?;
    let created = !path.exists();
    let mut project = if created {
        Project::create(path)
    } else {
        Project::open(path)
    }
    .map_err(|e| format!("{}: {e}", path.display()))?;
    let outcome = work(&mut project);
    drop(project);
    let imported = matches!(outcome, Ok((_, true)));
    if created && !imported {
        for suffix in ["", "-wal", "-shm", "-journal"] {
            let mut file = path.as_os_str().to_owned();
            file.push(suffix);
            // Best effort: the files may not exist.
            let _ = std::fs::remove_file(file);
        }
    }
    match outcome {
        Ok((report, true)) => Ok(report),
        Ok((report, false)) => Err(format!("nothing was imported: {report}")),
        Err(message) => Err(message),
    }
}

/// The window `[start, end)` of `total` rows, and how it was served.
fn page(offset: usize, limit: usize, total: usize) -> (usize, usize, Value) {
    let limit = limit.clamp(1, MAX_ROWS);
    let start = offset.min(total);
    let end = (start + limit).min(total);
    let served = json!({ "offset": start, "limit": limit, "returned": end - start });
    (start, end, served)
}

fn counts(level: &LevelCounts) -> Value {
    json!({
        "targets": level.targets, "decoys": level.decoys, "lowest_q": level.lowest_q,
        "estimated_fdr": level.estimated_fdr(),
    })
}

fn validation_json(summary: &ValidationSummary) -> Value {
    let t = &summary.settings.thresholds;
    json!({
        "current": summary.validated_at.is_some(),
        "rescore": summary.settings.rescore,
        "thresholds": {
            "by": t.mode.as_str(),
            "psm": t.psm, "peptide": t.peptide, "protein": t.protein,
            "min_peptides": t.min_peptides,
        },
        "passing": {
            "spectra": counts(&summary.spectra),
            "peptides": counts(&summary.peptides),
            "protein_groups": counts(&summary.protein_groups),
        },
        "engines": summary.engines.iter().map(|e| json!({
            "name": e.name, "psms": e.psms, "decoys": e.decoys,
            "method": e.method, "detail": e.detail,
        })).collect::<Vec<_>>(),
        "peptide_scoring": summary.peptide_scoring.as_ref().map(|p| json!({
            "method": p.method, "detail": p.detail,
        })),
        "localization": {
            "probability_threshold": summary.localization.probability_threshold,
            "scored_psms": summary.localization.scored_psms,
            "localized_psms": summary.localization.localized_psms,
            "total_target_sites": summary.localization.total_target_sites,
            "localized_sites": summary.localization.localized_sites,
            "estimated_flr": summary.localization.estimated_flr,
            "missing_peaks_psms": summary.localization.missing_peaks_psms,
            "capped_psms": summary.localization.capped_psms,
            "unscorable_psms": summary.localization.unscorable_psms,
        },
        "notes": summary.notes.iter().map(|n| json!({
            "code": n.code.as_str(), "message": n.message,
        })).collect::<Vec<_>>(),
    })
}

fn plex_json(p: &tuffold_store::PlexDefinition) -> Value {
    json!({
        "name": p.name,
        "label": p.plex.name(),
        "runs": p.runs,
        "impurity_corrected": p.impurities.is_some(),
        "references": p.references,
        "channels": p.channels.iter().map(|(c, s)| json!({
            "channel": c, "biosample": s,
        })).collect::<Vec<_>>(),
    })
}

/// Contents, imports, and validation state of a project.
pub fn project_summary(path: &Path) -> Result<Value> {
    let project = open(path)?;
    let s = project.summary().map_err(|e| e.to_string())?;
    let imports = project.imports().map_err(|e| e.to_string())?;
    let validation = project.validation_summary().map_err(|e| e.to_string())?;
    let plexes = project.plexes().map_err(|e| e.to_string())?;
    let biosamples = project
        .protein_list()
        .map_err(|e| e.to_string())?
        .biosamples;
    Ok(json!({
        "project": path.display().to_string(),
        "counts": {
            "imports": s.imports, "biosamples": s.biosamples, "ms_runs": s.ms_runs,
            "spectra": s.spectra, "psms": s.psms, "decoy_psms": s.decoy_psms,
            "peptides": s.peptides, "proteins": s.proteins,
            "decoy_proteins": s.decoy_proteins,
            "proteins_with_sequence": s.proteins_with_sequence,
        },
        "runs": s.runs.iter().map(|r| json!({
            "name": r.name, "biosample": r.biosample, "spectra": r.spectra, "psms": r.psms,
        })).collect::<Vec<_>>(),
        "biosamples": biosamples.iter().map(|b| json!({
            "name": b.name, "category": b.category,
        })).collect::<Vec<_>>(),
        "plexes": plexes.iter().map(plex_json).collect::<Vec<_>>(),
        "imports": imports.iter().map(|i| json!({
            "source": i.source_path, "format": i.format, "software": i.software,
            "version": i.software_version, "psms": i.psms,
        })).collect::<Vec<_>>(),
        "validation": validation.as_ref().map(validation_json),
    }))
}

/// Protein groups passing the thresholds, best first; `filter` keeps those
/// whose accession (or, with `match_descriptions`, description) contains it
/// (case-insensitive).
pub fn list_protein_groups(
    path: &Path,
    filter: Option<&str>,
    match_descriptions: bool,
    offset: usize,
    limit: usize,
) -> Result<Value> {
    let list = open(path)?.protein_list().map_err(|e| e.to_string())?;
    let needle = filter.map(str::to_lowercase).filter(|n| !n.is_empty());
    let matching: Vec<_> = list
        .rows
        .iter()
        .filter(|row| {
            needle.as_ref().is_none_or(|n| {
                row.accession.to_lowercase().contains(n)
                    || match_descriptions
                        && row
                            .description
                            .as_ref()
                            .is_some_and(|d| d.to_lowercase().contains(n))
            })
        })
        .collect();
    let (start, end, served) = page(offset, limit, matching.len());
    let samples: Vec<&str> = list.biosamples.iter().map(|b| b.name.as_str()).collect();
    Ok(json!({
        "total": matching.len(),
        "page": served,
        "biosamples": samples,
        "groups": matching[start..end].iter().map(|row| json!({
            "group_id": row.group_id, "accession": row.accession,
            "description": row.description, "equivalents": row.equivalents,
            "peptides": row.peptides, "spectra": row.spectra,
            "spectra_by_biosample": row.spectra_by_sample,
            "exclusive_spectra_by_biosample": row.exclusive_spectra_by_sample,
            "exclusive_peptides_by_biosample": row.exclusive_peptides_by_sample,
            "go_slim": row.go,
            "q_value": row.q_value, "probability": row.probability,
            "coverage": row.coverage,
        })).collect::<Vec<_>>(),
    }))
}

/// One protein group: members, peptides (with roles), and coverage.
pub fn protein_group(path: &Path, group_id: i64) -> Result<Value> {
    let detail = open(path)?
        .group_detail(group_id)
        .map_err(|e| e.to_string())?
        .ok_or_else(|| format!("no current protein group has id {group_id}"))?;
    let length = detail.sequence.as_ref().map(|s| s.chars().count());
    Ok(json!({
        "group_id": detail.group_id,
        "members": detail.members.iter().map(|m| json!({
            "accession": m.accession, "description": m.description, "role": m.role,
        })).collect::<Vec<_>>(),
        "sequence_length": length,
        "coverage": length.map(|length| covered_fraction(length, &detail.coverage)),
        "peptides": detail.peptides.iter().map(|p| json!({
            "peptide_id": p.peptide_id, "proforma": p.proforma, "role": p.role,
            "q_value": p.q_value, "probability": p.probability, "spectra": p.spectra,
            "localized_peptide": p.localized_peptide,
            "site_probabilities": p.site_probabilities,
            "best_ascore": p.best_ascore,
            "min_site_probability": p.min_site_probability,
        })).collect::<Vec<_>>(),
    }))
}

/// A peptide's best spectra (at most [`MAX_ROWS`]; `total` says how many
/// it has).
pub fn peptide_spectra(path: &Path, peptide_id: i64, limit: usize) -> Result<Value> {
    let limit = limit.clamp(1, MAX_ROWS);
    let page = open(path)?
        .peptide_spectra(peptide_id, limit)
        .map_err(|e| e.to_string())?;
    Ok(json!({
        "total": page.total,
        "page": { "offset": 0, "limit": limit, "returned": page.rows.len() },
        "spectra": page.rows.iter().map(|s| json!({
            "spectrum_id": s.spectrum_id, "run": s.run, "biosample": s.biosample,
            "spectrum": s.key, "charge": s.charge, "precursor_mz": s.precursor_mz,
            "rt_seconds": s.rt_seconds, "q_value": s.q_value,
            "probability": s.probability, "engines": s.engines,
            "passes": s.passes, "has_peaks": s.has_peaks,
            "localized_peptide": s.localized_peptide,
            "site_probabilities": s.site_probabilities,
            "best_ascore": s.best_ascore,
            "min_site_probability": s.min_site_probability,
        })).collect::<Vec<_>>(),
    }))
}

/// Compare two biosample categories (versus against baseline), most
/// significant first by the primary test.
pub fn compare_categories(
    path: &Path,
    baseline: &str,
    versus: &str,
    offset: usize,
    limit: usize,
) -> Result<Value> {
    let quantities = quantify_project(&open(path)?).map_err(|e| e.to_string())?;
    let comparison = quantities
        .compare(baseline, versus)
        .map_err(|e| e.to_string())?;
    let primary = comparison.primary;
    let q = |row: &tuffold_quant::ComparisonRow| match primary {
        PrimaryTest::Welch => row.welch_q,
        PrimaryTest::Fisher => Some(row.fisher_q),
        PrimaryTest::MannWhitney => row.mann_whitney_q,
    };
    let mut rows: Vec<_> = comparison.rows.iter().collect();
    // Untestable rows (q unknown) go last.
    rows.sort_by(|a, b| match (q(a), q(b)) {
        (Some(x), Some(y)) => x.total_cmp(&y),
        (x, y) => x.is_none().cmp(&y.is_none()),
    });
    let (start, end, served) = page(offset, limit, rows.len());
    Ok(json!({
        "baseline": comparison.a, "versus": comparison.b,
        "baseline_samples": comparison.samples_a, "versus_samples": comparison.samples_b,
        "primary_test": primary.describe(),
        "total": rows.len(),
        "page": served,
        "significant_at_5pct": rows.iter().filter(|r| q(r).is_some_and(|q| q <= 0.05)).count(),
        "rows": rows[start..end].iter().map(|r| json!({
            "group_id": r.group_id, "accession": r.accession,
            "mean_baseline": r.mean_a, "mean_versus": r.mean_b,
            "log2_fold": r.log2_fold, "q_value": q(r),
            "fisher_p": r.fisher_p, "welch_p": r.welch_p,
            "mann_whitney_p": r.mann_whitney_p,
        })).collect::<Vec<_>>(),
    }))
}

/// Compare two categories by MS1 intensity (median polish, Welch's t-test), most
/// significant first; groups that couldn't be tested last.
pub fn compare_intensities(
    path: &Path,
    baseline: &str,
    versus: &str,
    offset: usize,
    limit: usize,
) -> Result<Value> {
    let quant = tuffold_engine::quantify_intensity(&open(path)?).map_err(|e| e.to_string())?;
    let comparison = quant
        .quantities
        .compare(baseline, versus)
        .map_err(|e| e.to_string())?;
    let rows = &comparison.rows;
    let (start, end, served) = page(offset, limit, rows.len());
    Ok(json!({
        "method": "ms1_intensity",
        "test": "Welch's t-test on log2 median-polish intensities (needs two values per category)",
        "baseline": comparison.a, "versus": comparison.b,
        "baseline_samples": comparison.samples_a, "versus_samples": comparison.samples_b,
        "normalized_on_peptides": quant.quantities.shared_peptides,
        "normalized": quant.quantities.samples.iter().zip(&quant.quantities.normalized)
            .map(|(s, &n)| json!({"biosample": s.name, "normalized": n})).collect::<Vec<_>>(),
        "left_out_incomplete_ms1": quant.incomplete,
        "values_from_match_between_runs": quant.transferred,
        "runs": quant.runs.iter().map(|r| json!({
            "run": r.run, "ms1_scans": r.ms1_scans, "precursors_traced": r.traced,
            "precursors_measured": r.measured, "precursors_untraced": r.untraced,
        })).collect::<Vec<_>>(),
        "total": rows.len(),
        "tested": rows.iter().filter(|r| r.p.is_some()).count(),
        "page": served,
        "significant_at_5pct": rows.iter().filter(|r| r.q.is_some_and(|q| q <= 0.05)).count(),
        "rows": rows[start..end].iter().map(|r| json!({
            "group_id": r.group_id, "accession": r.accession,
            "mean_log2_baseline": r.mean_a, "mean_log2_versus": r.mean_b,
            "measured_baseline": r.measured_a, "measured_versus": r.measured_b,
            "log2_fold": r.log2_fold, "q_value": r.q, "welch_p": r.p,
        })).collect::<Vec<_>>(),
    }))
}

fn reporter_channels_json(quant: &tuffold_engine::ReporterQuant) -> Vec<Value> {
    let r = &quant.quantities;
    let q = &r.quantities;
    r.channels
        .iter()
        .zip(&q.samples)
        .zip(&q.shifts)
        .map(|((channel, sample), &shift)| {
            json!({
                "channel": channel,
                "biosample": sample.name,
                "category": sample.category,
                "log2_shift": -shift,
            })
        })
        .collect()
}

fn better_reagents_json(quant: &tuffold_engine::ReporterQuant) -> Vec<Value> {
    quant
        .better_reagents
        .iter()
        .map(|(plex, other, complete)| {
            json!({
                "plex": plex,
                "label": other.name(),
                "complete_psms": complete,
                "note": format!(
                    "in plex {plex}, at the {other} reagents' reporter m/z, {complete} PSMs have \
                     every channel; redefine it with label \"{other}\" if it was labelled with \
                     those"
                ),
            })
        })
        .collect()
}

/// Compare two categories by TMT reporter ions (Welch's t-test on log2
/// median-normalized reporter intensities), most significant first.
pub fn compare_reporters(
    path: &Path,
    plex: Option<&str>,
    baseline: &str,
    versus: &str,
    offset: usize,
    limit: usize,
) -> Result<Value> {
    let quant =
        tuffold_engine::quantify_reporters(&open(path)?, plex).map_err(|e| e.to_string())?;
    let r = &quant.quantities;
    let q = &r.quantities;
    let comparison = q.compare(baseline, versus).map_err(|e| e.to_string())?;
    let by_id: std::collections::HashMap<i64, &_> =
        q.rows.iter().map(|row| (row.group_id, row)).collect();
    let rows = &comparison.rows;
    let (start, end, served) = page(offset, limit, rows.len());
    Ok(json!({
        "method": "reporters",
        "plexes": quant.plexes.iter().map(plex_json).collect::<Vec<_>>(),
        "bridged": quant.plexes.len() > 1,
        "values_left_out_unbridged": quant.unbridged,
        "psms_used": r.psms_used,
        "psms_incomplete": r.psms_incomplete,
        "spectra_without_peaks": quant.without_peaks,
        "psms_from_ms3": quant.from_ms3,
        "min_precursor_purity": quant.min_purity,
        "spectra_below_purity": quant.impure,
        "spectra_purity_unknown": quant.purity_unknown,
        "better_reagents": better_reagents_json(&quant),
        "channels": reporter_channels_json(&quant),
        "test": "Welch's t-test on log2 reporter intensities (needs two values per category)",
        "baseline": comparison.a, "versus": comparison.b,
        "baseline_samples": comparison.samples_a, "versus_samples": comparison.samples_b,
        "total": rows.len(),
        "tested": rows.iter().filter(|r| r.p.is_some()).count(),
        "page": served,
        "significant_at_5pct": rows.iter().filter(|r| r.q.is_some_and(|q| q <= 0.05)).count(),
        "rows": rows[start..end].iter().map(|s| {
            let row = by_id.get(&s.group_id);
            json!({
                "group_id": s.group_id, "accession": s.accession,
                "psms": row.and_then(|r| r.peptides.first().copied()).unwrap_or_default(),
                "intensities": row.map(|r| &r.intensities),
                "log2": row.map(|r| &r.log2),
                "mean_log2_baseline": s.mean_a, "mean_log2_versus": s.mean_b,
                "measured_baseline": s.measured_a, "measured_versus": s.measured_b,
                "log2_fold": s.log2_fold, "q_value": s.q, "welch_p": s.p,
            })
        }).collect::<Vec<_>>(),
    }))
}

/// Quantify a TMT plex by reporter ions (median-normalized per channel),
/// optionally comparing two categories (`baseline` and `versus`).
pub fn quantify_reporters(
    path: &Path,
    plex: Option<&str>,
    baseline: Option<&str>,
    versus: Option<&str>,
    offset: usize,
    limit: usize,
) -> Result<Value> {
    match (baseline, versus) {
        (Some(baseline), Some(versus)) => {
            return compare_reporters(path, plex, baseline, versus, offset, limit);
        }
        (Some(_), None) | (None, Some(_)) => {
            return Err(
                "give both baseline and versus to compare two categories, or neither to list \
                 reporter intensities"
                    .into(),
            );
        }
        (None, None) => {}
    }
    let quant =
        tuffold_engine::quantify_reporters(&open(path)?, plex).map_err(|e| e.to_string())?;
    let r = &quant.quantities;
    let q = &r.quantities;
    let (start, end, served) = page(offset, limit, q.rows.len());
    Ok(json!({
        "plexes": quant.plexes.iter().map(plex_json).collect::<Vec<_>>(),
        "bridged": quant.plexes.len() > 1,
        "values_left_out_unbridged": quant.unbridged,
        "psms_used": r.psms_used,
        "psms_incomplete": r.psms_incomplete,
        "spectra_without_peaks": quant.without_peaks,
        "psms_from_ms3": quant.from_ms3,
        "min_precursor_purity": quant.min_purity,
        "spectra_below_purity": quant.impure,
        "spectra_purity_unknown": quant.purity_unknown,
        "better_reagents": better_reagents_json(&quant),
        "channels": reporter_channels_json(&quant),
        "total": q.rows.len(),
        "page": served,
        "rows": q.rows[start..end].iter().map(|row| json!({
            "group_id": row.group_id,
            "accession": row.accession,
            "psms": row.peptides.first().copied().unwrap_or_default(),
            "intensities": row.intensities,
            "log2": row.log2,
        })).collect::<Vec<_>>(),
    }))
}

/// Sample QC by `counts`, MS1 `intensity`, or `reporters`.
pub fn sample_qc(path: &Path, method: &str, plex: Option<&str>) -> Result<Value> {
    let project = open(path)?;
    let (qc, values) = match method {
        "counts" => (
            quantify_project(&project)
                .map_err(|e| e.to_string())?
                .sample_qc(),
            "log2(normalized counts + 1)",
        ),
        "intensity" => (
            tuffold_engine::quantify_intensity(&project)
                .map_err(|e| e.to_string())?
                .quantities
                .sample_qc(),
            "log2 median-polish intensities",
        ),
        "reporters" => (
            tuffold_engine::quantify_reporters(&project, plex)
                .map_err(|e| e.to_string())?
                .quantities
                .quantities
                .sample_qc(),
            "log2 reporter intensities",
        ),
        other => {
            return Err(format!(
                "method is counts, intensity, or reporters, not '{other}'"
            ));
        }
    };
    Ok(json!({
        "values": values,
        "groups": qc.groups,
        "samples": qc.samples.iter().enumerate().map(|(s, sample)| json!({
            "name": sample.name, "category": sample.category,
            "measured": sample.measured, "median_log2": sample.median_log2,
            "pc1": qc.pca.as_ref().map(|p| p.scores[s].0),
            "pc2": qc.pca.as_ref().map(|p| p.scores[s].1),
            "correlation": qc.correlation[s],
        })).collect::<Vec<_>>(),
        "pca": qc.pca.as_ref().map(|p| json!({
            "groups_measured_in_all": p.groups,
            "variance_explained": p.explained,
        })),
    }))
}

/// Compare every category at once (one-way ANOVA), by `counts`, MS1
/// `intensity`, or TMT `reporters`, most significant first.
pub fn compare_all_categories(
    path: &Path,
    method: &str,
    plex: Option<&str>,
    offset: usize,
    limit: usize,
) -> Result<Value> {
    let project = open(path)?;
    let (anova, values) = match method {
        "counts" => (
            quantify_project(&project)
                .map_err(|e| e.to_string())?
                .anova(),
            "log2(normalized counts + 1)",
        ),
        "intensity" => (
            tuffold_engine::quantify_intensity(&project)
                .map_err(|e| e.to_string())?
                .quantities
                .anova(),
            "log2 median-polish intensities",
        ),
        "reporters" => (
            tuffold_engine::quantify_reporters(&project, plex)
                .map_err(|e| e.to_string())?
                .quantities
                .quantities
                .anova(),
            "log2 reporter intensities",
        ),
        other => {
            return Err(format!(
                "method is counts, intensity, or reporters, not '{other}'"
            ));
        }
    };
    let anova = anova.map_err(|e| e.to_string())?;
    let (start, end, served) = page(offset, limit, anova.rows.len());
    Ok(json!({
        "test": format!("one-way ANOVA on {values}"),
        "categories": anova.categories.iter().zip(&anova.members)
            .map(|(c, m)| json!({"category": c, "samples": m})).collect::<Vec<_>>(),
        "total": anova.rows.len(),
        "tested": anova.rows.iter().filter(|r| r.p.is_some()).count(),
        "significant_at_5pct": anova.rows.iter().filter(|r| r.q.is_some_and(|q| q <= 0.05)).count(),
        "page": served,
        "rows": anova.rows[start..end].iter().map(|r| json!({
            "group_id": r.group_id, "accession": r.accession,
            "mean_log2": r.means, "measured": r.measured,
            "log2_range": r.range, "f": r.f, "p": r.p, "q_value": r.q,
        })).collect::<Vec<_>>(),
    }))
}

/// Assign biosamples to categories (an empty category clears it).
pub fn set_categories(path: &Path, assignments: &[(String, String)]) -> Result<Value> {
    let mut project = open(path)?;
    let pairs: Vec<(&str, Option<&str>)> = assignments
        .iter()
        .map(|(name, category)| {
            (
                name.as_str(),
                Some(category.trim()).filter(|c| !c.is_empty()),
            )
        })
        .collect();
    tuffold_engine::set_categories(&mut project, &pairs).map_err(|e| e.to_string())?;
    let biosamples = project
        .protein_list()
        .map_err(|e| e.to_string())?
        .biosamples;
    Ok(json!({ "biosamples": biosamples.iter().map(|b| json!({
        "name": b.name, "category": b.category,
    })).collect::<Vec<_>>() }))
}

/// Define (or redefine) an isobaric plex: its reagent set, runs (all runs not
/// in another plex when empty), channel -> biosample assignments, and the
/// lot's impurity certificate (a local file, or "none" to remove it), and
/// the reference channels that bridge it to other plexes (when either is
/// omitted, the plex keeps its own).
pub fn set_plex(
    path: &Path,
    name: &str,
    label: &str,
    runs: &[String],
    channels: &[(String, String)],
    impurities: Option<&str>,
    references: Option<&[String]>,
) -> Result<Value> {
    let plex = tuffold_core::reporter::Plex::parse(label)?;
    let mut project = open(path)?;
    let impurities = match impurities {
        None => tuffold_engine::kept_impurities(&project, name, plex).map_err(|e| e.to_string())?,
        Some("none") => None,
        Some(file) => {
            let text = std::fs::read_to_string(file).map_err(|e| format!("{file}: {e}"))?;
            Some(
                tuffold_core::reporter::Impurities::parse(&text, plex)
                    .map_err(|e| format!("{file}: {e}"))?,
            )
        }
    };
    let definition = tuffold_store::PlexDefinition {
        name: name.to_string(),
        plex,
        runs: runs.to_vec(),
        channels: channels.to_vec(),
        impurities,
        references: match references {
            Some(given) => given.to_vec(),
            None => {
                tuffold_engine::kept_references(&project, name, plex).map_err(|e| e.to_string())?
            }
        },
    };
    tuffold_engine::set_plex(&mut project, &definition).map_err(|e| e.to_string())?;
    plexes_json(&project)
}

/// Remove a named plex (its channel biosamples stay).
pub fn remove_plex(path: &Path, name: &str) -> Result<Value> {
    let mut project = open(path)?;
    tuffold_engine::remove_plex(&mut project, name).map_err(|e| e.to_string())?;
    plexes_json(&project)
}

fn plexes_json(project: &Project) -> Result<Value> {
    let plexes = project.plexes().map_err(|e| e.to_string())?;
    let biosamples = project
        .protein_list()
        .map_err(|e| e.to_string())?
        .biosamples;
    Ok(json!({
        "plexes": plexes.iter().map(plex_json).collect::<Vec<_>>(),
        "biosamples": biosamples.iter().map(|b| json!({
            "name": b.name, "category": b.category,
        })).collect::<Vec<_>>(),
    }))
}

/// Match between runs in the runs of the given spectra files.
pub fn match_between_runs(path: &Path, files: &[PathBuf]) -> Result<Value> {
    let mut project = open(path)?;
    let reports =
        tuffold_engine::match_between_runs(&mut project, files).map_err(|e| e.to_string())?;
    Ok(json!({
        "fdr": tuffold_store::TRANSFER_FDR,
        "runs": reports.iter().map(|r| json!({
            "run": r.run,
            "other_runs_aligned": r.aligned,
            "other_runs": r.donors,
            "looked_for": r.summary.candidates,
            "transferred": r.summary.accepted,
            "pairs_won_by_targets": r.summary.targets_won,
            "pairs_won_by_decoys": r.summary.decoys_won,
            "known_peaks": r.summary.known,
            "model_trained": r.summary.trained,
            "ms1_scans_damaged": r.ms1_scans_damaged,
        })).collect::<Vec<_>>(),
    }))
}

/// Import search results and spectra files (creating the project if
/// needed). Importing clears validation results; validate afterwards.
pub fn import_data(path: &Path, inputs: &[PathBuf]) -> Result<Value> {
    importing(path, |project| {
        let outcomes = import_inputs(project, inputs, None);
        let imported = outcomes.iter().any(|o| !o.is_failure());
        Ok((outcomes_json(&outcomes), imported))
    })
}

fn outcomes_json(outcomes: &[tuffold_engine::Outcome]) -> Value {
    json!({ "outcomes": outcomes.iter().map(|o| {
        let input = o.input.display().to_string();
        match &o.kind {
            OutcomeKind::Imported { software, report, warnings, notes, .. } => json!({
                "input": input, "status": "imported", "software": software,
                "psms": report.psms, "new_spectra": report.new_spectra,
                "new_peptides": report.new_peptides, "new_proteins": report.new_proteins,
                "warnings": warnings.iter().map(|w| w.message.clone())
                    .chain(duplicate_note(report)).collect::<Vec<_>>(),
                "notes": notes,
            }),
            OutcomeKind::Attached(report) => json!({
                "input": input, "status": "attached", "run": report.run,
                "read": report.read, "attached": report.attached,
                "identified_without_peaks": report.identified_without_peaks,
                "ms1_scans": report.ms1_scans,
                "ms1_scans_damaged": report.ms1_scans_damaged,
                "precursors_traced": report.precursors_traced,
                "precursors_measured": report.precursors_measured,
                "ms1": tuffold_engine::ms1_note(report),
                "ms3_scans": report.ms3_scans,
                "ms3_attached": report.ms3_attached,
                "ms3": tuffold_engine::ms3_note(report),
            }),
            OutcomeKind::Failed(message) => json!({
                "input": input, "status": "failed", "error": message,
            }),
        }
    }).collect::<Vec<_>>() })
}

/// Validate with the given thresholds (rescoring each engine, or ranking by
/// its best single score).
pub fn validate(path: &Path, request: &ThresholdRequest, rescore: bool) -> Result<Value> {
    let mut project = open(path)?;
    let thresholds = saved_thresholds(&project)?.adjusted(request)?;
    let summary = validate_project(
        &mut project,
        &Settings {
            rescore,
            thresholds,
        },
    )
    .map_err(|e| e.to_string())?;
    Ok(validation_json(&summary))
}

/// Re-apply thresholds without validating again (results must be current).
pub fn set_thresholds(path: &Path, request: &ThresholdRequest) -> Result<Value> {
    let mut project = open(path)?;
    let thresholds = saved_thresholds(&project)?.adjusted(request)?;
    let summary = apply_thresholds(&mut project, &thresholds).map_err(|e| e.to_string())?;
    Ok(validation_json(&summary))
}

/// The project's saved thresholds (the defaults before any validation), so a
/// request changes only what it names.
fn saved_thresholds(project: &Project) -> Result<Thresholds> {
    Ok(project
        .validation_summary()
        .map_err(|e| e.to_string())?
        .map(|summary| summary.settings.thresholds)
        .unwrap_or_default())
}

/// Search spectra files against a FASTA with the built-in search (Sage's
/// core) and import the PSMs; validate afterwards.
pub fn search_spectra(
    path: &Path,
    fasta: &Path,
    spectra: &[PathBuf],
    settings: SearchSettings,
) -> Result<Value> {
    importing(path, |project| {
        let described = settings.describe();
        let report = search_files(project, fasta, spectra, settings).map_err(|e| e.to_string())?;
        let imported = report.runs.iter().any(|run| {
            run.outcome
                .as_ref()
                .is_ok_and(|done| done.imported.is_some())
        });
        let mut json = search_json(&report);
        json["search"] = Value::String(described);
        Ok((json, imported))
    })
}

fn search_json(report: &tuffold_engine::SearchReport) -> Value {
    json!({
        "target_proteins": report.targets,
        "peptide_forms": report.peptides,
        "index_estimate_gb": report.index_bytes as f64 / 1e9,
        "memory_budget_gb": report.memory_budget as f64 / 1e9,
        "fasta_decoys_replaced": report.dropped_decoys,
        "fasta_entries_left_out": report.excluded_targets,
        "runs": report.runs.iter().map(|run| match &run.outcome {
            Ok(done) => json!({
                "input": run.input.display().to_string(), "run": run.run,
                "status": if done.imported.is_some() { "searched" } else { "no_matches" },
                "searched": done.searched,
                "not_searched": done.skipped.total(),
                "not_searched_reasons": done.skipped.describe(),
                "profile_spectra_centroided": done.centroided,
                "psms_imported": done.imported.as_ref().map_or(0, |i| i.import.psms),
                "peaks_attached": done.imported.as_ref().map_or(0, |i| i.attach.attached),
                "ms1": done.imported.as_ref().map(|i| tuffold_engine::ms1_note(&i.attach)),
            }),
            Err(message) => json!({
                "input": run.input.display().to_string(), "run": run.run,
                "status": "failed", "error": message,
            }),
        }).collect::<Vec<_>>(),
    })
}

/// Undo the last import, search, or sample-group change; called again, redo
/// it.
pub fn undo_last_change(path: &Path) -> Result<Value> {
    let path = absolute(path)?;
    let note = tuffold_engine::undo(path).map_err(|e| e.to_string())?;
    Ok(json!({
        "done": if note.redo { "redid" } else { "undid" },
        "change": note.label,
        "next": if note.redo { "calling undo again undoes it" } else { "calling undo again redoes it" },
    }))
}

/// Create the synthetic demo project (labelled as such) to explore with.
pub fn create_demo(path: &Path) -> Result<Value> {
    let path = absolute(path)?;
    if path.exists() {
        return Err(format!("{} already exists", path.display()));
    }
    let summary = create_demo_project(path).map_err(|e| e.to_string())?;
    Ok(json!({
        "project": path.display().to_string(),
        "synthetic": true,
        "validation": validation_json(&summary),
    }))
}

/// Write a validated project's Excel workbook to `to` (absolute, .xlsx,
/// not yet existing).
pub fn export_workbook(path: &Path, to: &Path) -> Result<Value> {
    let to = absolute(to)?;
    if !to
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("xlsx"))
    {
        return Err(format!(
            "{}: the workbook needs an .xlsx name",
            to.display()
        ));
    }
    if to.exists() {
        return Err(format!(
            "{} already exists; choose a new file name",
            to.display()
        ));
    }
    let written = tuffold_engine::write_workbook(&open(path)?, to).map_err(|e| e.to_string())?;
    Ok(json!({
        "workbook": to.display().to_string(),
        "sheets": written.sheets.iter().map(|(name, rows)| json!({"sheet": name, "rows": rows}))
            .collect::<Vec<_>>(),
        "cut_at_excel_row_limit": written.cut.iter()
            .map(|(name, total)| json!({"sheet": name, "rows_in_project": total}))
            .collect::<Vec<_>>(),
    }))
}

/// Write the results as mzIdentML 1.2 (.mzid, absolute path, not yet
/// existing).
pub fn export_mzid(path: &Path, to: &Path) -> Result<Value> {
    let to = absolute(to)?;
    if !to
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("mzid"))
    {
        return Err(format!("{}: the results need an .mzid name", to.display()));
    }
    if to.exists() {
        return Err(format!(
            "{} already exists; choose a new file name",
            to.display()
        ));
    }
    let project = open(path)?;
    let file = std::fs::File::create(to).map_err(|e| format!("{}: {e}", to.display()))?;
    let mut out = std::io::BufWriter::new(file);
    let written = tuffold_engine::write_mzid(&project, &mut out)
        .map_err(|e| e.to_string())
        .and_then(|report| {
            std::io::Write::flush(&mut out).map_err(|e| e.to_string())?;
            Ok(report)
        });
    match written {
        Ok(report) => Ok(json!({
            "mzidentml": to.display().to_string(),
            "spectra": report.spectra,
            "passing_target_spectra": report.passing_spectra,
            "protein_groups": report.protein_groups,
            "passing_target_groups": report.passing_groups,
        })),
        Err(message) => {
            drop(out);
            let _ = std::fs::remove_file(to);
            Err(message)
        }
    }
}

/// Write the spectral library (.msp, absolute path, not yet existing).
pub fn export_library(path: &Path, to: &Path) -> Result<Value> {
    let to = absolute(to)?;
    if !to
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("msp"))
    {
        return Err(format!("{}: the library needs an .msp name", to.display()));
    }
    if to.exists() {
        return Err(format!(
            "{} already exists; choose a new file name",
            to.display()
        ));
    }
    let project = open(path)?;
    let file = std::fs::File::create(to).map_err(|e| format!("{}: {e}", to.display()))?;
    let mut out = std::io::BufWriter::new(file);
    let written = tuffold_engine::write_msp(&project, &mut out).and_then(|report| {
        std::io::Write::flush(&mut out)?;
        Ok(report)
    });
    match written {
        Ok(report) => Ok(json!({
            "library": to.display().to_string(),
            "peptide_ions": report.entries,
            "left_out_without_peaks": report.without_peaks,
            "left_out_unnamed_modifications": report.unnamed_modifications,
        })),
        Err(error) => {
            drop(out);
            let _ = std::fs::remove_file(to);
            Err(error.to_string())
        }
    }
}

/// A methods paragraph for publication, from what the project did.
pub fn methods_text(path: &Path) -> Result<Value> {
    let methods = tuffold_engine::methods(&open(path)?).map_err(|e| e.to_string())?;
    Ok(json!({
        "paragraphs": methods.paragraphs,
        "references": methods.references,
        "text": methods.text(),
    }))
}

#[cfg(test)]
mod tests;
