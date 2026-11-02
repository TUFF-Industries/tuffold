//! `tuffold validate` and `tuffold thresholds`.

use std::path::Path;

use clap::Args;
use tuffold_engine::{DeepLearning, apply_thresholds, validate_project_with};
use tuffold_store::{LevelCounts, Project, ValidationSummary};
use tuffold_validate::{
    Level, Settings, ThresholdMode, ThresholdRequest, Thresholds, format_rate as percent,
};

/// Threshold flags shared by `validate` and `thresholds`; unset flags keep
/// the project's saved value (1% FDR and 1 peptide by default). FDR flags
/// and probability flags choose the mode; they can't be mixed.
#[derive(Args, Debug, Default)]
pub struct ThresholdArgs {
    /// Spectrum (PSM) FDR, as a percentage ("1%") or fraction ("0.01").
    #[arg(long, value_parser = parse_rate)]
    psm_fdr: Option<f64>,
    /// Peptide FDR; protein groups are built from peptides passing it.
    #[arg(long, value_parser = parse_rate)]
    peptide_fdr: Option<f64>,
    /// Protein group FDR (picked target-decoy competition).
    #[arg(long, value_parser = parse_rate)]
    protein_fdr: Option<f64>,
    /// Spectrum (PSM) probability, 1 − PEP ("95%" or "0.95").
    #[arg(long, value_parser = parse_rate)]
    psm_prob: Option<f64>,
    /// Peptide probability.
    #[arg(long, value_parser = parse_rate)]
    peptide_prob: Option<f64>,
    /// Protein group probability.
    #[arg(long, value_parser = parse_rate)]
    protein_prob: Option<f64>,
    /// Threshold by `fdr` (q-values) or `probability` (1 − PEP; defaults 95%
    /// spectra and peptides, 99% protein groups).
    #[arg(long, value_name = "MODE")]
    by: Option<String>,
    /// Distinct peptide sequences (unique or razor) a protein group needs.
    #[arg(long)]
    min_peptides: Option<u32>,
}

impl ThresholdArgs {
    fn over(&self, saved: Thresholds) -> Result<Thresholds, String> {
        let by = match self.by.as_deref() {
            None => None,
            Some(text) => Some(
                ThresholdMode::parse(text)
                    .ok_or_else(|| format!("--by takes fdr or probability, not '{text}'"))?,
            ),
        };
        saved.adjusted(&ThresholdRequest {
            by,
            fdr: [self.psm_fdr, self.peptide_fdr, self.protein_fdr],
            probability: [self.psm_prob, self.peptide_prob, self.protein_prob],
            min_peptides: self.min_peptides,
        })
    }
}

/// "1%" or "0.01" -> 0.01; anything outside (0, 1) is refused.
pub fn parse_rate(text: &str) -> Result<f64, String> {
    let trimmed = text.trim();
    let (number, scale) = match trimmed.strip_suffix('%') {
        Some(percent) => (percent.trim(), 100.0),
        None => (trimmed, 1.0),
    };
    let value = number
        .parse::<f64>()
        .map_err(|_| format!("'{text}' is not a rate; write 1% or 0.01"))?
        / scale;
    if value > 0.0 && value < 1.0 {
        Ok(value)
    } else {
        Err(format!(
            "'{text}' is not between 0% and 100%; write 1% or 0.01"
        ))
    }
}

pub fn validate(
    path: &Path,
    no_rescore: bool,
    no_deep: bool,
    args: &ThresholdArgs,
) -> Result<(), String> {
    let mut project = Project::open(path).map_err(|error| error.to_string())?;
    let saved = project
        .validation_summary()
        .map_err(|error| error.to_string())?
        .map(|summary| summary.settings)
        .unwrap_or_default();
    let settings = Settings {
        rescore: !no_rescore,
        thresholds: args.over(saved.thresholds)?,
    };
    settings.thresholds.check()?;
    let deep = if no_deep {
        DeepLearning::Off
    } else {
        DeepLearning::Auto
    };
    let summary =
        validate_project_with(&mut project, &settings, deep).map_err(|error| error.to_string())?;
    print_validation(&summary);
    Ok(())
}

pub fn thresholds(path: &Path, args: &ThresholdArgs) -> Result<(), String> {
    let mut project = Project::open(path).map_err(|error| error.to_string())?;
    let saved = project
        .validation_summary()
        .map_err(|error| error.to_string())?
        .map(|summary| summary.settings.thresholds)
        .unwrap_or_default();
    let summary =
        apply_thresholds(&mut project, &args.over(saved)?).map_err(|error| error.to_string())?;
    print_validation(&summary);
    Ok(())
}

pub fn print_validation(summary: &ValidationSummary) {
    if summary.validated_at.is_none() {
        out!("validation: stale (data was imported since); run `tuffold validate` to update it");
        return;
    }
    let rescore = if summary.settings.rescore {
        "on"
    } else {
        "off"
    };
    out!("validation (rescoring {rescore}):");
    for engine in &summary.engines {
        out!(
            "  engine {}: {} rank-1 PSMs ({} decoy), {} targets at 1% FDR alone",
            engine.name,
            engine.psms,
            engine.decoys,
            engine.accepted
        );
        out!("    {}", engine.detail);
        if engine.method == "rescored" {
            out!("    features: {}", engine.features.join(", "));
        }
    }
    let thresholds = summary.settings.thresholds;
    let kind = match thresholds.mode {
        ThresholdMode::Fdr => "FDR",
        ThresholdMode::Probability => "probability",
    };
    out!(
        "  thresholds: PSM {kind} {}, peptide {kind} {}, protein {kind} {}, at least {} \
         peptide(s) per group",
        percent(thresholds.psm),
        percent(thresholds.peptide),
        percent(thresholds.protein),
        thresholds.min_peptides
    );
    print_level("spectra", &summary.spectra, &thresholds, Level::Psm);
    print_level("peptides", &summary.peptides, &thresholds, Level::Peptide);
    if let Some(scoring) = &summary.peptide_scoring {
        out!("    {}", scoring.detail);
    }
    print_level(
        "protein groups",
        &summary.protein_groups,
        &thresholds,
        Level::Protein,
    );
    print_localization(&summary.localization);
    for note in &summary.notes {
        out!("note [{}] {}", note.code.as_str(), note.message);
    }
}

pub fn print_localization(loc: &tuffold_validate::LocalizationSummary) {
    if loc.scored_psms == 0
        && loc.missing_peaks_psms == 0
        && loc.capped_psms == 0
        && loc.unscorable_psms == 0
    {
        return;
    }
    let flr = loc.estimated_flr.map_or_else(String::new, |f| {
        format!(", decoy-AA estimated FLR {}", percent(f))
    });
    out!(
        "  PTM localization: {} of {} ambiguous passing PSMs ({} of {} sites) localized at >= {} site probability{flr}",
        loc.localized_psms,
        loc.scored_psms,
        loc.localized_sites,
        loc.total_target_sites,
        percent(loc.probability_threshold),
    );
    if loc.missing_peaks_psms > 0 {
        out!(
            "    {} ambiguous passing PSMs without peaks (attach spectra files, then validate again, to localize them)",
            loc.missing_peaks_psms
        );
    }
    if loc.capped_psms > 0 {
        out!(
            "    {} ambiguous passing PSMs over 256 isoforms (left unscored)",
            loc.capped_psms
        );
    }
    if loc.unscorable_psms > 0 {
        out!(
            "    {} ambiguous passing PSMs whose isoforms can't be fragmented (left unscored)",
            loc.unscorable_psms
        );
    }
}

pub fn localize(path: &Path, site_prob: Option<f64>, export: Option<&Path>) -> Result<(), String> {
    let project = Project::open(path).map_err(|error| error.to_string())?;
    let threshold = site_prob.unwrap_or(tuffold_validate::DEFAULT_SITE_PROBABILITY);
    let loc = project
        .localization_summary(threshold)
        .map_err(|error| error.to_string())?
        .ok_or_else(|| {
            "this project has not been validated since its data last changed; validate it first"
                .to_string()
        })?;
    if loc.scored_psms == 0
        && loc.missing_peaks_psms == 0
        && loc.capped_psms == 0
        && loc.unscorable_psms == 0
    {
        out!("PTM localization: no passing PSMs carry ambiguous variable modifications");
    } else {
        print_localization(&loc);
        // Which fragment tolerance was used is a validation note.
        let summary = project
            .validation_summary()
            .map_err(|error| error.to_string())?;
        for note in summary.iter().flat_map(|summary| &summary.notes) {
            if note.code == tuffold_validate::NoteCode::Localization {
                out!("    {}", note.message);
            }
        }
    }
    if let Some(to) = export {
        let rows = crate::write_new(to, |out| {
            tuffold_engine::write_localization_tsv(out, &project)
        })?;
        out!("wrote {}: {rows} localized PSMs", to.display());
    }
    Ok(())
}

fn print_level(name: &str, counts: &LevelCounts, thresholds: &Thresholds, level: Level) {
    out!(
        "  {name}: {} targets pass ({} decoys)",
        counts.targets,
        counts.decoys
    );
    if thresholds.mode == ThresholdMode::Probability {
        // What a probability threshold means in FDR terms, from the decoys.
        if let Some(fdr) = counts.estimated_fdr() {
            out!(
                "    estimated FDR {} ((decoys + 1) / targets passing)",
                percent(fdr)
            );
        }
        return;
    }
    let threshold = thresholds.value(level);
    if let Some(lowest) = counts.lowest_q.filter(|lowest| *lowest > threshold) {
        out!(
            "    note: no {name} can pass {}; the lowest q-value reached is {} (with few \
             entries, target-decoy estimates cannot go lower)",
            percent(threshold),
            percent(lowest)
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn threshold_flags_become_a_request() {
        let saved = Thresholds::default();
        let flags = ThresholdArgs {
            protein_prob: Some(0.9),
            ..ThresholdArgs::default()
        };
        let switched = flags.over(saved).unwrap();
        assert_eq!(
            (switched.mode, switched.protein),
            (ThresholdMode::Probability, 0.9)
        );
        let unknown = ThresholdArgs {
            by: Some("vibes".into()),
            ..ThresholdArgs::default()
        };
        assert!(
            unknown
                .over(saved)
                .unwrap_err()
                .contains("fdr or probability")
        );
    }

    #[test]
    fn rates_accept_percentages_and_fractions() {
        assert_eq!(parse_rate("1%"), Ok(0.01));
        assert_eq!(parse_rate(" 5 % "), Ok(0.05));
        assert_eq!(parse_rate("0.01"), Ok(0.01));
        for bad in ["1", "0", "100%", "-1%", "abc", "%"] {
            assert!(parse_rate(bad).is_err(), "{bad}");
        }
    }
}
