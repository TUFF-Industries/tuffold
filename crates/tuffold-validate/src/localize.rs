//! PTM site localization: per-site probabilities (phosphoRS / ptmRS style),
//! pairwise Ascores (Beausoleil et al. 2006), and decoy-amino-acid false
//! localization rate (FLR) estimation (Ramsbottom et al. 2022).
//!
//! For each passing rank-1 PSM with one or more variable modifications that
//! could sit on more than one candidate residue (`m > k`), candidate isoforms
//! are enumerated up to [`MAX_ISOFORMS`].
//!
//! ## Scoring model
//! 1. **Precursor and precursor neutral-loss exclusion**: peaks matching the
//!    precursor $[M + z\text{H}]^{z+}$ or its neutral losses ($-\text{H}_2\text{O}$,
//!    $-\text{NH}_3$, $-\text{H}_3\text{PO}_4$, $-\text{HPO}_3$,
//!    $-2\text{H}_3\text{PO}_4$) within fragment tolerance are removed before
//!    windowing so intense precursor losses do not displace backbone ions.
//!    Fragment neutral losses ($-\text{H}_3\text{PO}_4$, $-97.9769\text{ Da}$)
//!    are deliberately not scored as site-determining: $\text{Phospho} -
//!    \text{H}_3\text{PO}_4 = -18.0106\text{ Da}$ is mass-identical to ordinary
//!    water loss ($-\text{H}_2\text{O}$) from an unmodified residue.
//! 2. **Peak-depth windowing**: the spectrum is partitioned into $100\text{ }m/z$
//!    windows and peaks in each window are ranked $1 \dots 10$ by descending
//!    intensity (ties by ascending $m/z$).
//! 3. **Pairwise site-determining binomial score**: for any two isoforms $S_a$
//!    and $S_b$, let $I(S_a, S_b)$ be the theoretical $b$ and $y$ ions
//!    ([`tuffold_core::fragments::fragment_ions`]) whose $m/z$ differs between
//!    $S_a$ and $S_b$ ($n = |I(S_a, S_b)|$). At depth $d \in 1..=10$, with trial
//!    probability $p_d = 0.01 \times d$ ($d$ peaks per 100 nominal-mass slots in
//!    a $100\text{ }m/z$ window) and $k_d$ matched site-determining ions:
//!    $$s_d = -10 \log_{10} \sum_{j=k_d}^{n} \binom{n}{j} p_d^j (1 - p_d)^{n - j}$$
//!    Weighted across depths $d = 1 \dots 10$ with Beausoleil weights
//!    $w = [0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 0.75, 0.5, 0.25, 0.25]$ ($\sum w_d = 7$):
//!    $$\text{score}_{a,b}(S) = \frac{1}{7} \sum_{d=1}^{10} w_d s_d(S)$$
//! 4. **Site probabilities and Ascores**: relative to the winning isoform $S^*$,
//!    each isoform $S_i$ has pairwise score difference
//!    $\Delta(S^*, S_i) = \text{score}_{*,i}(S^*) - \text{score}_{*,i}(S_i)$,
//!    isoform probability
//!    $$P_{\text{iso}}(S_i) = \frac{10^{-\Delta(S^*, S_i)/10}}{\sum_j 10^{-\Delta(S^*, S_j)/10}}$$
//!    and per-site probability $P_{\text{site}}(r) = \sum_{S \ni r} P_{\text{iso}}(S)$.
//!    The site's Ascore is $\max(0, \text{score}_{+,-}(S_{+r}^*) - \text{score}_{+,-}(S_{-r}^*))$
//!    between the best isoform carrying $r$ and the best isoform without $r$.
//!    For 2 candidate sites and 1 modification, $P_{\text{site}} = 1 / (1 + 10^{-\text{Ascore}/10})$
//!    holds as an exact identity ($\text{Ascore} = 13 \iff 95.2\%$,
//!    $\text{Ascore} = 20 \iff 99.0\%$).

use std::fmt::Write as _;

use tuffold_core::chemistry::{PROTON, WATER_MONO, peptide_mass};
use tuffold_core::fragments::{FragmentIon, Tolerance, fragment_ions};
use tuffold_core::peptide::{DELTA_DECIMALS, ModSite, Modification, PeptideForm};
use tuffold_core::unimod;

/// Maximum number of candidate isoforms enumerated for a single PSM.
pub const MAX_ISOFORMS: usize = 256;

/// Default probability cut for calling a site localized (95%, Ascore ≈ 13).
pub const DEFAULT_SITE_PROBABILITY: f64 = 0.95;

/// Number of peak-depth levels per 100 m/z window (Beausoleil et al. 2006).
const MAX_DEPTH: usize = 10;

/// Beausoleil et al. (2006) weights for peak depths 1..=10 (sum = 7.0).
const DEPTH_WEIGHTS: [f64; MAX_DEPTH] = [0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 0.75, 0.5, 0.25, 0.25];
const DEPTH_WEIGHT_SUM: f64 = 7.0;

/// Depth used to report matched/total site-determining ion counts.
const REPORT_DEPTH: usize = 6;

const AMMONIA_MONO: f64 = 17.026_549;
const H3PO4_MONO: f64 = 97.976_896;
const HPO3_MONO: f64 = 79.966_331;

/// Input for localizing modifications on one spectrum.
#[derive(Debug, Clone)]
pub struct LocalizationInput<'a> {
    pub spectrum_id: i64,
    pub form: &'a PeptideForm,
    pub charge: Option<u8>,
    pub mz: &'a [f64],
    pub intensity: &'a [f32],
}

/// Confidence score for one candidate residue site and modification mass.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteScore {
    /// 1-based residue index in `form.sequence`.
    pub position: u32,
    /// Single-letter amino acid code at `position`.
    pub residue: char,
    /// Snapped or rounded modification mass delta.
    pub mass_delta: f64,
    /// Canonical modification name if known (e.g. `"Phospho"`, `"Oxidation"`).
    pub modification_name: Option<String>,
    /// Posterior probability in `[0.0, 1.0]` that this site carries the modification.
    pub probability: f64,
    /// Beausoleil Ascore (`-10 log10 P` difference on pairwise site-determining
    /// ions between the best isoform with this site modified vs unmodified).
    pub ascore: f64,
    /// Number of pairwise site-determining ions matched for the winning isoform
    /// at depth 6 out of the total pairwise site-determining ions.
    pub determining_ions_matched: u32,
    pub determining_ions_total: u32,
    /// Whether this site is modified in the top-scoring isoform (`localized_peptide`).
    pub assigned: bool,
    /// Whether this entry is a decoy-amino-acid control site used for FLR estimation.
    pub is_decoy: bool,
}

/// Outcome of localizing one PSM.
#[derive(Debug, Clone, PartialEq)]
pub enum SpectrumLocalization {
    /// No variable modifications, or only 1 possible isoform (`m_T == k`).
    Unambiguous,
    /// Variable modification could sit on multiple sites, but the spectrum has
    /// no stored MS2 peaks.
    MissingPeaks { isoform_count: usize },
    /// Number of candidate isoforms exceeded [`MAX_ISOFORMS`].
    Capped { isoform_count: usize },
    /// Ambiguous, but an isoform could not be built or fragmented (a residue
    /// or modification without a mass, e.g. `X`); counted, never scored.
    Unscorable { isoform_count: usize },
    /// Successfully scored across candidate isoforms.
    Scored(ScoredLocalization),
}

/// Scored localization result for one ambiguous PSM.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoredLocalization {
    pub spectrum_id: i64,
    pub isoform_count: usize,
    /// ProForma key of the highest-scoring isoform.
    pub localized_peptide: String,
    /// Highest Ascore among assigned sites on the best isoform.
    pub best_ascore: f64,
    /// Lowest site probability among assigned sites on the best isoform.
    pub min_probability: f64,
    /// Number of false target site opportunities (`sum(m_T - k)`) on this PSM.
    pub false_target_sites: u32,
    /// Target and decoy candidate site scores, ordered by `(is_decoy, position, mass_delta)`.
    pub sites: Vec<SiteScore>,
}

impl ScoredLocalization {
    /// Target candidate sites only (`!is_decoy`).
    pub fn target_sites(&self) -> impl Iterator<Item = &SiteScore> {
        self.sites.iter().filter(|s| !s.is_decoy)
    }

    /// Compact summary of target site probabilities, e.g. `"S2: 98.8%, T4: 1.2%"`.
    pub fn format_site_probabilities(&self) -> String {
        format_sites(self.sites.iter().filter(|s| !s.is_decoy))
    }
}

/// Format target site scores as `"S2: 98.8%, T4: 1.2%"` (or `"S2(+79.9663): 98.8%"`
/// when more than one distinct modification mass is localized on the same peptide).
pub fn format_sites<'a, I>(sites: I) -> String
where
    I: IntoIterator<Item = &'a SiteScore>,
{
    let target: Vec<&SiteScore> = sites.into_iter().filter(|s| !s.is_decoy).collect();
    if target.is_empty() {
        return String::new();
    }
    let first_delta = target[0].mass_delta;
    let multi_mod = target
        .iter()
        .any(|s| (s.mass_delta - first_delta).abs() > unimod::SNAP_TOLERANCE);
    let mut out = String::new();
    for (idx, site) in target.iter().enumerate() {
        if idx > 0 {
            out.push_str(", ");
        }
        let pct = (site.probability * 100.0).clamp(0.0, 100.0);
        if multi_mod {
            if let Some(name) = &site.modification_name {
                let _ = write!(
                    out,
                    "{}{}({}): {:.1}%",
                    site.residue, site.position, name, pct
                );
            } else {
                let _ = write!(
                    out,
                    "{}{}({:+.prec$}): {:.1}%",
                    site.residue,
                    site.position,
                    site.mass_delta,
                    pct,
                    prec = DELTA_DECIMALS
                );
            }
        } else {
            let _ = write!(out, "{}{}: {:.1}%", site.residue, site.position, pct);
        }
    }
    out
}

/// Dataset-wide summary of PTM site localization and decoy-amino-acid FLR.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LocalizationSummary {
    /// Passing target PSMs with ambiguous PTM sites that were scored.
    pub scored_psms: usize,
    /// Passing target PSMs with ambiguous PTM sites that had at least every
    /// assigned site passing `probability_threshold`.
    pub localized_psms: usize,
    /// Passing target PSMs with ambiguous PTM sites that lacked MS2 peaks.
    pub missing_peaks_psms: usize,
    /// Passing target PSMs whose candidate isoform count exceeded [`MAX_ISOFORMS`].
    pub capped_psms: usize,
    /// Passing target PSMs with ambiguous sites whose isoforms could not be
    /// fragmented (see [`SpectrumLocalization::Unscorable`]).
    pub unscorable_psms: usize,
    /// Target candidate sites with `probability >= probability_threshold`.
    pub localized_sites: usize,
    /// Total target candidate sites across scored passing PSMs.
    pub total_target_sites: usize,
    /// Site probability cut used for `localized_psms`, `localized_sites`, and `estimated_flr`.
    pub probability_threshold: f64,
    /// Estimated false localization rate (FLR) at `probability_threshold` via
    /// the decoy amino acid method, or `None` when no decoy sites were tested
    /// or no target sites pass.
    pub estimated_flr: Option<f64>,
}

/// Estimate the false localization rate (FLR) at `probability_threshold` from
/// scored PSMs (Ramsbottom et al., *J Proteome Res* 2022):
///
/// $$\widehat{\text{FLR}}(\tau) = \min\left(1,\; \frac{T_{\text{false}} \cdot (N_X(\tau) + 1)}{X_{\text{tested}} \cdot N_T(\tau)}\right)$$
///
/// where $T_{\text{false}} = \sum (m_T - k)$ is the number of false target site
/// opportunities, $X_{\text{tested}}$ is the number of decoy sites evaluated,
/// $N_T(\tau)$ is the number of target sites with $P \ge \tau$, and $N_X(\tau)$
/// is the number of decoy sites with $P \ge \tau$, counted as $N_X(\tau) + 1$
/// (the +1 every other Tuffold FDR carries): with no passing decoy the estimate
/// is the method's resolution, not 0.
pub fn estimate_flr(
    false_target_sites: usize,
    decoy_sites_tested: usize,
    passing_target_sites: usize,
    passing_decoy_sites: usize,
) -> Option<f64> {
    if decoy_sites_tested == 0 || passing_target_sites == 0 {
        return None;
    }
    let num = (false_target_sites as f64) * ((passing_decoy_sites + 1) as f64);
    let den = (decoy_sites_tested as f64) * (passing_target_sites as f64);
    Some((num / den).clamp(0.0, 1.0))
}

/// Summarize a batch of [`SpectrumLocalization`] outputs at `probability_threshold`.
pub fn summarize_localizations<'a, I>(results: I, probability_threshold: f64) -> LocalizationSummary
where
    I: IntoIterator<Item = &'a SpectrumLocalization>,
{
    let mut summary = LocalizationSummary {
        probability_threshold,
        ..LocalizationSummary::default()
    };
    let mut false_target_sites = 0usize;
    let mut decoy_sites_tested = 0usize;
    let mut passing_decoy_sites = 0usize;

    for result in results {
        match result {
            SpectrumLocalization::Unambiguous => {}
            SpectrumLocalization::MissingPeaks { .. } => {
                summary.missing_peaks_psms += 1;
            }
            SpectrumLocalization::Capped { .. } => {
                summary.capped_psms += 1;
            }
            SpectrumLocalization::Unscorable { .. } => {
                summary.unscorable_psms += 1;
            }
            SpectrumLocalization::Scored(scored) => {
                summary.scored_psms += 1;
                if scored.min_probability >= probability_threshold {
                    summary.localized_psms += 1;
                }
                false_target_sites += scored.false_target_sites as usize;
                for site in &scored.sites {
                    if site.is_decoy {
                        decoy_sites_tested += 1;
                        if site.probability >= probability_threshold {
                            passing_decoy_sites += 1;
                        }
                    } else {
                        summary.total_target_sites += 1;
                        if site.probability >= probability_threshold {
                            summary.localized_sites += 1;
                        }
                    }
                }
            }
        }
    }

    summary.estimated_flr = estimate_flr(
        false_target_sites,
        decoy_sites_tested,
        summary.localized_sites,
        passing_decoy_sites,
    );
    summary
}

/// Returns `Some(isoform_count)` (with `isoform_count >= 2`) when `form` has
/// variable modifications that could sit on more than one candidate residue,
/// or `None` when the peptide is unmodified or has only one possible placement.
pub fn ambiguous_isoform_count(form: &PeptideForm) -> Option<usize> {
    if form.modifications.is_empty() {
        return None;
    }
    let residues: Vec<char> = form.sequence.chars().collect();
    let plan = build_modification_plan(form, &residues)?;
    let count = count_joint_isoforms(&plan.groups);
    (count >= 2).then_some(count)
}

/// Localize variable modifications on a single spectrum.
pub fn localize_spectrum(
    input: &LocalizationInput<'_>,
    tolerance: Tolerance,
) -> SpectrumLocalization {
    assert_eq!(
        input.mz.len(),
        input.intensity.len(),
        "m/z and intensity slices must have equal length"
    );
    let residues: Vec<char> = input.form.sequence.chars().collect();
    if residues.is_empty() {
        return SpectrumLocalization::Unambiguous;
    }

    let Some(plan) = build_modification_plan(input.form, &residues) else {
        return SpectrumLocalization::Unambiguous;
    };

    let isoform_count = count_joint_isoforms(&plan.groups);
    if isoform_count <= 1 {
        return SpectrumLocalization::Unambiguous;
    }
    if isoform_count > MAX_ISOFORMS {
        return SpectrumLocalization::Capped { isoform_count };
    }
    if input.mz.is_empty() {
        return SpectrumLocalization::MissingPeaks { isoform_count };
    }

    let peaks = window_peaks(
        input.form,
        input.charge,
        input.mz,
        input.intensity,
        tolerance,
    );

    // Enumerate and score target isoforms.
    let mut placements = Vec::with_capacity(isoform_count);
    enumerate_placements(&plan.groups, 0, &mut Vec::new(), &mut placements);
    if placements.len() <= 1 {
        return SpectrumLocalization::Unambiguous;
    }

    let Some(isoforms) = build_isoforms(
        &input.form.sequence,
        &plan.fixed_modifications,
        &plan.groups,
        &placements,
        input.charge,
    ) else {
        return SpectrumLocalization::Unscorable { isoform_count };
    };

    let evaluation = evaluate_isoform_pool(&isoforms, &plan.groups, &residues, &peaks, tolerance);

    // Evaluate decoy amino acid sites for each ambiguous modification group.
    let mut sites = evaluation.sites;
    let mut false_target_sites = 0u32;
    for (group_idx, group) in plan.groups.iter().enumerate() {
        let m_t = group.target_positions.len();
        let k = group.templates.len();
        if m_t <= k {
            continue;
        }
        // Explicit decoy residues (the search allowed the modification on
        // the decoy residue) compete for the k modifications but are no
        // false *target* opportunities.
        let real_targets = m_t - group.explicit_decoy_positions.len();
        false_target_sites += u32::try_from(real_targets.saturating_sub(k)).unwrap_or(0);

        if group.explicit_decoy_positions.is_empty() {
            // Find the lowest-scoring non-assigned target position in this group
            // so replacing it with a decoy residue keeps the isoform pool size equal.
            let drop_pos = evaluation
                .group_lowest_unassigned
                .get(group_idx)
                .copied()
                .flatten();
            let Some(drop_pos) = drop_pos else {
                continue;
            };
            for &decoy_pos in &group.synthetic_decoy_positions {
                let mut decoy_groups = plan.groups.clone();
                let target_positions = &mut decoy_groups[group_idx].target_positions;
                if let Some(slot) = target_positions.iter_mut().find(|&&mut p| p == drop_pos) {
                    *slot = decoy_pos;
                    target_positions.sort_unstable();
                } else {
                    continue;
                }
                let mut decoy_placements = Vec::with_capacity(isoform_count);
                enumerate_placements(&decoy_groups, 0, &mut Vec::new(), &mut decoy_placements);
                let Some(decoy_isoforms) = build_isoforms(
                    &input.form.sequence,
                    &plan.fixed_modifications,
                    &decoy_groups,
                    &decoy_placements,
                    input.charge,
                ) else {
                    continue;
                };
                let decoy_eval = evaluate_isoform_pool(
                    &decoy_isoforms,
                    &decoy_groups,
                    &residues,
                    &peaks,
                    tolerance,
                );
                if let Some(mut decoy_site) = decoy_eval.sites.into_iter().find(|s| {
                    s.position == decoy_pos
                        && (s.mass_delta - group.snapped_delta).abs() <= unimod::SNAP_TOLERANCE
                }) {
                    decoy_site.is_decoy = true;
                    decoy_site.assigned = false;
                    sites.push(decoy_site);
                }
            }
        }
    }

    sites.sort_by(|a, b| {
        a.is_decoy
            .cmp(&b.is_decoy)
            .then(a.position.cmp(&b.position))
            .then(a.mass_delta.total_cmp(&b.mass_delta))
    });

    SpectrumLocalization::Scored(ScoredLocalization {
        spectrum_id: input.spectrum_id,
        isoform_count,
        localized_peptide: evaluation.localized_peptide,
        best_ascore: evaluation.best_ascore,
        min_probability: evaluation.min_probability,
        false_target_sites,
        sites,
    })
}

#[derive(Debug, Clone)]
struct ModGroup {
    snapped_delta: f64,
    modification_name: Option<String>,
    /// Modification templates (count = `k`).
    templates: Vec<Modification>,
    /// 1-based candidate target positions in `sequence` (sorted ascending).
    target_positions: Vec<u32>,
    /// Positions in `target_positions` that were originally modified in `form`.
    reported_positions: Vec<u32>,
    /// Positions in `target_positions` that are already decoy residues (when a
    /// search explicitly allowed the decoy residue as a variable modification).
    explicit_decoy_positions: Vec<u32>,
    /// Unmodified decoy residue positions in `sequence` for synthetic FLR evaluation.
    synthetic_decoy_positions: Vec<u32>,
}

#[derive(Debug, Clone)]
struct ModificationPlan {
    fixed_modifications: Vec<Modification>,
    groups: Vec<ModGroup>,
}

fn build_modification_plan(form: &PeptideForm, residues: &[char]) -> Option<ModificationPlan> {
    let mut fixed_modifications = Vec::new();
    let mut fixed_positions = Vec::new();
    let mut variable_mods: Vec<(u32, Modification)> = Vec::new();

    for modification in &form.modifications {
        match modification.site {
            ModSite::NTerm | ModSite::CTerm => {
                fixed_modifications.push(modification.clone());
            }
            ModSite::Residue(pos) => {
                let idx = usize::try_from(pos).ok()?.checked_sub(1)?;
                let residue = *residues.get(idx)?;
                if is_fixed_residue_modification(modification, pos, residue) {
                    fixed_modifications.push(modification.clone());
                    fixed_positions.push(pos);
                } else {
                    variable_mods.push((pos, modification.clone()));
                }
            }
        }
    }

    if variable_mods.is_empty() {
        return None;
    }

    // Group variable modifications by snapped mass delta.
    let mut groups: Vec<ModGroup> = Vec::new();
    for (pos, modification) in variable_mods {
        let snapped = snap_or_round(modification.mass_delta);
        if let Some(group) = groups
            .iter_mut()
            .find(|g| (g.snapped_delta - snapped).abs() <= unimod::SNAP_TOLERANCE)
        {
            if group.modification_name.is_none() {
                group.modification_name = canonical_mod_name(&modification);
            }
            group.templates.push(modification);
            group.reported_positions.push(pos);
        } else {
            groups.push(ModGroup {
                snapped_delta: snapped,
                modification_name: canonical_mod_name(&modification),
                templates: vec![modification],
                target_positions: Vec::new(),
                reported_positions: vec![pos],
                explicit_decoy_positions: Vec::new(),
                synthetic_decoy_positions: Vec::new(),
            });
        }
    }

    let mut any_ambiguous = false;
    for group in &mut groups {
        group.reported_positions.sort_unstable();
        group.reported_positions.dedup();
        let canonical_residues =
            canonical_target_residues(group.snapped_delta, group.modification_name.as_deref());
        let decoy_char = choose_decoy_residue(canonical_residues);

        let mut target_positions = Vec::new();
        let mut explicit_decoy_positions = Vec::new();
        let mut synthetic_decoy_positions = Vec::new();

        for (idx, &residue) in residues.iter().enumerate() {
            let pos = u32::try_from(idx + 1).expect("short peptide");
            if fixed_positions.contains(&pos) && !group.reported_positions.contains(&pos) {
                continue;
            }
            let is_reported = group.reported_positions.contains(&pos);
            let is_canonical = match canonical_residues {
                Some(list) => list.contains(&residue),
                None => group
                    .reported_positions
                    .iter()
                    .any(|&rp| residues.get(rp as usize - 1) == Some(&residue)),
            };

            if is_canonical || is_reported {
                target_positions.push(pos);
                if residue == decoy_char
                    && canonical_residues.is_some_and(|c| !c.contains(&residue))
                {
                    explicit_decoy_positions.push(pos);
                }
            } else if residue == decoy_char {
                synthetic_decoy_positions.push(pos);
            }
        }

        // If the search itself placed a modification on the decoy residue `decoy_char`,
        // any other occurrence of `decoy_char` in the peptide is also an explicit decoy site.
        if !explicit_decoy_positions.is_empty() {
            for &pos in &synthetic_decoy_positions {
                if !target_positions.contains(&pos) {
                    target_positions.push(pos);
                }
                if !explicit_decoy_positions.contains(&pos) {
                    explicit_decoy_positions.push(pos);
                }
            }
            target_positions.sort_unstable();
            explicit_decoy_positions.sort_unstable();
            synthetic_decoy_positions.clear();
        }

        if target_positions.len() > group.templates.len() {
            any_ambiguous = true;
        }
        group.target_positions = target_positions;
        group.explicit_decoy_positions = explicit_decoy_positions;
        group.synthetic_decoy_positions = synthetic_decoy_positions;
    }

    if !any_ambiguous {
        return None;
    }

    Some(ModificationPlan {
        fixed_modifications,
        groups,
    })
}

fn is_fixed_residue_modification(modification: &Modification, pos: u32, residue: char) -> bool {
    let known = modification
        .name
        .as_deref()
        .and_then(unimod::by_name)
        .or_else(|| unimod::by_mass(modification.mass_delta));
    let Some(known) = known else {
        return false;
    };
    match known.name {
        // Alkylation of Cys is a fixed modification in standard proteomics workflows.
        "Carbamidomethyl" | "Methylthio" => residue == 'C',
        // Isobaric labeling reagents are fixed modifications on K (or N-term).
        "iTRAQ4plex" | "iTRAQ8plex" | "TMT6plex" | "TMTpro" => true,
        // Pyro-Glu cyclization requires the free N-terminal amine at position 1.
        "Glu->pyro-Glu" | "Gln->pyro-Glu" => pos == 1,
        _ => false,
    }
}

fn canonical_mod_name(modification: &Modification) -> Option<String> {
    if let Some(name) = &modification.name {
        if let Some(known) = unimod::by_name(name) {
            return Some(known.name.to_string());
        }
        return Some(name.clone());
    }
    if let Some(acc) = &modification.accession
        && let Some(delta) = unimod::mass_delta(acc)
        && let Some(known) = unimod::by_mass(delta)
    {
        return Some(known.name.to_string());
    }
    unimod::by_mass(modification.mass_delta).map(|k| k.name.to_string())
}

fn canonical_target_residues(snapped_delta: f64, name: Option<&str>) -> Option<&'static [char]> {
    let known_name = name
        .and_then(unimod::by_name)
        .or_else(|| unimod::by_mass(snapped_delta))
        .map(|k| k.name)
        .or(name);
    match known_name {
        Some("Phospho") => Some(&['S', 'T', 'Y']),
        Some("Oxidation") => Some(&['M']),
        Some("Acetyl" | "Carbamyl" | "GG") => Some(&['K']),
        Some("Deamidated") => Some(&['N', 'Q']),
        Some("Methyl" | "Dimethyl") => Some(&['K', 'R']),
        _ => None,
    }
}

fn choose_decoy_residue(canonical: Option<&[char]>) -> char {
    if canonical.is_some_and(|list| list.contains(&'A')) {
        'L'
    } else {
        'A'
    }
}

fn snap_or_round(delta: f64) -> f64 {
    if let Some(snapped) = unimod::snap(delta) {
        snapped
    } else {
        let factor = 10f64.powi(DELTA_DECIMALS as i32);
        (delta * factor).round() / factor
    }
}

fn count_joint_isoforms(groups: &[ModGroup]) -> usize {
    fn recurse(groups: &[ModGroup], idx: usize, used: &[u32]) -> usize {
        if idx == groups.len() {
            return 1;
        }
        let group = &groups[idx];
        let k = group.templates.len();
        let available: Vec<u32> = group
            .target_positions
            .iter()
            .copied()
            .filter(|p| !used.contains(p))
            .collect();
        if available.len() < k {
            return 0;
        }
        if idx + 1 == groups.len() {
            return binomial_coeff(available.len(), k);
        }
        let mut combos = Vec::new();
        choose_positions(&available, k, 0, &mut Vec::new(), &mut combos);
        let mut total = 0usize;
        for combo in combos {
            let mut next_used = used.to_vec();
            next_used.extend_from_slice(&combo);
            total = total.saturating_add(recurse(groups, idx + 1, &next_used));
            if total > MAX_ISOFORMS {
                return total;
            }
        }
        total
    }
    recurse(groups, 0, &[])
}

fn binomial_coeff(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1usize;
    for i in 0..k {
        result = result.saturating_mul(n - i) / (i + 1);
    }
    result
}

fn choose_positions(
    available: &[u32],
    k: usize,
    start: usize,
    current: &mut Vec<u32>,
    out: &mut Vec<Vec<u32>>,
) {
    if current.len() == k {
        out.push(current.clone());
        return;
    }
    if available.len() - start < k - current.len() {
        return;
    }
    for i in start..available.len() {
        current.push(available[i]);
        choose_positions(available, k, i + 1, current, out);
        current.pop();
        if out.len() > MAX_ISOFORMS {
            return;
        }
    }
}

fn enumerate_placements(
    groups: &[ModGroup],
    idx: usize,
    current: &mut Vec<Vec<u32>>,
    out: &mut Vec<Vec<Vec<u32>>>,
) {
    if idx == groups.len() {
        out.push(current.clone());
        return;
    }
    let group = &groups[idx];
    let k = group.templates.len();
    let available: Vec<u32> = group
        .target_positions
        .iter()
        .copied()
        .filter(|p| !current.iter().flatten().any(|u| u == p))
        .collect();
    let mut combos = Vec::new();
    choose_positions(&available, k, 0, &mut Vec::new(), &mut combos);
    for combo in combos {
        current.push(combo);
        enumerate_placements(groups, idx + 1, current, out);
        current.pop();
        if out.len() > MAX_ISOFORMS {
            return;
        }
    }
}

#[derive(Debug, Clone)]
struct CandidateIsoform {
    form: PeptideForm,
    ions: Vec<FragmentIon>,
    /// Per-group chosen 1-based positions.
    group_positions: Vec<Vec<u32>>,
    /// True if this isoform matches the positions originally reported on `form`.
    is_reported: bool,
}

fn build_isoforms(
    sequence: &str,
    fixed_modifications: &[Modification],
    groups: &[ModGroup],
    placements: &[Vec<Vec<u32>>],
    charge: Option<u8>,
) -> Option<Vec<CandidateIsoform>> {
    let mut isoforms = Vec::with_capacity(placements.len());
    for placement in placements {
        let mut mods = fixed_modifications.to_vec();
        let mut is_reported = true;
        for (group, chosen) in groups.iter().zip(placement.iter()) {
            if chosen != &group.reported_positions {
                is_reported = false;
            }
            for (template, &pos) in group.templates.iter().zip(chosen.iter()) {
                let mut m = template.clone();
                m.site = ModSite::Residue(pos);
                mods.push(m);
            }
        }
        let form = PeptideForm::new(sequence, mods).ok()?;
        let ions = fragment_ions(&form, charge)?;
        isoforms.push(CandidateIsoform {
            form,
            ions,
            group_positions: placement.clone(),
            is_reported,
        });
    }
    Some(isoforms)
}

#[derive(Debug, Clone, Copy)]
struct WindowedPeak {
    mz: f64,
    intensity: f32,
    /// 1-based rank within its 100 m/z window (`1..=MAX_DEPTH`).
    rank: u8,
}

/// Filter out precursor and precursor neutral-loss peaks, then rank remaining
/// peaks within each 100 m/z window up to depth [`MAX_DEPTH`].
fn window_peaks(
    form: &PeptideForm,
    charge: Option<u8>,
    mz: &[f64],
    intensity: &[f32],
    tolerance: Tolerance,
) -> Vec<WindowedPeak> {
    let exclusions = precursor_exclusion_mzs(form, charge);
    let mut by_window: std::collections::BTreeMap<i64, Vec<(f64, f32)>> =
        std::collections::BTreeMap::new();

    for (&peak_mz, &peak_int) in mz.iter().zip(intensity.iter()) {
        if !peak_mz.is_finite() || peak_mz <= 0.0 || !peak_int.is_finite() || peak_int <= 0.0 {
            continue;
        }
        if exclusions
            .iter()
            .any(|&excl| (peak_mz - excl).abs() <= tolerance.window(excl).max(0.05))
        {
            continue;
        }
        let window_idx = (peak_mz / 100.0).floor() as i64;
        by_window
            .entry(window_idx)
            .or_default()
            .push((peak_mz, peak_int));
    }

    let mut kept = Vec::new();
    for (_, mut window_peaks) in by_window {
        window_peaks.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.total_cmp(&b.0)));
        for (idx, (peak_mz, peak_int)) in window_peaks.into_iter().take(MAX_DEPTH).enumerate() {
            kept.push(WindowedPeak {
                mz: peak_mz,
                intensity: peak_int,
                rank: (idx + 1) as u8,
            });
        }
    }
    kept.sort_by(|a, b| a.mz.total_cmp(&b.mz));
    kept
}

fn precursor_exclusion_mzs(form: &PeptideForm, charge: Option<u8>) -> Vec<f64> {
    let Some(z) = charge.filter(|&z| z >= 1) else {
        return Vec::new();
    };
    let Some(unmodified) = peptide_mass(&form.sequence) else {
        return Vec::new();
    };
    let mod_sum: f64 = form.modifications.iter().map(|m| m.mass_delta).sum();
    let neutral = unmodified + mod_sum;
    let zf = f64::from(z);
    [
        0.0,
        AMMONIA_MONO,
        WATER_MONO,
        HPO3_MONO,
        H3PO4_MONO,
        2.0 * H3PO4_MONO,
    ]
    .into_iter()
    .filter_map(|loss| {
        let mass = neutral - loss;
        (mass > 0.0).then_some((mass + zf * PROTON) / zf)
    })
    .collect()
}

/// For each ion of each isoform, find the best (lowest window rank, highest
/// intensity) matching peak within `tolerance`.
fn match_isoform_ions(
    isoforms: &[CandidateIsoform],
    peaks: &[WindowedPeak],
    tolerance: Tolerance,
) -> Vec<Vec<Option<u8>>> {
    isoforms
        .iter()
        .map(|iso| {
            iso.ions
                .iter()
                .map(|ion| match_ion_rank(ion.mz, peaks, tolerance))
                .collect()
        })
        .collect()
}

fn match_ion_rank(ion_mz: f64, peaks: &[WindowedPeak], tolerance: Tolerance) -> Option<u8> {
    let win = tolerance.window(ion_mz);
    let first = peaks.partition_point(|p| p.mz < ion_mz - win);
    let last = peaks.partition_point(|p| p.mz <= ion_mz + win);
    peaks[first..last]
        .iter()
        .max_by(|a, b| {
            a.intensity
                .total_cmp(&b.intensity)
                .then_with(|| (b.mz - ion_mz).abs().total_cmp(&(a.mz - ion_mz).abs()))
        })
        .map(|p| p.rank)
}

#[derive(Debug, Clone, Copy)]
struct PairwiseComparison {
    score_a: f64,
    score_b: f64,
    matched_a_depth6: u32,
    total_variant: u32,
}

/// Compare two isoforms `a` and `b` specifically on the theoretical ions whose
/// m/z differs between `a` and `b`.
fn compare_isoform_pair(
    a: &CandidateIsoform,
    ranks_a: &[Option<u8>],
    b: &CandidateIsoform,
    ranks_b: &[Option<u8>],
) -> PairwiseComparison {
    let mut counts_a = [0usize; MAX_DEPTH];
    let mut counts_b = [0usize; MAX_DEPTH];
    let mut total_variant = 0usize;

    for i in 0..a.ions.len() {
        if (a.ions[i].mz - b.ions[i].mz).abs() <= 1e-6 {
            continue;
        }
        total_variant += 1;
        if let Some(r) = ranks_a[i] {
            for count in &mut counts_a[usize::from(r - 1)..MAX_DEPTH] {
                *count += 1;
            }
        }
        if let Some(r) = ranks_b[i] {
            for count in &mut counts_b[usize::from(r - 1)..MAX_DEPTH] {
                *count += 1;
            }
        }
    }

    if total_variant == 0 {
        return PairwiseComparison {
            score_a: 0.0,
            score_b: 0.0,
            matched_a_depth6: 0,
            total_variant: 0,
        };
    }

    let score_a = weighted_binomial_score(&counts_a, total_variant);
    let score_b = weighted_binomial_score(&counts_b, total_variant);
    PairwiseComparison {
        score_a,
        score_b,
        matched_a_depth6: counts_a[REPORT_DEPTH - 1] as u32,
        total_variant: total_variant as u32,
    }
}

/// Score an isoform across the union of all variant ions in the pool (used for
/// initial ranking before pairwise refinement).
fn score_isoform_on_mask(ranks: &[Option<u8>], variant_mask: &[bool], total_variant: usize) -> f64 {
    if total_variant == 0 {
        return 0.0;
    }
    let mut counts = [0usize; MAX_DEPTH];
    for (&rank, &is_var) in ranks.iter().zip(variant_mask.iter()) {
        if is_var && let Some(r) = rank {
            for count in &mut counts[usize::from(r - 1)..MAX_DEPTH] {
                *count += 1;
            }
        }
    }
    weighted_binomial_score(&counts, total_variant)
}

fn weighted_binomial_score(counts: &[usize; MAX_DEPTH], n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut weighted_sum = 0.0;
    for (idx, (&k, &weight)) in counts.iter().zip(DEPTH_WEIGHTS.iter()).enumerate() {
        if k > 0 {
            let p = 0.01 * ((idx + 1) as f64);
            let tail_p = cumulative_binomial_tail(k, n, p);
            let depth_score = -10.0 * tail_p.log10();
            weighted_sum += weight * depth_score;
        }
    }
    weighted_sum / DEPTH_WEIGHT_SUM
}

/// Upper tail $P(X \ge k \mid n, p) = \sum_{j=k}^n \binom{n}{j} p^j (1 - p)^{n - j}$,
/// computed in log-space when $n$ or $k$ is large to prevent underflow.
pub fn cumulative_binomial_tail(k: usize, n: usize, p: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if k > n || p <= 0.0 {
        return f64::MIN_POSITIVE;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let ln_p = p.ln();
    let ln_q = (1.0 - p).ln();
    let mut terms = Vec::with_capacity(n - k + 1);
    let mut max_ln = f64::NEG_INFINITY;
    for j in k..=n {
        let ln_term = ln_choose(n, j) + (j as f64) * ln_p + ((n - j) as f64) * ln_q;
        if ln_term > max_ln {
            max_ln = ln_term;
        }
        terms.push(ln_term);
    }
    let sum_exp: f64 = terms.into_iter().map(|t| (t - max_ln).exp()).sum();
    (max_ln + sum_exp.ln()).exp().clamp(f64::MIN_POSITIVE, 1.0)
}

fn ln_choose(n: usize, k: usize) -> f64 {
    if k == 0 || k == n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut acc = 0.0;
    for i in 0..k {
        acc += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
    }
    acc
}

struct PoolEvaluation {
    localized_peptide: String,
    best_ascore: f64,
    min_probability: f64,
    sites: Vec<SiteScore>,
    group_lowest_unassigned: Vec<Option<u32>>,
}

fn evaluate_isoform_pool(
    isoforms: &[CandidateIsoform],
    groups: &[ModGroup],
    residues: &[char],
    peaks: &[WindowedPeak],
    tolerance: Tolerance,
) -> PoolEvaluation {
    let ion_count = isoforms[0].ions.len();
    let mut variant_mask = vec![false; ion_count];
    let mut total_variant = 0usize;
    for (q, is_variant) in variant_mask.iter_mut().enumerate() {
        let m0 = isoforms[0].ions[q].mz;
        if isoforms[1..]
            .iter()
            .any(|iso| (iso.ions[q].mz - m0).abs() > 1e-6)
        {
            *is_variant = true;
            total_variant += 1;
        }
    }

    let matched_ranks = match_isoform_ions(isoforms, peaks, tolerance);
    let global_scores: Vec<f64> = matched_ranks
        .iter()
        .map(|ranks| score_isoform_on_mask(ranks, &variant_mask, total_variant))
        .collect();

    // Pick initial winning isoform (highest global score; ties keep reported
    // placement first, then lexicographical positions).
    let mut best_idx = 0usize;
    for idx in 1..isoforms.len() {
        if isoform_better(
            global_scores[idx],
            &isoforms[idx],
            global_scores[best_idx],
            &isoforms[best_idx],
        ) {
            best_idx = idx;
        }
    }

    // Refine `best_idx` so no rival beats it in direct pairwise comparison on
    // their pairwise site-determining ions.
    for _ in 0..3 {
        let mut challenger = None;
        let mut best_margin = 0.0;
        for idx in 0..isoforms.len() {
            if idx == best_idx {
                continue;
            }
            let cmp = compare_isoform_pair(
                &isoforms[best_idx],
                &matched_ranks[best_idx],
                &isoforms[idx],
                &matched_ranks[idx],
            );
            let margin = cmp.score_b - cmp.score_a;
            if margin > best_margin + 1e-9 {
                best_margin = margin;
                challenger = Some(idx);
            }
        }
        if let Some(next) = challenger {
            best_idx = next;
        } else {
            break;
        }
    }

    // Compute pairwise score deltas relative to `best_idx`, and transform to
    // isoform posterior probabilities.
    let deltas_from_best: Vec<f64> = (0..isoforms.len())
        .map(|idx| {
            if idx == best_idx {
                0.0
            } else {
                let cmp = compare_isoform_pair(
                    &isoforms[best_idx],
                    &matched_ranks[best_idx],
                    &isoforms[idx],
                    &matched_ranks[idx],
                );
                cmp.score_a - cmp.score_b
            }
        })
        .collect();

    let min_delta = deltas_from_best.iter().copied().fold(0.0f64, f64::min);
    let weights: Vec<f64> = deltas_from_best
        .iter()
        .map(|&d| 10f64.powf(-(d - min_delta) / 10.0))
        .collect();
    let weight_sum: f64 = weights.iter().sum();
    let isoform_probs: Vec<f64> = if weight_sum > 0.0 {
        weights.into_iter().map(|w| w / weight_sum).collect()
    } else {
        vec![1.0 / (isoforms.len() as f64); isoforms.len()]
    };

    let best_isoform = &isoforms[best_idx];
    let mut sites = Vec::new();
    let mut group_lowest_unassigned = Vec::with_capacity(groups.len());

    for (group_idx, group) in groups.iter().enumerate() {
        let assigned_positions = &best_isoform.group_positions[group_idx];
        let mut lowest_unassigned: Option<(f64, u32)> = None;

        for &pos in &group.target_positions {
            let site_prob: f64 = isoforms
                .iter()
                .zip(isoform_probs.iter())
                .filter(|(iso, _)| iso.group_positions[group_idx].contains(&pos))
                .map(|(_, &prob)| prob)
                .sum::<f64>()
                .clamp(0.0, 1.0);

            // Best isoform with `pos` modified vs best isoform with `pos` unmodified.
            let best_with =
                best_isoform_matching(isoforms, &isoform_probs, &global_scores, best_idx, |iso| {
                    iso.group_positions[group_idx].contains(&pos)
                });
            let best_without =
                best_isoform_matching(isoforms, &isoform_probs, &global_scores, best_idx, |iso| {
                    !iso.group_positions[group_idx].contains(&pos)
                });

            let (ascore, matched, total) = match (best_with, best_without) {
                (Some(w_idx), Some(wo_idx)) => {
                    let cmp = compare_isoform_pair(
                        &isoforms[w_idx],
                        &matched_ranks[w_idx],
                        &isoforms[wo_idx],
                        &matched_ranks[wo_idx],
                    );
                    let ascore = (cmp.score_a - cmp.score_b).max(0.0);
                    (ascore, cmp.matched_a_depth6, cmp.total_variant)
                }
                _ => (0.0, 0, 0),
            };

            let assigned = assigned_positions.contains(&pos);
            if !assigned
                && !group.explicit_decoy_positions.contains(&pos)
                && lowest_unassigned.is_none_or(|(best_p, _)| site_prob < best_p)
            {
                lowest_unassigned = Some((site_prob, pos));
            }

            let residue = residues
                .get(usize::try_from(pos).unwrap_or(1) - 1)
                .copied()
                .unwrap_or('?');
            sites.push(SiteScore {
                position: pos,
                residue,
                mass_delta: group.snapped_delta,
                modification_name: group.modification_name.clone(),
                probability: site_prob,
                ascore,
                determining_ions_matched: matched,
                determining_ions_total: total,
                assigned,
                is_decoy: group.explicit_decoy_positions.contains(&pos),
            });
        }

        group_lowest_unassigned.push(lowest_unassigned.map(|(_, pos)| pos));
    }

    let best_ascore = sites
        .iter()
        .filter(|s| s.assigned && !s.is_decoy)
        .map(|s| s.ascore)
        .fold(0.0f64, f64::max);
    // A best isoform carrying the modification on a decoy residue localizes
    // nothing: the PSM counts as not localized.
    let min_probability = if sites.iter().any(|s| s.assigned && s.is_decoy) {
        0.0
    } else {
        sites
            .iter()
            .filter(|s| s.assigned)
            .map(|s| s.probability)
            .reduce(f64::min)
            .unwrap_or(1.0)
    };

    PoolEvaluation {
        localized_peptide: best_isoform.form.proforma(),
        best_ascore,
        min_probability,
        sites,
        group_lowest_unassigned,
    }
}

fn isoform_better(
    score_a: f64,
    iso_a: &CandidateIsoform,
    score_b: f64,
    iso_b: &CandidateIsoform,
) -> bool {
    let diff = score_a - score_b;
    if diff.abs() > 1e-9 {
        return diff > 0.0;
    }
    if iso_a.is_reported != iso_b.is_reported {
        return iso_a.is_reported;
    }
    iso_a.group_positions < iso_b.group_positions
}

fn best_isoform_matching<F>(
    isoforms: &[CandidateIsoform],
    isoform_probs: &[f64],
    global_scores: &[f64],
    best_idx: usize,
    predicate: F,
) -> Option<usize>
where
    F: Fn(&CandidateIsoform) -> bool,
{
    if predicate(&isoforms[best_idx]) {
        return Some(best_idx);
    }
    (0..isoforms.len())
        .filter(|&i| predicate(&isoforms[i]))
        .max_by(|&a, &b| {
            isoform_probs[a]
                .total_cmp(&isoform_probs[b])
                .then_with(|| global_scores[a].total_cmp(&global_scores[b]))
                .then_with(|| isoforms[a].is_reported.cmp(&isoforms[b].is_reported))
                .then_with(|| {
                    isoforms[b]
                        .group_positions
                        .cmp(&isoforms[a].group_positions)
                })
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mod_at(pos: u32, delta: f64, name: Option<&str>) -> Modification {
        Modification {
            site: ModSite::Residue(pos),
            mass_delta: delta,
            name: name.map(str::to_string),
            accession: None,
        }
    }

    #[test]
    fn unmodified_and_single_site_peptides_are_unambiguous() {
        let unmodified = PeptideForm::unmodified("PEPTIDEK");
        assert_eq!(
            localize_spectrum(
                &LocalizationInput {
                    spectrum_id: 1,
                    form: &unmodified,
                    charge: Some(2),
                    mz: &[200.0],
                    intensity: &[100.0],
                },
                Tolerance::Ppm(20.0),
            ),
            SpectrumLocalization::Unambiguous
        );

        // Only one S/T/Y in PEPTIDEK (T4), so Phospho has 1 candidate site = Unambiguous.
        let single_site =
            PeptideForm::new("PEPTIDEK", vec![mod_at(4, 79.966_331, Some("Phospho"))]).unwrap();
        assert_eq!(
            localize_spectrum(
                &LocalizationInput {
                    spectrum_id: 2,
                    form: &single_site,
                    charge: Some(2),
                    mz: &[200.0],
                    intensity: &[100.0],
                },
                Tolerance::Ppm(20.0),
            ),
            SpectrumLocalization::Unambiguous
        );
    }

    #[test]
    fn hand_derived_two_site_localization_scores_match_binomial_math() {
        // Peptide SPTK (2+ precursor -> 1+ b1..b3 and y1..y3 ions), 1 Phospho (+79.966331).
        // Candidate sites: S1 and T3.
        // Residue masses:
        //   S = 87.0320284, P = 97.0527639, T = 101.0476785, K = 128.0949630
        //   Phospho = 79.966331, H+ = 1.0072765, H2O = 18.0105647
        // Isoform 1 (S[+79.9663]PTK):
        //   b1 = 87.0320284 + 79.966331 + 1.0072765 = 168.0056359 (differs from T3: 88.0393049)
        //   b2 = 168.0056359 + 97.0527639 = 265.0583998          (differs from T3: 185.0920688)
        //   b3 = shared (holds both S1 and T3)
        //   y1 = 128.0949630 + 18.0105647 + 1.0072765 = 147.1128042 (shared: holds neither)
        //   y2 = 147.1128042 + 101.0476785 = 248.1604827         (differs from T3: 328.1268137)
        //   y3 = 248.1604827 + 97.0527639 = 345.2132466          (differs from T3: 425.1795776)
        // Total pairwise site-determining ions: n = 4 (b1, b2, y2, y3).
        let form = PeptideForm::new("SPTK", vec![mod_at(1, 79.966_331, Some("Phospho"))]).unwrap();

        // Case 1: Only shared y1 (147.1128042) present -> k = 0/4 for both isoforms.
        // P(X >= 0) = 1.0 -> score = 0.0 -> P(S1) = 0.5, P(T3) = 0.5, Ascore = 0.0.
        let SpectrumLocalization::Scored(no_determining) = localize_spectrum(
            &LocalizationInput {
                spectrum_id: 10,
                form: &form,
                charge: Some(2),
                mz: &[147.112_804_2],
                intensity: &[1000.0],
            },
            Tolerance::Ppm(20.0),
        ) else {
            panic!("expected Scored");
        };
        assert_eq!(no_determining.sites.len(), 2);
        assert!((no_determining.sites[0].probability - 0.5).abs() < 1e-6);
        assert!((no_determining.sites[1].probability - 0.5).abs() < 1e-6);
        assert_eq!(no_determining.best_ascore, 0.0);

        // Case 2: 1 site-determining ion (b1 of S1 = 168.0056359) at window rank 1.
        // A rank-1 peak counts at every depth, so k_d = 1 of n = 4 for S1 and
        // k_d = 0 (s_d = 0) for T3, whose four ions match nothing. At depth d
        // (p_d = 0.01 * d): P_d(X >= 1) = 1 - (1 - p_d)^4 and
        // s_d = -10 * log10(1 - (1 - p_d)^4), by hand (d = 1..10):
        //   14.0446, 11.0996, 9.4041, 8.2202, 7.3167, 6.5906, 5.9869, 5.4728, 5.0272, 4.6357
        // Weighted average with [0.5,0.75,1,1,1,1,0.75,0.5,0.25,0.25]/7:
        //   score(S1) = 8.0744, score(T3) = 0, so Ascore = 8.0744 and
        //   P(S1) = 1 / (1 + 10^(-0.80744)) = 0.8652.
        const ASCORE_1ION: f64 = 8.074_415;
        const PROBABILITY_1ION: f64 = 0.865_204;

        let SpectrumLocalization::Scored(one_ion) = localize_spectrum(
            &LocalizationInput {
                spectrum_id: 11,
                form: &form,
                charge: Some(2),
                mz: &[168.005_635_9],
                intensity: &[1000.0],
            },
            Tolerance::Ppm(20.0),
        ) else {
            panic!("expected Scored");
        };
        assert!((one_ion.best_ascore - ASCORE_1ION).abs() < 1e-3);
        assert!((one_ion.sites[0].probability - PROBABILITY_1ION).abs() < 1e-3);
        assert!((one_ion.sites[1].probability - (1.0 - PROBABILITY_1ION)).abs() < 1e-3);
        assert_eq!(
            (
                one_ion.sites[0].determining_ions_matched,
                one_ion.sites[0].determining_ions_total
            ),
            (1, 4)
        );

        // Case 3: 2 site-determining ions for T3 (b1 of T3 = 88.0393049 in window 0,
        // y2 of T3 = 328.1268137 in window 3, so both are rank 1 in their 100 m/z window).
        // Even though the input form reported Phospho on S1, the spectrum evidence places
        // it on T3 (`SPT[+79.9663]K`) with k = 2/4:
        //   P_d(X >= 2) = 1 - (1 - p)^4 - 4*p*(1 - p)^3, k_d = 2 at every depth (both
        //   peaks rank 1), and k_d = 0 for S1. By hand, s_d (d = 1..10):
        //   32.2766, 26.3144, 22.8513, 20.4116, 18.5329, 17.0091, 15.7303, 14.6310, 13.6689, 12.8150
        //   weighted: Ascore = 20.0590, P(T3) = 1 / (1 + 10^(-2.00590)) = 0.99023 -> 99.0%.
        const ASCORE_2ION: f64 = 20.059_036;
        const PROBABILITY_2ION: f64 = 0.990_231;

        let SpectrumLocalization::Scored(two_ions_t3) = localize_spectrum(
            &LocalizationInput {
                spectrum_id: 12,
                form: &form,
                charge: Some(2),
                mz: &[88.039_304_9, 328.126_813_7],
                intensity: &[1000.0, 900.0],
            },
            Tolerance::Ppm(20.0),
        ) else {
            panic!("expected Scored");
        };
        assert_eq!(two_ions_t3.localized_peptide, "SPT[+79.9663]K");
        assert!(!two_ions_t3.sites[0].assigned);
        assert!(two_ions_t3.sites[1].assigned);
        assert!((two_ions_t3.best_ascore - ASCORE_2ION).abs() < 1e-3);
        assert!((two_ions_t3.sites[1].probability - PROBABILITY_2ION).abs() < 1e-3);
        assert_eq!(
            two_ions_t3.format_site_probabilities(),
            "S1: 1.0%, T3: 99.0%"
        );
    }

    #[test]
    fn precursor_neutral_loss_peaks_do_not_crowd_out_fragment_window() {
        // SPTK + Phospho (2+): neutral mass = 87.0320284 + 97.0527639 + 101.0476785
        // + 128.0949630 + 18.0105647 + 79.966331 = 511.2043295 Da.
        // Precursor - H3PO4 (97.976896 Da) at 2+:
        //   (511.2043295 - 97.976896 + 2 * 1.0072765) / 2 = 207.6209933 m/z (in window [200, 300)).
        // Fill window [200, 300) with 10 giant peaks at the precursor - H3PO4 m/z, plus 9
        // decoy peaks at 210..218, and put the real b2 ion of S1 (265.0583998) at 10th rank
        // among real fragment peaks. Because the precursor - H3PO4 peak is excluded before
        // windowing, b2 stays within depth 10 and is matched.
        let form = PeptideForm::new("SPTK", vec![mod_at(1, 79.966_331, Some("Phospho"))]).unwrap();
        let prec_loss_mz = 207.620_993_3;
        let b2_s1_mz = 265.058_399_8;
        let mut mz = vec![
            prec_loss_mz,
            210.0,
            211.0,
            212.0,
            213.0,
            214.0,
            215.0,
            216.0,
            217.0,
            218.0,
            b2_s1_mz,
        ];
        let mut intensity = vec![
            1_000_000.0,
            900.0,
            800.0,
            700.0,
            600.0,
            500.0,
            400.0,
            300.0,
            200.0,
            100.0,
            50.0,
        ];
        // Sort ascending by m/z as required.
        let mut pairs: Vec<(f64, f32)> = mz.into_iter().zip(intensity).collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        mz = pairs.iter().map(|p| p.0).collect();
        intensity = pairs.iter().map(|p| p.1).collect();

        let SpectrumLocalization::Scored(scored) = localize_spectrum(
            &LocalizationInput {
                spectrum_id: 20,
                form: &form,
                charge: Some(2),
                mz: &mz,
                intensity: &intensity,
            },
            Tolerance::Ppm(20.0),
        ) else {
            panic!("expected Scored");
        };
        assert!(
            scored.best_ascore > 0.0,
            "b2 should reach top-10 after precursor -H3PO4 exclusion"
        );
    }

    #[test]
    fn missing_peaks_and_isoform_cap_are_counted_never_fabricated() {
        let form = PeptideForm::new("SPTK", vec![mod_at(1, 79.966_331, Some("Phospho"))]).unwrap();
        assert_eq!(
            localize_spectrum(
                &LocalizationInput {
                    spectrum_id: 30,
                    form: &form,
                    charge: Some(2),
                    mz: &[],
                    intensity: &[],
                },
                Tolerance::Ppm(20.0),
            ),
            SpectrumLocalization::MissingPeaks { isoform_count: 2 }
        );

        // 12 Serines with 5 Phospho modifications: C(12, 5) = 792 > MAX_ISOFORMS (256).
        let many_sites = PeptideForm::new(
            "SSSSSSSSSSSSK",
            vec![
                mod_at(1, 79.966_331, Some("Phospho")),
                mod_at(2, 79.966_331, Some("Phospho")),
                mod_at(3, 79.966_331, Some("Phospho")),
                mod_at(4, 79.966_331, Some("Phospho")),
                mod_at(5, 79.966_331, Some("Phospho")),
            ],
        )
        .unwrap();
        assert_eq!(
            localize_spectrum(
                &LocalizationInput {
                    spectrum_id: 31,
                    form: &many_sites,
                    charge: Some(2),
                    mz: &[200.0],
                    intensity: &[100.0],
                },
                Tolerance::Ppm(20.0),
            ),
            SpectrumLocalization::Capped { isoform_count: 792 }
        );
    }

    #[test]
    fn synthetic_site_probabilities_and_flr_are_calibrated() {
        // Deterministic PRNG (xorshift64) generating 600 synthetic phosphopeptide
        // spectra with 2 or 3 S/T/Y candidate sites + 1 Alanine decoy site.
        // Each spectrum has random fragmentation coverage of the true isoform plus
        // random noise peaks, spanning the full confidence range from 50% to 99.9%.
        fn next_u64(rng: &mut u64) -> u64 {
            *rng ^= *rng << 13;
            *rng ^= *rng >> 7;
            *rng ^= *rng << 17;
            *rng
        }
        fn next_f64(rng: &mut u64) -> f64 {
            (next_u64(rng) as f64) / (u64::MAX as f64)
        }
        let mut rng = 0x5EED_CAFE_1234_5678u64;

        let templates = [
            ("ASGPTEYK", [2u32, 5, 7]),
            ("SAGPTVSK", [1, 5, 7]),
            ("GAPSTLDK", [4, 5, 0]),
            ("ALSTGPEK", [3, 4, 0]),
            ("SVAGTYPK", [1, 5, 6]),
        ];

        // Calibration bins for target site probabilities: [0.0, 0.25), [0.25, 0.60),
        // [0.60, 0.90), [0.90, 1.01].
        let mut bin_pred_sum = [0.0f64; 4];
        let mut bin_true_sum = [0usize; 4];
        let mut bin_count = [0usize; 4];
        let mut results = Vec::new();
        let mut true_passing_95 = 0usize;
        let mut false_passing_95 = 0usize;

        for spec_id in 0..600i64 {
            let (seq, sites) = templates[(spec_id as usize) % templates.len()];
            let valid_sites: Vec<u32> = sites.into_iter().filter(|&p| p > 0).collect();
            let true_pos = valid_sites[(next_u64(&mut rng) as usize) % valid_sites.len()];
            let reported_pos = valid_sites[(next_u64(&mut rng) as usize) % valid_sites.len()];

            let true_form =
                PeptideForm::new(seq, vec![mod_at(true_pos, 79.966_331, Some("Phospho"))]).unwrap();
            let reported_form =
                PeptideForm::new(seq, vec![mod_at(reported_pos, 79.966_331, Some("Phospho"))])
                    .unwrap();

            let ions = fragment_ions(&true_form, Some(2)).unwrap();
            // Vary ion survival probability from 0.05 to 0.70 across spectra.
            let survival = 0.05 + 0.65 * next_f64(&mut rng);
            let mut peaks: Vec<(f64, f32)> = Vec::new();
            for ion in &ions {
                if next_f64(&mut rng) < survival {
                    peaks.push((ion.mz, 200.0 + 800.0 * next_f64(&mut rng) as f32));
                }
            }
            // Noise: each distinct ion of the candidate isoforms and of the Ala
            // isoform, at 4%, so wrong sites and the decoy get chance matches
            // equally often. (Drawn per isoform instead, ions shared by several
            // wrong isoforms got several chances and the decoy estimate ran low:
            // see tests/localization_flr.rs.)
            let ala_pos = seq.find('A').unwrap() as u32 + 1;
            let mut noise_mz: Vec<f64> = Vec::new();
            for &alt_pos in valid_sites.iter().chain(std::iter::once(&ala_pos)) {
                let alt_form =
                    PeptideForm::new(seq, vec![mod_at(alt_pos, 79.966_331, Some("Phospho"))])
                        .unwrap();
                noise_mz.extend(
                    fragment_ions(&alt_form, Some(2))
                        .unwrap()
                        .iter()
                        .map(|i| i.mz),
                );
            }
            noise_mz.sort_by(f64::total_cmp);
            noise_mz.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
            for mz in noise_mz {
                if next_f64(&mut rng) < 0.04 {
                    peaks.push((mz, 50.0 + 300.0 * next_f64(&mut rng) as f32));
                }
            }

            if peaks.is_empty() {
                peaks.push((150.0, 100.0));
            }
            peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
            let intensity: Vec<f32> = peaks.iter().map(|p| p.1).collect();

            let loc = localize_spectrum(
                &LocalizationInput {
                    spectrum_id: spec_id,
                    form: &reported_form,
                    charge: Some(2),
                    mz: &mz,
                    intensity: &intensity,
                },
                Tolerance::Ppm(20.0),
            );
            if let SpectrumLocalization::Scored(ref scored) = loc {
                for site in scored.target_sites() {
                    let is_true = site.position == true_pos;
                    let bin = if site.probability < 0.25 {
                        0
                    } else if site.probability < 0.60 {
                        1
                    } else if site.probability < 0.90 {
                        2
                    } else {
                        3
                    };
                    bin_pred_sum[bin] += site.probability;
                    if is_true {
                        bin_true_sum[bin] += 1;
                    }
                    bin_count[bin] += 1;

                    if site.probability >= 0.95 {
                        if is_true {
                            true_passing_95 += 1;
                        } else {
                            false_passing_95 += 1;
                        }
                    }
                }
            }
            results.push(loc);
        }

        // 1. Empirical true fraction increases monotonically across the 4 probability bins
        // and stays close to the mean predicted probability in each bin.
        let mut prev_frac = -1.0;
        for b in 0..4 {
            assert!(bin_count[b] >= 20, "bin {b} too small: {}", bin_count[b]);
            let pred_mean = bin_pred_sum[b] / (bin_count[b] as f64);
            let true_frac = (bin_true_sum[b] as f64) / (bin_count[b] as f64);
            assert!(
                true_frac > prev_frac,
                "empirical accuracy must increase with probability bin ({prev_frac} -> {true_frac})"
            );
            assert!(
                (true_frac - pred_mean).abs() < 0.12,
                "bin {b}: empirical {true_frac:.3} vs predicted {pred_mean:.3}"
            );
            prev_frac = true_frac;
        }

        // 2. At the 95% probability cut, true empirical FLR is <= 5% and estimated FLR
        // from decoy Alanine sites is in the same low regime.
        let empirical_flr_95 =
            (false_passing_95 as f64) / ((true_passing_95 + false_passing_95) as f64);
        let summary = summarize_localizations(&results, 0.95);
        let est_flr_95 = summary.estimated_flr.expect("FLR estimated");
        assert!(
            empirical_flr_95 <= 0.05,
            "empirical FLR at 95% cut too high: {empirical_flr_95:.4}"
        );
        assert!(
            (est_flr_95 - empirical_flr_95).abs() < 0.03,
            "estimated FLR {est_flr_95:.4} vs empirical FLR {empirical_flr_95:.4}"
        );
    }
}
