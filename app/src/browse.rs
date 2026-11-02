//! View models for the browsing views: loading a project's protein groups and
//! drilling down (group -> peptides -> spectra -> one spectrum), the
//! annotated spectrum plot, and the coverage map. Pure data and store calls;
//! no GPUI, so it is unit-tested and runs off the UI thread.

use std::path::Path;

use tuffold_core::coverage::{SequenceCoverage, covered_mask};
use tuffold_core::fragments::{IonSeries, Tolerance, annotate, fragment_ions};
use tuffold_core::peptide::DELTA_DECIMALS;
use tuffold_core::{ModSite, unimod};
use tuffold_store::{GroupDetail, Project, ProteinList, SpectrumDetail, SpectrumPage, StoredPeaks};

/// Residues per line of the coverage map.
pub const COVERAGE_LINE: usize = 60;
/// Fragment tolerances offered for annotation.
pub const TOLERANCES: [Tolerance; 2] = [Tolerance::Ppm(20.0), Tolerance::Da(0.5)];

/// Everything the browser shows at once: the list and the current drill-down.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct BrowserData {
    pub list: ProteinList,
    pub drill: Drill,
}

/// The selected group, its peptides' spectra, and the selected spectrum.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Drill {
    pub group: Option<GroupDetail>,
    pub peptide_id: Option<i64>,
    pub spectra: SpectrumPage,
    pub spectrum: Option<SpectrumDetail>,
}

/// At most this many of a peptide's spectra are listed (the best); the
/// heading says how many there are in all.
pub const SPECTRA_SHOWN: usize = 500;

fn open(path: &Path) -> Result<Project, String> {
    Project::open(path).map_err(|error| error.to_string())
}

/// The protein list and the drill-down of `prefer` when it is still listed,
/// else of the first group.
pub fn load(path: &Path, prefer: Option<i64>) -> Result<BrowserData, String> {
    let project = open(path)?;
    let list = project.protein_list().map_err(|e| e.to_string())?;
    let chosen = prefer
        .filter(|id| list.rows.iter().any(|row| row.group_id == *id))
        .or_else(|| list.rows.first().map(|row| row.group_id));
    let drill = match chosen {
        Some(group_id) => drill_group(&project, group_id)?,
        None => Drill::default(),
    };
    Ok(BrowserData { list, drill })
}

/// A group, with its first peptide's spectra and the best of them.
pub fn load_group(path: &Path, group_id: i64) -> Result<Drill, String> {
    drill_group(&open(path)?, group_id)
}

/// A peptide's spectra and the best of them.
pub fn load_peptide(
    path: &Path,
    peptide_id: i64,
) -> Result<(SpectrumPage, Option<SpectrumDetail>), String> {
    drill_peptide(&open(path)?, peptide_id)
}

pub fn load_spectrum(path: &Path, spectrum_id: i64) -> Result<Option<SpectrumDetail>, String> {
    open(path)?
        .spectrum_detail(spectrum_id)
        .map_err(|e| e.to_string())
}

fn drill_group(project: &Project, group_id: i64) -> Result<Drill, String> {
    let group = project.group_detail(group_id).map_err(|e| e.to_string())?;
    let peptide_id = group
        .as_ref()
        .and_then(|g| g.peptides.first())
        .map(|p| p.peptide_id);
    let (spectra, spectrum) = match peptide_id {
        Some(peptide_id) => drill_peptide(project, peptide_id)?,
        None => (SpectrumPage::default(), None),
    };
    Ok(Drill {
        group,
        peptide_id,
        spectra,
        spectrum,
    })
}

fn drill_peptide(
    project: &Project,
    peptide_id: i64,
) -> Result<(SpectrumPage, Option<SpectrumDetail>), String> {
    let spectra = project
        .peptide_spectra(peptide_id, SPECTRA_SHOWN)
        .map_err(|e| e.to_string())?;
    let spectrum = match spectra.rows.first() {
        Some(row) => project
            .spectrum_detail(row.spectrum_id)
            .map_err(|e| e.to_string())?,
        None => None,
    };
    Ok((spectra, spectrum))
}

/// One stick of the plot.
#[derive(Debug, Clone, PartialEq)]
pub struct PlotPeak {
    pub mz: f64,
    /// Intensity relative to the tallest peak, 0..=1.
    pub relative: f32,
    /// The ion it explains (the lowest-charge, b-before-y match).
    pub ion: Option<(IonSeries, String)>,
}

/// One residue of the fragment ladder above the plot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LadderStep {
    pub residue: char,
    pub modified: bool,
    /// Modified by PTM localization (`assigned_on_best` in `site_scores`).
    pub localized: bool,
    /// A b ion ends after this residue.
    pub b: bool,
    /// A y ion starts after this residue.
    pub y: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectrumPlot {
    pub peaks: Vec<PlotPeak>,
    pub min_mz: f64,
    pub max_mz: f64,
    pub ladder: Vec<LadderStep>,
    /// Distinct theoretical ions matched, and how many there were.
    pub matched: usize,
    pub theoretical: usize,
    /// Share of the total intensity the matched peaks carry.
    pub explained: f32,
}

/// Why a spectrum cannot be plotted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoPlot {
    /// No spectra file was attached for its run.
    NoPeaks,
    /// Its stored peaks could not be read.
    Unreadable,
}

/// A mass delta as localization groups it: its Unimod table value when within
/// the snap tolerance, else rounded to the key precision.
fn snapped(delta: f64) -> f64 {
    unimod::snap(delta).unwrap_or_else(|| {
        let factor = 10f64.powi(DELTA_DECIMALS as i32);
        (delta * factor).round() / factor
    })
}

/// The assigned peptide with variable modifications moved to the best
/// localized isoform's sites when PTM localization ran on this spectrum.
fn localized_form(assigned: &tuffold_store::AssignedPeptide) -> tuffold_core::PeptideForm {
    if assigned.site_scores.is_empty() {
        return assigned.form.clone();
    }
    // Match by snapped mass delta, as `tuffold_validate`'s localization groups
    // them: the scores carry the canonical name and snapped mass, while the
    // form's modification may be unnamed with the engine's unsnapped mass.
    let same_modification =
        |a: f64, b: f64| (snapped(a) - snapped(b)).abs() <= unimod::SNAP_TOLERANCE;
    let mut modifications: Vec<tuffold_core::Modification> = assigned
        .form
        .modifications
        .iter()
        .filter(|m| {
            !matches!(m.site, ModSite::Residue(_))
                || !assigned
                    .site_scores
                    .iter()
                    .any(|s| same_modification(s.mass_delta, m.mass_delta))
        })
        .cloned()
        .collect();
    for score in &assigned.site_scores {
        if score.assigned {
            let accession = assigned
                .form
                .modifications
                .iter()
                .find(|m| {
                    matches!(m.site, ModSite::Residue(_))
                        && same_modification(score.mass_delta, m.mass_delta)
                })
                .and_then(|m| m.accession.clone());
            modifications.push(tuffold_core::Modification {
                site: ModSite::Residue(score.position),
                name: score.modification_name.clone(),
                accession,
                mass_delta: score.mass_delta,
            });
        }
    }
    modifications.sort_by_key(|m| m.site);
    tuffold_core::PeptideForm {
        sequence: assigned.form.sequence.clone(),
        modifications,
    }
}

/// Annotate a spectrum's peaks with its assigned peptide's b/y ions.
pub fn spectrum_plot(
    detail: &SpectrumDetail,
    tolerance: Tolerance,
) -> Result<SpectrumPlot, NoPlot> {
    let peaks = match &detail.peaks {
        StoredPeaks::Stored(peaks) => peaks,
        StoredPeaks::Missing => return Err(NoPlot::NoPeaks),
        StoredPeaks::Unreadable => return Err(NoPlot::Unreadable),
    };
    let form = detail.assigned.as_ref().map(localized_form);
    let ions = detail
        .assigned
        .as_ref()
        .zip(form.as_ref())
        .and_then(|(assigned, form)| fragment_ions(form, assigned.charge))
        .unwrap_or_default();
    let mut matches = annotate(&ions, &peaks.mz, &peaks.intensity, tolerance);
    // Lowest charge first, then b before y, so each peak shows its simplest ion.
    matches.sort_by_key(|m| (m.ion.charge, m.ion.series, m.ion.ordinal));
    let tallest = peaks
        .intensity
        .iter()
        .copied()
        .fold(0.0_f32, f32::max)
        .max(f32::MIN_POSITIVE);
    let mut plot_peaks: Vec<PlotPeak> = peaks
        .mz
        .iter()
        .zip(&peaks.intensity)
        .map(|(&mz, &intensity)| PlotPeak {
            mz,
            relative: intensity / tallest,
            ion: None,
        })
        .collect();
    for found in &matches {
        let peak = &mut plot_peaks[found.peak];
        if peak.ion.is_none() {
            peak.ion = Some((found.ion.series, found.ion.label()));
        }
    }
    let total: f64 = peaks.intensity.iter().map(|&i| f64::from(i)).sum();
    let explained: f64 = plot_peaks
        .iter()
        .zip(&peaks.intensity)
        .filter(|(peak, _)| peak.ion.is_some())
        .map(|(_, &intensity)| f64::from(intensity))
        .sum();
    let (min_mz, max_mz) = mz_range(detail, &peaks.mz);
    Ok(SpectrumPlot {
        peaks: plot_peaks,
        min_mz,
        max_mz,
        ladder: ladder(detail, form.as_ref(), &matches),
        matched: matches.len(),
        theoretical: ions.len(),
        explained: if total > 0.0 {
            (explained / total) as f32
        } else {
            0.0
        },
    })
}

/// The plotted m/z range: the peaks, padded 5%, from 0 up.
fn mz_range(detail: &SpectrumDetail, mz: &[f64]) -> (f64, f64) {
    let highest = mz.last().copied().or(detail.precursor_mz).unwrap_or(1000.0);
    (0.0, (highest * 1.05).max(100.0))
}

fn ladder(
    detail: &SpectrumDetail,
    form: Option<&tuffold_core::PeptideForm>,
    matches: &[tuffold_core::fragments::Match],
) -> Vec<LadderStep> {
    let (Some(assigned), Some(form)) = (&detail.assigned, form) else {
        return Vec::new();
    };
    let length = form.sequence.chars().count();
    let has = |series: IonSeries, ordinal: usize| {
        matches
            .iter()
            .any(|m| m.ion.series == series && usize::try_from(m.ion.ordinal) == Ok(ordinal))
    };
    form.sequence
        .chars()
        .enumerate()
        .map(|(index, residue)| {
            let position = index + 1;
            let pos_u32 = u32::try_from(position).unwrap_or(0);
            let localized = assigned
                .site_scores
                .iter()
                .any(|s| s.assigned && s.position == pos_u32);
            LadderStep {
                residue,
                modified: form
                    .modifications
                    .iter()
                    .any(|m| m.site == ModSite::Residue(pos_u32)),
                localized,
                b: position < length && has(IonSeries::B, position),
                y: position < length && has(IonSeries::Y, length - position),
            }
        })
        .collect()
}

/// One line of the coverage map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoverageLine {
    /// 1-based position of the first residue.
    pub start: usize,
    pub residues: Vec<(char, bool)>,
}

/// The sequence in lines of `width`, each residue marked covered or not.
pub fn coverage_lines(
    sequence: &str,
    spans: &[SequenceCoverage],
    width: usize,
) -> Vec<CoverageLine> {
    let residues: Vec<char> = sequence.chars().collect();
    let mask = covered_mask(residues.len(), spans);
    residues
        .chunks(width.max(1))
        .zip(mask.chunks(width.max(1)))
        .enumerate()
        .map(|(line, (chars, covered))| CoverageLine {
            start: line * width.max(1) + 1,
            residues: chars.iter().copied().zip(covered.iter().copied()).collect(),
        })
        .collect()
}

/// `20 ppm`, `0.5 Da`.
pub fn tolerance_label(tolerance: Tolerance) -> String {
    match tolerance {
        Tolerance::Ppm(ppm) => format!("{ppm} ppm"),
        Tolerance::Da(da) => format!("{da} Da"),
    }
}

/// Which protein groups the Proteins and Samples views list. Empty fields
/// keep everything; a group must pass every one that is set.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ProteinFilter {
    /// Text the accession or description contains (any case).
    pub text: String,
    /// Passing spectra the group needs in all.
    pub min_spectra: usize,
    /// The organism (UniProt's `OS=`) its leader must come from.
    pub organism: Option<String>,
    /// A modification one of its unique or razor peptides must carry.
    pub modification: Option<String>,
    /// Only groups the user starred.
    pub starred_only: bool,
    /// A GO slim term (`name (aspect)`) the leader's annotations reach.
    pub go: Option<String>,
}

impl ProteinFilter {
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }

    pub fn matches(&self, row: &tuffold_store::ProteinListRow) -> bool {
        let text = self.text.trim().to_lowercase();
        let text_matches = text.is_empty()
            || row.accession.to_lowercase().contains(&text)
            || row
                .description
                .as_ref()
                .is_some_and(|d| d.to_lowercase().contains(&text));
        text_matches
            && row.spectra >= self.min_spectra
            && self
                .organism
                .as_deref()
                .is_none_or(|wanted| organism(row.description.as_deref()) == Some(wanted))
            && self
                .modification
                .as_ref()
                .is_none_or(|wanted| row.modifications.contains(wanted))
            && (!self.starred_only || row.starred)
            && self
                .go
                .as_ref()
                .is_none_or(|wanted| row.go.contains(wanted))
    }

    /// The rows of `list` that pass, in order.
    pub fn apply(
        &self,
        rows: &[tuffold_store::ProteinListRow],
    ) -> Vec<tuffold_store::ProteinListRow> {
        rows.iter()
            .filter(|row| self.matches(row))
            .cloned()
            .collect()
    }
}

/// The organism in a UniProt description: what follows ` OS=` up to the
/// next ` XX=` field. `None` for other descriptions.
pub fn organism(description: Option<&str>) -> Option<&str> {
    let after = description?.split_once("OS=")?.1;
    let end = after
        .match_indices('=')
        .filter_map(|(at, _)| {
            // A field starts with a space and two capitals before its "=".
            let key = after.get(at.checked_sub(3)?..at)?;
            let mut chars = key.chars();
            (chars.next() == Some(' ') && chars.all(|c| c.is_ascii_uppercase())).then_some(at - 3)
        })
        .next()
        .unwrap_or(after.len());
    Some(after[..end].trim()).filter(|name| !name.is_empty())
}

/// Choices with how many groups have each.
pub type Counted = Vec<(String, usize)>;

/// The filter menus' choices: what the listed groups have, each with how
/// many groups have it, most common first.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Facets {
    pub organisms: Counted,
    pub modifications: Counted,
    pub go: Counted,
}

pub fn facets(rows: &[tuffold_store::ProteinListRow]) -> Facets {
    let mut organisms: std::collections::HashMap<&str, usize> = Default::default();
    let mut modifications: std::collections::HashMap<&str, usize> = Default::default();
    let mut go: std::collections::HashMap<&str, usize> = Default::default();
    for row in rows {
        if let Some(name) = organism(row.description.as_deref()) {
            *organisms.entry(name).or_default() += 1;
        }
        for modification in &row.modifications {
            *modifications.entry(modification).or_default() += 1;
        }
        for term in &row.go {
            *go.entry(term).or_default() += 1;
        }
    }
    let ranked = |counts: std::collections::HashMap<&str, usize>| {
        let mut ranked: Vec<(String, usize)> = counts
            .into_iter()
            .map(|(k, n)| (k.to_string(), n))
            .collect();
        ranked.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        ranked
    };
    Facets {
        organisms: ranked(organisms),
        modifications: ranked(modifications),
        go: ranked(go),
    }
}

/// What the Samples view shows for each group in each biosample, as
/// Scaffold's display options do.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SampleMeasure {
    /// Passing spectra of unique and razor peptides.
    #[default]
    Spectra,
    /// Passing spectra of unique peptides only.
    ExclusiveSpectra,
    /// Distinct unique peptide sequences.
    ExclusivePeptides,
    /// Normalized spectral abundance factor.
    Nsaf,
    /// Exponentially modified protein abundance index.
    Empai,
}

impl SampleMeasure {
    pub const ALL: [Self; 5] = [
        Self::Spectra,
        Self::ExclusiveSpectra,
        Self::ExclusivePeptides,
        Self::Nsaf,
        Self::Empai,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::Spectra => "Spectra",
            Self::ExclusiveSpectra => "Exclusive spectra",
            Self::ExclusivePeptides => "Exclusive peptides",
            Self::Nsaf => "NSAF",
            Self::Empai => "emPAI",
        }
    }

    /// One line under the heading saying what the numbers are.
    pub fn explain(self) -> &'static str {
        match self {
            Self::Spectra => "Passing spectra of each group's unique and razor peptides",
            Self::ExclusiveSpectra => "Passing spectra of peptides only this group has",
            Self::ExclusivePeptides => "Distinct peptide sequences only this group has",
            Self::Nsaf => "Spectra per residue, as a share of the sample's total",
            Self::Empai => "10^(peptides seen / tryptic peptides possible) − 1",
        }
    }

    /// Whole numbers (counts) rather than fractions.
    pub fn is_count(self) -> bool {
        matches!(
            self,
            Self::Spectra | Self::ExclusiveSpectra | Self::ExclusivePeptides
        )
    }
}

/// Each biosample's NSAF denominator: Σ spectra / length over the groups
/// whose leader's length is known.
pub fn saf_totals(list: &ProteinList) -> Vec<f64> {
    (0..list.biosamples.len())
        .map(|sample| {
            list.rows
                .iter()
                .filter_map(|row| {
                    let length = row.length.filter(|&length| length > 0)?;
                    Some(*row.spectra_by_sample.get(sample)? as f64 / length as f64)
                })
                .sum()
        })
        .collect()
}

/// One group's value in one biosample; `None` when it can't be known (NSAF
/// or emPAI without the leader's sequence).
pub fn sample_value(
    row: &tuffold_store::ProteinListRow,
    sample: usize,
    measure: SampleMeasure,
    saf_totals: &[f64],
) -> Option<f64> {
    let count = |values: &[usize]| values.get(sample).map(|&n| n as f64);
    match measure {
        SampleMeasure::Spectra => count(&row.spectra_by_sample),
        SampleMeasure::ExclusiveSpectra => count(&row.exclusive_spectra_by_sample),
        SampleMeasure::ExclusivePeptides => count(&row.exclusive_peptides_by_sample),
        SampleMeasure::Nsaf => {
            let length = row.length.filter(|&length| length > 0)?;
            let total = *saf_totals.get(sample)?;
            let saf = *row.spectra_by_sample.get(sample)? as f64 / length as f64;
            Some(if total > 0.0 { saf / total } else { 0.0 })
        }
        SampleMeasure::Empai => tuffold_quant::measures::empai(
            *row.peptides_by_sample.get(sample)?,
            row.observable_peptides?,
        ),
    }
}

/// The largest value of `measure` anywhere in the list (the heat scale's
/// top), at least the smallest positive value so the scale never divides by
/// zero.
pub fn highest_value<'a>(
    rows: impl IntoIterator<Item = &'a tuffold_store::ProteinListRow>,
    samples: usize,
    measure: SampleMeasure,
    saf_totals: &[f64],
) -> f64 {
    rows.into_iter()
        .flat_map(|row| {
            (0..samples).filter_map(move |sample| sample_value(row, sample, measure, saf_totals))
        })
        .fold(f64::MIN_POSITIVE, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tuffold_core::PeptideForm;
    use tuffold_store::{AssignedPeptide, Peaks};

    fn detail(peaks: StoredPeaks) -> SpectrumDetail {
        SpectrumDetail {
            spectrum_id: 1,
            run: "r".into(),
            key: "scan=1".into(),
            precursor_mz: Some(400.0),
            precursor_charge: Some(2),
            rt_seconds: None,
            peaks,
            assigned: Some(AssignedPeptide {
                peptide_id: 1,
                form: PeptideForm::unmodified("PEPTIDE"),
                localized_peptide: None,
                site_probabilities: None,
                best_ascore: None,
                min_site_probability: None,
                site_scores: Vec::new(),
                is_decoy: false,
                q_value: 0.001,
                probability: 0.999,
                score: 5.0,
                charge: Some(2),
            }),
            psms: Vec::new(),
        }
    }

    #[test]
    fn localized_sites_move_fragment_annotation_and_mark_the_ladder() {
        // Engine reported Oxidation at M1 of MAMK (M[Oxidation]AMK), but PTM
        // localization placed it at M3 (MAM[Oxidation]K).
        // For MAM[Oxidation]K:
        //   b2 (MA) = 131.040485 + 71.037114 + 1.007276 = 203.084875 (unoxidized)
        //   y2 (M[Oxidation]K) = 131.040485 + 15.994915 + 128.094963 + 18.010565 + 1.007276
        //                      = 294.148204
        let mut spec = detail(StoredPeaks::Stored(Peaks {
            mz: vec![203.084_875, 294.148_204],
            intensity: vec![100.0, 80.0],
        }));
        let assigned = spec.assigned.as_mut().unwrap();
        assigned.form = PeptideForm {
            sequence: "MAMK".into(),
            modifications: vec![tuffold_core::Modification {
                site: ModSite::Residue(1),
                name: Some("Oxidation".into()),
                accession: None,
                mass_delta: 15.994_915,
            }],
        };
        assigned.localized_peptide = Some("MAM[Oxidation]K".into());
        assigned.site_scores = vec![
            tuffold_validate::SiteScore {
                position: 1,
                residue: 'M',
                mass_delta: 15.994_915,
                modification_name: Some("Oxidation".into()),
                probability: 0.01,
                ascore: 0.0,
                determining_ions_matched: 0,
                determining_ions_total: 2,
                assigned: false,
                is_decoy: false,
            },
            tuffold_validate::SiteScore {
                position: 3,
                residue: 'M',
                mass_delta: 15.994_915,
                modification_name: Some("Oxidation".into()),
                probability: 0.99,
                ascore: 25.0,
                determining_ions_matched: 2,
                determining_ions_total: 2,
                assigned: true,
                is_decoy: false,
            },
        ];
        let plot = spectrum_plot(&spec, Tolerance::Ppm(20.0)).unwrap();
        assert_eq!(plot.matched, 2);
        let localized: Vec<(char, bool, bool)> = plot
            .ladder
            .iter()
            .map(|s| (s.residue, s.modified, s.localized))
            .collect();
        assert_eq!(
            localized,
            [
                ('M', false, false),
                ('A', false, false),
                ('M', true, true),
                ('K', false, false),
            ]
        );
    }

    #[test]
    fn localization_replaces_an_unnamed_unsnapped_modification_instead_of_doubling_it() {
        // The engine's Phospho on S2 is unnamed with mass 79.96633; the scores
        // carry the snapped mass 79.966331 and the canonical name. Moving it
        // to T3 must leave one modification, keeping the form's accession.
        let mut spec = detail(StoredPeaks::Missing);
        let assigned = spec.assigned.as_mut().unwrap();
        assigned.form = PeptideForm {
            sequence: "PSTK".into(),
            modifications: vec![tuffold_core::Modification {
                site: ModSite::Residue(2),
                name: None,
                accession: Some("UNIMOD:21".into()),
                mass_delta: 79.966_33,
            }],
        };
        let site = |position, residue, probability, assigned| tuffold_validate::SiteScore {
            position,
            residue,
            mass_delta: 79.966_331,
            modification_name: Some("Phospho".into()),
            probability,
            ascore: 0.0,
            determining_ions_matched: 0,
            determining_ions_total: 2,
            assigned,
            is_decoy: false,
        };
        assigned.site_scores = vec![site(2, 'S', 0.01, false), site(3, 'T', 0.99, true)];
        let form = localized_form(assigned);
        assert_eq!(form.modifications.len(), 1, "{:?}", form.modifications);
        let modification = &form.modifications[0];
        assert_eq!(modification.site, ModSite::Residue(3));
        assert_eq!(modification.name.as_deref(), Some("Phospho"));
        assert_eq!(modification.accession.as_deref(), Some("UNIMOD:21"));
    }

    #[test]
    fn the_plot_labels_matched_peaks_and_builds_the_ladder() {
        // PEPTIDE b2 = 227.102633 and y1 = 148.060434 (see fragments tests);
        // 500.0 is noise.
        let peaks = Peaks {
            mz: vec![148.060_434, 227.102_633, 500.0],
            intensity: vec![50.0, 100.0, 25.0],
        };
        let plot =
            spectrum_plot(&detail(StoredPeaks::Stored(peaks)), Tolerance::Ppm(20.0)).unwrap();
        let labels: Vec<Option<&str>> = plot
            .peaks
            .iter()
            .map(|p| p.ion.as_ref().map(|(_, label)| label.as_str()))
            .collect();
        assert_eq!(labels, [Some("y1"), Some("b2"), None]);
        assert_eq!(plot.peaks[1].relative, 1.0);
        assert_eq!((plot.matched, plot.theoretical), (2, 12));
        // (50 + 100) / 175 of the intensity is explained.
        assert!((plot.explained - 150.0 / 175.0).abs() < 1e-6);
        // b2 ends after residue 2 (E); y1 starts after residue 6 (D).
        let marks: Vec<(char, bool, bool)> =
            plot.ladder.iter().map(|s| (s.residue, s.b, s.y)).collect();
        assert_eq!(
            marks,
            [
                ('P', false, false),
                ('E', true, false),
                ('P', false, false),
                ('T', false, false),
                ('I', false, false),
                ('D', false, true),
                ('E', false, false),
            ]
        );
        assert_eq!((plot.min_mz, plot.max_mz), (0.0, 525.0));
    }

    #[test]
    fn spectra_without_peaks_say_why() {
        assert_eq!(
            spectrum_plot(&detail(StoredPeaks::Missing), TOLERANCES[0]),
            Err(NoPlot::NoPeaks)
        );
        assert_eq!(
            spectrum_plot(&detail(StoredPeaks::Unreadable), TOLERANCES[0]),
            Err(NoPlot::Unreadable)
        );
    }

    #[test]
    fn coverage_lines_wrap_and_mark_covered_residues() {
        let spans = [SequenceCoverage {
            start: 2,
            end: 4,
            peptide_id: "p".into(),
        }];
        let lines = coverage_lines("MPEPTIDE", &spans, 5);
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1].start, 6);
        let covered: Vec<bool> = lines
            .iter()
            .flat_map(|l| l.residues.iter().map(|r| r.1))
            .collect();
        assert_eq!(
            covered,
            [false, true, true, true, false, false, false, false]
        );
    }

    #[test]
    fn tolerances_read_naturally() {
        assert_eq!(tolerance_label(TOLERANCES[0]), "20 ppm");
        assert_eq!(tolerance_label(TOLERANCES[1]), "0.5 Da");
    }
}

#[cfg(test)]
mod filter_tests {
    use tuffold_store::ProteinListRow;

    use super::*;

    fn row(
        accession: &str,
        description: &str,
        spectra: usize,
        mods: &[&str],
        starred: bool,
    ) -> ProteinListRow {
        ProteinListRow {
            group_id: 1,
            accession: accession.into(),
            description: Some(description.into()),
            equivalents: 0,
            peptides: 1,
            spectra,
            spectra_by_sample: vec![spectra],
            exclusive_spectra_by_sample: vec![spectra],
            exclusive_peptides_by_sample: vec![1],
            peptides_by_sample: vec![1],
            q_value: None,
            probability: None,
            score: 0.0,
            coverage: None,
            length: None,
            observable_peptides: None,
            modifications: mods.iter().map(|m| m.to_string()).collect(),
            starred,
            go: Vec::new(),
        }
    }

    #[test]
    fn organisms_come_from_uniprot_descriptions() {
        assert_eq!(
            organism(Some(
                "Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2"
            )),
            Some("Homo sapiens")
        );
        assert_eq!(
            organism(Some(
                "Enolase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"
            )),
            Some("Saccharomyces cerevisiae (strain ATCC 204508 / S288c)")
        );
        assert_eq!(organism(Some("Demo protein 2 (synthetic)")), None);
        assert_eq!(organism(None), None);
    }

    #[test]
    fn a_row_must_pass_every_filter_set() {
        let human = row(
            "sp|P02768|ALBU_HUMAN",
            "Serum albumin OS=Homo sapiens OX=9606",
            12,
            &["Oxidation"],
            true,
        );
        let yeast = row(
            "sp|P00924|ENO1_YEAST",
            "Enolase 1 OS=Saccharomyces cerevisiae OX=559292",
            3,
            &[],
            false,
        );
        let rows = vec![human.clone(), yeast.clone()];
        let pick = |filter: ProteinFilter| -> Vec<String> {
            filter
                .apply(&rows)
                .into_iter()
                .map(|r| r.accession)
                .collect()
        };
        assert_eq!(pick(ProteinFilter::default()).len(), 2);
        assert_eq!(
            pick(ProteinFilter {
                text: "ENOLASE".into(),
                ..Default::default()
            }),
            [yeast.accession.as_str()]
        );
        assert_eq!(
            pick(ProteinFilter {
                text: "p02768".into(),
                ..Default::default()
            }),
            [human.accession.as_str()]
        );
        assert_eq!(
            pick(ProteinFilter {
                min_spectra: 4,
                ..Default::default()
            }),
            [human.accession.as_str()]
        );
        assert_eq!(
            pick(ProteinFilter {
                organism: Some("Saccharomyces cerevisiae".into()),
                ..Default::default()
            }),
            [yeast.accession.as_str()]
        );
        assert_eq!(
            pick(ProteinFilter {
                modification: Some("Oxidation".into()),
                ..Default::default()
            }),
            [human.accession.as_str()]
        );
        assert_eq!(
            pick(ProteinFilter {
                starred_only: true,
                ..Default::default()
            }),
            [human.accession.as_str()]
        );
        // Both set: nothing passes both.
        assert!(
            pick(ProteinFilter {
                starred_only: true,
                min_spectra: 20,
                ..Default::default()
            })
            .is_empty()
        );
        let found = facets(&rows);
        assert_eq!(found.organisms.len(), 2);
        assert_eq!(found.modifications, [("Oxidation".to_string(), 1)]);
        assert!(found.go.is_empty());
    }
}

#[cfg(test)]
mod measure_tests {
    use tuffold_store::{Biosample, ProteinListRow};

    use super::*;

    fn row(group_id: i64, spectra: [usize; 2], length: Option<usize>) -> ProteinListRow {
        ProteinListRow {
            group_id,
            accession: format!("P{group_id}"),
            description: None,
            equivalents: 0,
            peptides: 2,
            spectra: spectra.iter().sum(),
            spectra_by_sample: spectra.to_vec(),
            exclusive_spectra_by_sample: vec![spectra[0], 0],
            exclusive_peptides_by_sample: vec![1, 0],
            peptides_by_sample: vec![2, 1],
            q_value: Some(0.001),
            probability: Some(0.999),
            score: 1.0,
            coverage: None,
            length,
            observable_peptides: length.map(|_| 4),
            modifications: Vec::new(),
            starred: false,
            go: Vec::new(),
        }
    }

    fn list() -> ProteinList {
        let sample = |id, name: &str| Biosample {
            id,
            name: name.into(),
            category: None,
        };
        ProteinList {
            biosamples: vec![sample(1, "a"), sample(2, "b")],
            // P1: 4 and 2 spectra, length 100; P2: 6 and 0, length 300;
            // P3: 5 and 5, no sequence.
            rows: vec![
                row(1, [4, 2], Some(100)),
                row(2, [6, 0], Some(300)),
                row(3, [5, 5], None),
            ],
            ..ProteinList::default()
        }
    }

    #[test]
    fn nsaf_divides_by_length_then_by_the_samples_total() {
        let list = list();
        // Sample a: 4/100 + 6/300 = 0.06; sample b: 2/100 + 0 = 0.02.
        let totals = saf_totals(&list);
        assert!((totals[0] - 0.06).abs() < 1e-12 && (totals[1] - 0.02).abs() < 1e-12);
        let nsaf = |row: usize, sample| {
            sample_value(&list.rows[row], sample, SampleMeasure::Nsaf, &totals)
        };
        // P1 in a: 0.04 / 0.06 = 2/3; P2 in a: 0.02 / 0.06 = 1/3; P1 in b: 1.
        assert!((nsaf(0, 0).unwrap() - 2.0 / 3.0).abs() < 1e-12);
        assert!((nsaf(1, 0).unwrap() - 1.0 / 3.0).abs() < 1e-12);
        assert!((nsaf(0, 1).unwrap() - 1.0).abs() < 1e-12);
        // No sequence, no NSAF.
        assert_eq!(nsaf(2, 0), None);
    }

    #[test]
    fn counts_and_empai_come_from_the_row() {
        let list = list();
        let totals = saf_totals(&list);
        let value = |measure| sample_value(&list.rows[0], 0, measure, &totals);
        assert_eq!(value(SampleMeasure::Spectra), Some(4.0));
        assert_eq!(value(SampleMeasure::ExclusiveSpectra), Some(4.0));
        assert_eq!(value(SampleMeasure::ExclusivePeptides), Some(1.0));
        // 2 of 4 observable peptides: 10^0.5 − 1.
        assert!((value(SampleMeasure::Empai).unwrap() - (10_f64.sqrt() - 1.0)).abs() < 1e-12);
        assert_eq!(
            sample_value(&list.rows[2], 0, SampleMeasure::Empai, &totals),
            None
        );
        assert_eq!(
            highest_value(&list.rows, 2, SampleMeasure::Spectra, &totals),
            6.0
        );
        assert_eq!(
            highest_value(&list.rows, 2, SampleMeasure::Nsaf, &totals),
            1.0
        );
    }
}
