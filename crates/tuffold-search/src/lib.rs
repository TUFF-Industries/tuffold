//! Built-in database search: Sage's search core (Michael Lazear, MIT) —
//! its fragment index and hyperscore — over Tuffold's own records.
//!
//! A [`Searcher`] digests a FASTA once (targets only; Sage generates
//! reversed decoys tagged [`DECOY_TAG`], so every decoy follows one
//! convention), then searches one run's spectra at a time so only one run
//! is held in memory. Each run becomes an ordinary [`ImportBatch`]
//! (software "Sage"): rank-1 PSMs with Sage's sub-scores, decoys included.
//! Sage's own FDR and rescoring are not used — Tuffold validates the PSMs
//! like any other engine's.

use std::collections::{BTreeSet, HashMap};

use rayon::prelude::*;
use sage_core::database::IndexedDatabase;
use sage_core::fasta::Fasta;
use sage_core::mass::Tolerance;
use sage_core::scoring::{Feature, Scorer};
use sage_core::spectrum::{Precursor, RawSpectrum, Representation, SpectrumProcessor};
use tuffold_core::centroid::centroid;
use tuffold_core::chemistry::{mz, peptide_mass};
use tuffold_core::{
    BioSample, Evidence, ImportBatch, MsRun, PeakSpectrum, Protein, Psm, Score, Software, Spectrum,
};
use tuffold_io::fasta::FastaRecord;

mod budget;
mod mods;
mod settings;

pub use budget::{IndexEstimate, default_memory_budget};
use settings::FRAGMENT_MZ;
pub use settings::{
    Enzyme, FragmentTolerance, Label, ModRule, SearchOptions, SearchSettings, Site,
};

/// The Sage release whose search core is built in (tag v0.14.7; its
/// `sage-core` crate manifest still reads 0.14.6).
pub const SAGE_VERSION: &str = "0.14.7";
/// Prefix of the decoy proteins Sage generates (reversed sequences).
pub const DECOY_TAG: &str = "rev_";

/// The most intense peaks kept per spectrum.
const TOP_PEAKS: usize = 150;
/// The fewest peaks a spectrum needs to be searched.
pub const MIN_PEAKS: usize = 15;
/// Matched fragments a PSM needs to be reported at all.
const MIN_MATCHED: u16 = 4;
/// Precursor isotope errors considered (monoisotopic peak picked 1 low to 3
/// high).
const ISOTOPE_ERRORS: (i8, i8) = (-1, 3);
const MAX_FRAGMENT_CHARGE: u8 = 2;

#[derive(Debug, thiserror::Error)]
pub enum SearchError {
    #[error("the FASTA has no target proteins with sequences to search")]
    NoTargets,
    #[error("invalid search settings: {0}")]
    Settings(String),
    #[error("{}", too_large(.estimate, *.budget, .settings))]
    TooLarge {
        estimate: IndexEstimate,
        budget: u64,
        settings: Box<SearchSettings>,
    },
}

/// "This search's index needs at least 24.1 GB ...": the numbers, and what
/// in the settings makes it big.
fn too_large(estimate: &IndexEstimate, budget: u64, settings: &SearchSettings) -> String {
    let gb = |bytes: u64| bytes as f64 / 1e9;
    let needs = if estimate.complete {
        format!(
            "about {:.1} GB ({:.1} million peptide forms with decoys)",
            gb(estimate.peak_bytes),
            estimate.forms as f64 / 1e6
        )
    } else {
        format!(
            "over {:.1} GB (counting stopped there, at {:.1} million peptide forms)",
            gb(estimate.peak_bytes),
            estimate.forms as f64 / 1e6
        )
    };
    let mut advice = Vec::new();
    if settings.semi_specific {
        advice.push("search fully specific instead of semi-specific".to_string());
    }
    if settings.missed_cleavages > 2 {
        advice.push(format!(
            "allow fewer than {} missed cleavages",
            settings.missed_cleavages
        ));
    }
    if settings.variable.len() > 1 || settings.max_variable_mods > 2 {
        advice.push("search fewer variable modifications".to_string());
    }
    advice.push("search a smaller FASTA".to_string());
    advice.push("raise the memory budget only if that much memory is really free".to_string());
    format!(
        "this search's index needs {needs}, but its memory budget is {:.1} GB, so it was not \
         started (it would make the computer swap or stop it); {}",
        gb(budget),
        advice.join(", or ")
    )
}

/// One run's search: the batch to import and what was (not) searched.
#[derive(Debug, Clone)]
pub struct RunSearch {
    pub batch: ImportBatch,
    /// MS/MS spectra scored.
    pub searched: usize,
    /// Of those, profile spectra centroided first.
    pub centroided: usize,
    pub skipped: Skipped,
}

/// MS/MS spectra not searched, by reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Skipped {
    /// No usable precursor m/z.
    pub no_precursor: usize,
    /// Fewer than [`MIN_PEAKS`] peaks (after centroiding, for profile data).
    pub too_few_peaks: usize,
    /// A spectrum key seen earlier in the same run (the first one is used).
    pub repeated_key: usize,
}

impl Skipped {
    pub fn total(&self) -> usize {
        self.no_precursor + self.too_few_peaks + self.repeated_key
    }

    /// "3 without a precursor m/z, 2 with fewer than 15 peaks", or `None` when nothing
    /// was skipped.
    pub fn describe(&self) -> Option<String> {
        let parts: Vec<String> = [
            (self.no_precursor, "without a precursor m/z".to_string()),
            (
                self.too_few_peaks,
                format!("with fewer than {MIN_PEAKS} peaks"),
            ),
            (
                self.repeated_key,
                "repeating an earlier spectrum's scan".to_string(),
            ),
        ]
        .into_iter()
        .filter(|(count, _)| *count > 0)
        .map(|(count, reason)| format!("{count} {reason}"))
        .collect();
        (!parts.is_empty()).then(|| parts.join(", "))
    }
}

/// A digested database, ready to search runs against.
pub struct Searcher {
    database: IndexedDatabase,
    settings: SearchSettings,
    /// Target accession -> (description, sequence).
    targets: HashMap<String, (Option<String>, String)>,
    dropped_decoys: usize,
    excluded_targets: Vec<String>,
    index: IndexEstimate,
    budget: u64,
}

impl Searcher {
    /// Digest the FASTA's targets. Decoys already in it are dropped (and
    /// counted): Sage generates its own, tagged [`DECOY_TAG`].
    pub fn new(records: &[FastaRecord], settings: SearchSettings) -> Result<Self, SearchError> {
        settings.check()?;
        let dropped_decoys = records.iter().filter(|r| r.is_decoy).count();
        let mut targets: HashMap<String, (Option<String>, String)> = HashMap::new();
        let mut excluded_targets = Vec::new();
        for record in records
            .iter()
            .filter(|r| !r.is_decoy && !r.sequence.is_empty())
        {
            let Some(accession) = record.identifiers.first() else {
                continue;
            };
            // Sage treats any accession containing the decoy tag as a decoy,
            // and one accession can only be one protein: report, never drop
            // silently.
            if accession.contains(DECOY_TAG) || targets.contains_key(accession) {
                excluded_targets.push(accession.clone());
                continue;
            }
            targets.insert(
                accession.clone(),
                (record.description.clone(), record.sequence.clone()),
            );
        }
        if targets.is_empty() {
            return Err(SearchError::NoTargets);
        }
        // Sorted, so the digest (and the search) is deterministic.
        let fasta_text: String = targets
            .iter()
            .collect::<std::collections::BTreeMap<_, _>>()
            .into_iter()
            .map(|(accession, (_, sequence))| format!(">{accession}\n{sequence}\n"))
            .collect();
        let parameters = settings.builder()?.make_parameters();
        let fasta = Fasta::parse(fasta_text, DECOY_TAG, true);
        let budget = settings.memory_budget.unwrap_or_else(default_memory_budget);
        let index = budget::estimate(
            &fasta.targets,
            &parameters,
            budget.saturating_mul(budget::COUNT_PAST_BUDGET),
        );
        if index.peak_bytes > budget {
            return Err(SearchError::TooLarge {
                estimate: index,
                budget,
                settings: Box::new(settings),
            });
        }
        let database = parameters.build(fasta);
        Ok(Self {
            database,
            settings,
            targets,
            dropped_decoys,
            excluded_targets,
            index,
            budget,
        })
    }

    /// Target entries left out: their accession contains [`DECOY_TAG`] or
    /// repeats an earlier entry's.
    pub fn excluded_targets(&self) -> &[String] {
        &self.excluded_targets
    }

    /// What the index was estimated to need before it was built, and the
    /// budget it was checked against.
    pub fn index_estimate(&self) -> (IndexEstimate, u64) {
        (self.index, self.budget)
    }

    /// Target and decoy peptide forms in the index.
    pub fn peptides(&self) -> usize {
        self.database.peptides.len()
    }

    pub fn targets(&self) -> usize {
        self.targets.len()
    }

    /// Decoys found in the FASTA and replaced by Sage's.
    pub fn dropped_decoys(&self) -> usize {
        self.dropped_decoys
    }

    /// Search one run's MS/MS spectra. `run` names the run, its biosample,
    /// and the spectra file, exactly as the batch will record them.
    pub fn search_run(&self, run: &MsRun, spectra: &[PeakSpectrum]) -> RunSearch {
        let scorer = Scorer {
            db: &self.database,
            precursor_tol: Tolerance::Ppm(
                -self.settings.precursor_ppm,
                self.settings.precursor_ppm,
            ),
            fragment_tol: self.settings.fragment.sage(),
            min_matched_peaks: MIN_MATCHED,
            min_isotope_err: ISOTOPE_ERRORS.0,
            max_isotope_err: ISOTOPE_ERRORS.1,
            min_precursor_charge: self.settings.charges.0,
            max_precursor_charge: self.settings.charges.1,
            max_fragment_charge: Some(MAX_FRAGMENT_CHARGE),
            min_fragment_mass: FRAGMENT_MZ.0,
            max_fragment_mass: FRAGMENT_MZ.1,
            chimera: false,
            report_psms: 1,
            wide_window: false,
            annotate_matches: false,
        };
        let processor = SpectrumProcessor::new(TOP_PEAKS, FRAGMENT_MZ.0, FRAGMENT_MZ.1, false);
        let (searchable, skipped) = searchable(spectra);
        let centroided = searchable.iter().filter(|s| s.centroids.is_some()).count();
        let mut hits: Vec<(usize, Psm)> = searchable
            .par_iter()
            .flat_map_iter(|prepared| {
                let ix = prepared.ix;
                let spectrum = &spectra[ix];
                let processed = processor.process(raw_spectrum(spectrum, prepared.peaks(spectrum)));
                scorer
                    .score(&processed)
                    .into_iter()
                    .filter(|feature| feature.rank == 1)
                    .map(move |feature| (ix, self.psm(&run.name, spectrum, &feature)))
            })
            .collect();
        // Parallel collection order varies; the batch must not.
        hits.sort_by(|a, b| a.1.spectrum.cmp(&b.1.spectrum));
        RunSearch {
            batch: self.batch(run, spectra, hits),
            searched: searchable.len(),
            centroided,
            skipped,
        }
    }

    fn psm(&self, run: &str, spectrum: &PeakSpectrum, feature: &Feature) -> Psm {
        let peptide = &self.database[feature.peptide_idx];
        let form = mods::peptide_form(peptide, &self.settings);
        let charge = i16::from(feature.charge);
        let calculated_mz = peptide_mass(&form.sequence).map(|mass| {
            let deltas: f64 = form.modifications.iter().map(|m| m.mass_delta).sum();
            // peptide_mass includes water; mz() adds the protons.
            mz(mass + deltas, feature.charge)
        });
        let scores = [
            ("hyperscore", feature.hyperscore),
            // With a single candidate Sage scores the runner-up as 0, so
            // delta_next is the hyperscore itself (its definition).
            ("delta_next", feature.delta_next),
            ("matched_peaks", f64::from(feature.matched_peaks)),
            ("longest_b", f64::from(feature.longest_b)),
            ("longest_y", f64::from(feature.longest_y)),
            (
                "matched_intensity_pct",
                f64::from(feature.matched_intensity_pct),
            ),
        ]
        .into_iter()
        // A score Sage could not compute is left out, never invented.
        .filter(|(_, value)| value.is_finite())
        .map(|(name, value)| Score {
            name: name.to_string(),
            accession: None,
            value,
        })
        .collect();
        Psm {
            ms_run: run.to_string(),
            spectrum: spectrum.key.clone(),
            rank: 1,
            peptide: form,
            charge,
            calculated_mz,
            experimental_mz: spectrum.precursor_mz,
            pass_threshold: None,
            scores,
            evidence: peptide
                .proteins
                .iter()
                .map(|accession| Evidence {
                    // Sage keeps the target's accession on a decoy peptide
                    // and labels it on output, as here.
                    accession: if peptide.decoy {
                        format!("{DECOY_TAG}{accession}")
                    } else {
                        accession.to_string()
                    },
                    start: None,
                    end: None,
                    pre: None,
                    post: None,
                })
                .collect(),
        }
    }

    fn batch(&self, run: &MsRun, spectra: &[PeakSpectrum], hits: Vec<(usize, Psm)>) -> ImportBatch {
        let accessions: BTreeSet<&str> = hits
            .iter()
            .flat_map(|(_, psm)| psm.evidence.iter().map(|e| e.accession.as_str()))
            .collect();
        let proteins = accessions
            .into_iter()
            .map(|accession| {
                let target = self.targets.get(accession);
                Protein {
                    accession: accession.to_string(),
                    description: target.and_then(|(description, _)| description.clone()),
                    sequence: target.map(|(_, sequence)| sequence.clone()),
                    is_decoy: target.is_none(),
                    decoy_explicit: false,
                }
            })
            .collect();
        ImportBatch {
            software: Software {
                name: "Sage".into(),
                version: Some(SAGE_VERSION.into()),
            },
            biosamples: vec![BioSample {
                name: run.biosample.clone(),
                category: None,
            }],
            ms_runs: vec![run.clone()],
            spectra: hits
                .iter()
                .map(|&(ix, _)| &spectra[ix])
                .map(|spectrum| Spectrum {
                    ms_run: run.name.clone(),
                    key: spectrum.key.clone(),
                    source_id: Some(spectrum.native_id.clone()),
                    precursor_mz: spectrum.precursor_mz,
                    precursor_charge: spectrum.precursor_charge,
                    rt_seconds: spectrum.rt_seconds,
                })
                .collect(),
            proteins,
            psms: hits.into_iter().map(|(_, psm)| psm).collect(),
        }
    }
}

/// Indices of the MS/MS spectra worth searching (first of each key), and
/// what was left out.
/// A spectrum to search, with its centroids when the file had profile data.
struct Prepared {
    ix: usize,
    centroids: Option<(Vec<f64>, Vec<f32>)>,
}

impl Prepared {
    fn peaks<'a>(&'a self, spectrum: &'a PeakSpectrum) -> (&'a [f64], &'a [f32]) {
        match &self.centroids {
            Some((mz, intensity)) => (mz, intensity),
            None => (&spectrum.mz, &spectrum.intensity),
        }
    }
}

/// The MS/MS spectra worth searching (the first of each key), centroiding
/// profile data, and what was left out.
fn searchable(spectra: &[PeakSpectrum]) -> (Vec<Prepared>, Skipped) {
    let mut seen = std::collections::HashSet::new();
    let mut skipped = Skipped::default();
    let mut keep = Vec::new();
    for (ix, spectrum) in spectra.iter().enumerate().filter(|(_, s)| s.ms_level == 2) {
        if !seen.insert(&spectrum.key) {
            skipped.repeated_key += 1;
            continue;
        }
        if !spectrum
            .precursor_mz
            .is_some_and(|mz| mz.is_finite() && mz > 0.0)
        {
            skipped.no_precursor += 1;
            continue;
        }
        let centroids = spectrum
            .profile
            .then(|| centroid(&spectrum.mz, &spectrum.intensity));
        let peaks = centroids
            .as_ref()
            .map_or(spectrum.mz.len(), |(mz, _)| mz.len());
        if peaks < MIN_PEAKS {
            skipped.too_few_peaks += 1;
        } else {
            keep.push(Prepared { ix, centroids });
        }
    }
    (keep, skipped)
}

fn raw_spectrum(spectrum: &PeakSpectrum, (mz, intensity): (&[f64], &[f32])) -> RawSpectrum {
    RawSpectrum {
        file_id: 0,
        ms_level: spectrum.ms_level,
        id: spectrum.native_id.clone(),
        precursors: vec![Precursor {
            mz: spectrum.precursor_mz.unwrap_or_default() as f32,
            intensity: None,
            charge: spectrum
                .precursor_charge
                .and_then(|z| u8::try_from(z).ok())
                .filter(|&z| z > 0),
            spectrum_ref: None,
            isolation_window: None,
            inverse_ion_mobility: None,
        }],
        representation: Representation::Centroid,
        scan_start_time: spectrum.rt_seconds.unwrap_or_default() as f32 / 60.0,
        ion_injection_time: 0.0,
        total_ion_current: intensity.iter().sum(),
        mz: mz.iter().map(|&mz| mz as f32).collect(),
        intensity: intensity.to_vec(),
    }
}

#[cfg(test)]
mod tests;
