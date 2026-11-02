//! A synthetic demo project, so Tuffold can be explored without lab data.
//!
//! Everything is generated from a fixed seed and labelled synthetic (source
//! format `synthetic`, software "Tuffold demo (synthetic)"):
//! - a random proteome with Swiss-Prot-like residue frequencies, reversed
//!   decoys, and isoforms that share most of a parent's peptides;
//! - four biosamples (two control, two treated, one run each); 150 proteins
//!   are expressed, 30 of them four times higher in treated and 20 four times
//!   lower;
//! - tryptic spectra: a share are correct (peaks from the reported peptide's
//!   real b/y fragment masses, at a per-spectrum quality), the rest report a
//!   random target or decoy peptide while their peaks come from an unseen one;
//! - scores computed from how many fragment ions the reported peptide
//!   explains, then validated like any import;
//! - MS1 scans, one before each MS/MS scan. Each correctly identified
//!   precursor elutes as a Gaussian peak (σ 6 s) centred on its
//!   identifications. Its height is the protein's abundance × the fold
//!   change × a fixed per-peptide ionization factor × a per-run loading
//!   factor (up to ±40%) × 15% noise, spread over averagine-like isotope
//!   peaks at 2 ppm error, among random noise peaks. Wrongly identified
//!   spectra have no signal at their reported m/z. The MS1 data comes from
//!   its own random stream, so the identifications above don't change.
//!
//! Nothing here is ever mixed into user data: it only creates new projects.

use std::collections::HashMap;
use std::path::Path;

use tuffold_core::fragments::{Tolerance, annotate, fragment_ions};
use tuffold_core::{
    BioSample, Evidence, ImportBatch, ModSite, Modification, MsRun, PeakSpectrum, PeptideForm,
    Protein, Psm, Score, Software, SourceFormat, Spectrum, SpectrumKey,
};
use tuffold_store::{ImportSource, Project, StoreError, ValidationSummary};
use tuffold_validate::Settings;

use crate::import::attach_scans;
use crate::validation::{ValidateError, validate_project};

pub const DEMO_SEED: u64 = 20_261_008;
const TARGETS: usize = 400;
const EXPRESSED: usize = 150;
const ISOFORMS: usize = 20;
const UP: std::ops::Range<usize> = 0..30;
const DOWN: std::ops::Range<usize> = 30..50;
const SPECTRA_PER_RUN: u32 = 1200;
const CORRECT_SHARE: f64 = 0.55;
const SAMPLES: [(&str, &str); 4] = [
    ("control-A", "control"),
    ("control-B", "control"),
    ("treated-A", "treated"),
    ("treated-B", "treated"),
];
/// Residues by approximate Swiss-Prot frequency (percent).
const RESIDUES: [(char, f64); 20] = [
    ('A', 8.25),
    ('R', 5.53),
    ('N', 4.06),
    ('D', 5.45),
    ('C', 1.37),
    ('Q', 3.93),
    ('E', 6.75),
    ('G', 7.07),
    ('H', 2.27),
    ('I', 5.96),
    ('L', 9.66),
    ('K', 5.84),
    ('M', 2.42),
    ('F', 3.86),
    ('P', 4.70),
    ('S', 6.56),
    ('T', 5.34),
    ('W', 1.08),
    ('Y', 2.92),
    ('V', 6.87),
];
const OXIDATION: f64 = 15.994_915;
/// Seconds between MS/MS scans; an MS1 scan sits halfway before each.
const SCAN_SECONDS: f64 = 2.5;
const ELUTION_SIGMA: f64 = 6.0;
/// Per-run loading factors, so normalization has something to correct.
const LOADING: [f64; 4] = [1.0, 1.4, 0.7, 1.1];
const NOISE_PEAKS: usize = 40;

#[derive(Debug, thiserror::Error)]
pub enum DemoError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error(transparent)]
    Validate(#[from] ValidateError),
    #[error("attaching the demo spectra failed: {0}")]
    Attach(String),
}

/// Create a new demo project at `path` (refusing to overwrite) and validate it.
pub fn create_demo_project(path: &Path) -> Result<ValidationSummary, DemoError> {
    let project = Project::create(path)?;
    let result = fill_demo(project);
    if result.is_err() {
        // Never leave a half-built demo behind.
        let _ = std::fs::remove_file(path);
    }
    result
}

fn fill_demo(mut project: Project) -> Result<ValidationSummary, DemoError> {
    let demo = generate(DEMO_SEED);
    let source = ImportSource {
        path: format!("Tuffold demo (synthetic data, seed {DEMO_SEED})"),
        format: SourceFormat::Synthetic,
    };
    project.import(&source, &demo.batch)?;
    for (run, scans) in demo.peaks {
        attach_scans(&mut project, &run, scans.into_iter().map(Ok)).map_err(DemoError::Attach)?;
    }
    Ok(validate_project(&mut project, &Settings::default())?)
}

/// The generated search and every run's scans (MS1 and MS/MS, in time
/// order).
pub struct Demo {
    pub batch: ImportBatch,
    pub peaks: Vec<(String, Vec<PeakSpectrum>)>,
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn normal(&mut self) -> f64 {
        let (u1, u2) = (self.uniform(), self.uniform());
        (-2.0 * (1.0 - u1).ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    fn index(&mut self, n: usize) -> usize {
        usize::try_from(self.next() % n as u64).unwrap_or_default()
    }

    fn residue(&mut self) -> char {
        let total: f64 = RESIDUES.iter().map(|(_, f)| f).sum();
        let mut pick = self.uniform() * total;
        for (residue, frequency) in RESIDUES {
            if pick < frequency {
                return residue;
            }
            pick -= frequency;
        }
        'A'
    }
}

/// Tryptic peptides of 7..=25 residues, no missed cleavages.
pub fn digest(sequence: &str) -> Vec<String> {
    let residues: Vec<char> = sequence.chars().collect();
    let mut peptides = Vec::new();
    let mut start = 0;
    for i in 0..residues.len() {
        let cleaves =
            matches!(residues[i], 'K' | 'R') && residues.get(i + 1).is_none_or(|&next| next != 'P');
        if cleaves || i + 1 == residues.len() {
            if (7..=25).contains(&(i + 1 - start)) {
                peptides.push(residues[start..=i].iter().collect());
            }
            start = i + 1;
        }
    }
    peptides
}

struct Proteome {
    proteins: Vec<Protein>,
    /// Per protein (targets then decoys), its digest.
    digests: Vec<Vec<String>>,
    /// Peptide -> accessions containing it.
    evidence: HashMap<String, Vec<String>>,
}

fn proteome(rng: &mut Rng) -> Proteome {
    let mut sequences: Vec<String> = (0..TARGETS)
        .map(|_| {
            let length = 150 + rng.index(550);
            std::iter::once('M')
                .chain((1..length).map(|_| rng.residue()))
                .collect()
        })
        .collect();
    // Isoforms: an unexpressed copy of a parent with a stretch replaced.
    for isoform in 0..ISOFORMS {
        let parent: Vec<char> = sequences[isoform * 3].chars().collect();
        let (from, to) = (parent.len() / 3, parent.len() / 3 + parent.len() / 4);
        let replaced: String = parent
            .iter()
            .enumerate()
            .map(|(i, &residue)| {
                if (from..to).contains(&i) {
                    rng.residue()
                } else {
                    residue
                }
            })
            .collect();
        sequences[EXPRESSED + isoform] = replaced;
    }
    let mut proteins = Vec::new();
    for (index, sequence) in sequences.iter().enumerate() {
        let accession = format!("sp|DEMO{index:04}|DEMO{index}_SYNTH");
        proteins.push(Protein {
            description: Some(format!("Demo protein {index} (synthetic)")),
            sequence: Some(sequence.clone()),
            is_decoy: false,
            decoy_explicit: false,
            accession,
        });
    }
    for (index, sequence) in sequences.iter().enumerate() {
        proteins.push(Protein {
            accession: format!("rev_sp|DEMO{index:04}|DEMO{index}_SYNTH"),
            description: None,
            sequence: Some(sequence.chars().rev().collect()),
            is_decoy: true,
            decoy_explicit: false,
        });
    }
    let digests: Vec<Vec<String>> = proteins
        .iter()
        .map(|p| digest(p.sequence.as_deref().unwrap_or_default()))
        .collect();
    let mut evidence: HashMap<String, Vec<String>> = HashMap::new();
    for (protein, peptides) in proteins.iter().zip(&digests) {
        for peptide in peptides {
            let accessions = evidence.entry(peptide.clone()).or_default();
            if !accessions.contains(&protein.accession) {
                accessions.push(protein.accession.clone());
            }
        }
    }
    Proteome {
        proteins,
        digests,
        evidence,
    }
}

/// Generate the demo search.
pub fn generate(seed: u64) -> Demo {
    generate_with(seed, false)
}

/// [`generate`], optionally with chromatography: each run's spectra ordered
/// by a fixed elution position per peptide (plus a little noise), so a
/// precursor elutes at about the same time in every run, as on a real
/// gradient. Off for the demo project: its numbers and pictures were
/// made without it. Match-between-runs needs it.
pub fn generate_with(seed: u64, gradient: bool) -> Demo {
    let mut rng = Rng(seed);
    let proteome = proteome(&mut rng);
    let abundance: Vec<f64> = (0..EXPRESSED).map(|_| (1.2 * rng.normal()).exp()).collect();
    let mut batch = ImportBatch {
        software: Software {
            name: "Tuffold demo (synthetic)".into(),
            version: Some(env!("CARGO_PKG_VERSION").into()),
        },
        proteins: proteome.proteins.clone(),
        ..ImportBatch::default()
    };
    // Peptide -> the expressed proteins containing it: a shared peptide's
    // MS1 signal is the sum of theirs.
    let mut expressed_in: HashMap<&str, Vec<usize>> = HashMap::new();
    for (protein, digest) in proteome.digests.iter().take(EXPRESSED).enumerate() {
        for peptide in digest {
            expressed_in
                .entry(peptide.as_str())
                .or_default()
                .push(protein);
        }
    }
    let mut peaks = Vec::new();
    for (run_index, (name, category)) in SAMPLES.into_iter().enumerate() {
        batch.biosamples.push(BioSample {
            name: name.into(),
            category: Some(category.into()),
        });
        batch.ms_runs.push(MsRun {
            name: name.into(),
            biosample: name.into(),
            spectra_file: None,
        });
        let weights: Vec<f64> = abundance
            .iter()
            .enumerate()
            .map(|(i, a)| a * fold(i, category))
            .collect();
        let mut msms = Vec::new();
        let mut eluting: Vec<Identified> = Vec::new();
        let mut drawn: Vec<_> = (1..=SPECTRA_PER_RUN)
            .map(|scan| spectrum(&mut rng, &proteome, &weights, name, scan))
            .collect();
        if gradient {
            let mut placed: Vec<(f64, _)> = drawn
                .into_iter()
                .map(|drawn| {
                    (
                        elution_position(&drawn.0.peptide.sequence) + 0.004 * rng.normal(),
                        drawn,
                    )
                })
                .collect();
            placed.sort_by(|a, b| a.0.total_cmp(&b.0));
            drawn = placed
                .into_iter()
                .zip(1..)
                .map(|((_, drawn), scan)| renumbered(drawn, scan))
                .collect();
        }
        for (psm, spectrum, measured, identified) in drawn {
            batch.spectra.push(spectrum);
            batch.psms.push(psm);
            msms.extend(measured);
            eluting.extend(identified);
        }
        let mut ms1_rng = Rng(seed ^ (0x4D53_3100 + run_index as u64));
        let elution = elution(
            &eluting,
            &expressed_in,
            &weights,
            LOADING[run_index],
            &mut ms1_rng,
        );
        peaks.push((
            name.to_string(),
            interleave_ms1(msms, &elution, &mut ms1_rng),
        ));
    }
    Demo { batch, peaks }
}

/// A peptide's place on the gradient, 0..1, fixed by its sequence (FNV-1a).
fn elution_position(sequence: &str) -> f64 {
    let hash = sequence.bytes().fold(0x811c_9dc5_u32, |h, b| {
        (h ^ u32::from(b)).wrapping_mul(0x0100_0193)
    });
    f64::from(hash) / f64::from(u32::MAX)
}

/// A drawn spectrum moved to scan `scan` (key, id, and time follow).
fn renumbered(
    (mut psm, mut spectrum, measured, identified): (
        Psm,
        Spectrum,
        Option<PeakSpectrum>,
        Option<Identified>,
    ),
    scan: u32,
) -> (Psm, Spectrum, Option<PeakSpectrum>, Option<Identified>) {
    let key = SpectrumKey::Scan(scan);
    let rt = f64::from(scan) * SCAN_SECONDS;
    psm.spectrum = key.clone();
    spectrum.key = key.clone();
    spectrum.source_id = Some(format!("scan={scan}"));
    spectrum.rt_seconds = Some(rt);
    let measured = measured.map(|peaks| PeakSpectrum {
        key,
        native_id: format!("scan={scan}"),
        rt_seconds: Some(rt),
        ..peaks
    });
    let identified = identified.map(|id| Identified { rt, ..id });
    (psm, spectrum, measured, identified)
}

fn fold(protein: usize, category: &str) -> f64 {
    match category {
        "treated" if UP.contains(&protein) => 4.0,
        "treated" if DOWN.contains(&protein) => 0.25,
        _ => 1.0,
    }
}

fn pick_weighted(rng: &mut Rng, weights: &[f64]) -> usize {
    let total: f64 = weights.iter().sum();
    let mut pick = rng.uniform() * total;
    for (index, weight) in weights.iter().enumerate() {
        if pick < *weight {
            return index;
        }
        pick -= weight;
    }
    weights.len() - 1
}

/// A random peptide of `protein`, if its digest has any.
fn random_peptide(rng: &mut Rng, proteome: &Proteome, protein: usize) -> Option<String> {
    let digest = &proteome.digests[protein];
    (!digest.is_empty()).then(|| digest[rng.index(digest.len())].clone())
}

fn form(rng: &mut Rng, sequence: &str) -> PeptideForm {
    let oxidized = sequence.find('M').filter(|_| rng.uniform() < 0.2);
    let modifications = oxidized
        .map(|position| Modification {
            site: ModSite::Residue(u32::try_from(position + 1).unwrap_or(1)),
            mass_delta: OXIDATION,
            name: Some("Oxidation".into()),
            accession: Some("UNIMOD:35".into()),
        })
        .into_iter()
        .collect();
    PeptideForm::new(sequence, modifications).unwrap_or_else(|_| PeptideForm::unmodified(sequence))
}

/// A correctly identified spectrum's precursor, for MS1.
struct Identified {
    form: PeptideForm,
    charge: u8,
    rt: f64,
}

/// One precursor's elution peak in a run.
struct Elution {
    mz: f64,
    charge: u8,
    apex: f64,
    height: f64,
    /// Relative isotope peak heights (M, M+1, M+2, M+3).
    isotopes: [f64; 4],
}

/// Elution peaks of the run's correctly identified precursors: centred on
/// the median of their identifications, height from the protein's weight.
fn elution(
    identified: &[Identified],
    expressed_in: &HashMap<&str, Vec<usize>>,
    weights: &[f64],
    loading: f64,
    rng: &mut Rng,
) -> Vec<Elution> {
    let mut by_precursor: HashMap<(String, u8), Vec<&Identified>> = HashMap::new();
    for id in identified {
        by_precursor
            .entry((id.form.proforma(), id.charge))
            .or_default()
            .push(id);
    }
    let mut keys: Vec<&(String, u8)> = by_precursor.keys().collect();
    keys.sort();
    keys.into_iter()
        .map(|key| {
            let ids = &by_precursor[key];
            let mut rts: Vec<f64> = ids.iter().map(|id| id.rt).collect();
            rts.sort_by(f64::total_cmp);
            let first = ids[0];
            let mz = precursor_mz(&first.form, first.charge);
            // Poisson isotope envelope, λ ≈ neutral mass / 1,800 Da.
            let lambda = (mz * f64::from(first.charge)) / 1_800.0;
            let mut isotopes = [0.0; 4];
            let mut term = 1.0;
            for (k, slot) in isotopes.iter_mut().enumerate() {
                if k > 0 {
                    term *= lambda / k as f64;
                }
                *slot = term;
            }
            Elution {
                mz: mz * (1.0 + 2e-6 * rng.normal()),
                charge: first.charge,
                apex: rts[rts.len() / 2],
                height: 1e6
                    * expressed_in
                        .get(first.form.sequence.as_str())
                        .map_or(0.0, |proteins| proteins.iter().map(|&p| weights[p]).sum())
                    * ionization(&first.form.sequence)
                    * loading
                    * (0.15 * rng.normal()).exp(),
                isotopes,
            }
        })
        .collect()
}

/// A fixed per-peptide factor (0.25×..4×), the same in every run, as
/// ionization efficiency is.
fn ionization(sequence: &str) -> f64 {
    let hash = sequence.bytes().fold(0xcbf2_9ce4_8422_2325_u64, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    });
    4_f64.powf((hash % 1_000) as f64 / 500.0 - 1.0)
}

/// The run's scans in time order: an MS1 scan before each MS/MS scan.
fn interleave_ms1(
    msms: Vec<PeakSpectrum>,
    elution: &[Elution],
    rng: &mut Rng,
) -> Vec<PeakSpectrum> {
    let mut scans = Vec::with_capacity(msms.len() * 2);
    for (index, spectrum) in msms.into_iter().enumerate() {
        let rt = spectrum.rt_seconds.unwrap_or_default() - SCAN_SECONDS / 2.0;
        let mut peaks: Vec<(f64, f32)> = (0..NOISE_PEAKS)
            .map(|_| {
                let intensity = 2e3 * (rng.normal()).exp();
                (400.0 + 1_200.0 * rng.uniform(), intensity as f32)
            })
            .collect();
        for peak in elution {
            let offset = (rt - peak.apex) / ELUTION_SIGMA;
            if offset.abs() > 5.0 {
                continue;
            }
            let height = peak.height * (-offset * offset / 2.0).exp();
            for (k, ratio) in peak.isotopes.iter().enumerate() {
                let mz = peak.mz + k as f64 * 1.003_355 / f64::from(peak.charge);
                peaks.push((mz, (height * ratio) as f32));
            }
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (mz, intensity) = peaks.into_iter().unzip();
        scans.push(PeakSpectrum {
            key: SpectrumKey::Index(u32::try_from(index).unwrap_or(u32::MAX)),
            native_id: format!("ms1={index}"),
            ms_level: 1,
            rt_seconds: Some(rt),
            precursor_mz: None,
            precursor_charge: None,
            mz,
            intensity,
            profile: false,
            parent: None,
            isolation: None,
            precursor_mobility: None,
            mobility: Vec::new(),
        });
        scans.push(spectrum);
    }
    scans
}

/// One spectrum: the reported PSM, the spectrum record, its peaks, and the
/// precursor it measured when it's a correct identification.
fn spectrum(
    rng: &mut Rng,
    proteome: &Proteome,
    weights: &[f64],
    run: &str,
    scan: u32,
) -> (Psm, Spectrum, Option<PeakSpectrum>, Option<Identified>) {
    let correct = rng.uniform() < CORRECT_SHARE;
    let charge: u8 = if rng.uniform() < 0.65 { 2 } else { 3 };
    let (reported, measured) = loop {
        let expressed = pick_weighted(rng, weights);
        let unseen = EXPRESSED + ISOFORMS + rng.index(TARGETS - EXPRESSED - ISOFORMS);
        let random = rng.index(TARGETS) + if rng.uniform() < 0.5 { TARGETS } else { 0 };
        let pair = if correct {
            random_peptide(rng, proteome, expressed).map(|p| (p.clone(), p))
        } else {
            random_peptide(rng, proteome, random).zip(random_peptide(rng, proteome, unseen))
        };
        if let Some(pair) = pair {
            break pair;
        }
    };
    let reported = form(rng, &reported);
    let measured = if correct {
        reported.clone()
    } else {
        form(rng, &measured)
    };
    let quality = 0.15 + 0.7 * rng.uniform();
    let (mz, intensity) = peaks(rng, &measured, charge, quality);
    let calculated = precursor_mz(&reported, charge);
    let experimental = calculated * (1.0 + 2e-6 * rng.normal());
    let (matched, explained) = explained(&reported, charge, &mz, &intensity);
    let hyperscore = 2.0 * matched as f64 + (1.0 + explained).ln() + rng.normal();
    let expect = 10_f64.powf(-(hyperscore - 12.0) / 4.0 + 0.3 * rng.normal());
    let accessions = proteome
        .evidence
        .get(&reported.sequence)
        .cloned()
        .unwrap_or_default();
    let key = SpectrumKey::Scan(scan);
    let psm = Psm {
        ms_run: run.into(),
        spectrum: key.clone(),
        rank: 1,
        peptide: reported,
        charge: i16::from(charge),
        calculated_mz: Some(calculated),
        experimental_mz: Some(experimental),
        pass_threshold: None,
        scores: vec![
            Score {
                name: "hyperscore".into(),
                accession: None,
                value: hyperscore,
            },
            Score {
                name: "expect".into(),
                accession: None,
                value: expect,
            },
            Score {
                name: "matched ions".into(),
                accession: None,
                value: matched as f64,
            },
        ],
        evidence: accessions
            .into_iter()
            .map(|accession| Evidence {
                accession,
                start: None,
                end: None,
                pre: None,
                post: None,
            })
            .collect(),
    };
    let spectrum = Spectrum {
        ms_run: run.into(),
        key: key.clone(),
        source_id: Some(format!("scan={scan}")),
        precursor_mz: Some(experimental),
        precursor_charge: Some(i16::from(charge)),
        rt_seconds: Some(f64::from(scan) * SCAN_SECONDS),
    };
    let measured_peaks = PeakSpectrum {
        key,
        native_id: format!("scan={scan}"),
        ms_level: 2,
        rt_seconds: Some(f64::from(scan) * SCAN_SECONDS),
        precursor_mz: Some(experimental),
        precursor_charge: Some(i16::from(charge)),
        mz,
        intensity,
        profile: false,
        parent: None,
        isolation: None,
        precursor_mobility: None,
        mobility: Vec::new(),
    };
    let identified = correct.then(|| Identified {
        form: measured.clone(),
        charge,
        rt: f64::from(scan) * SCAN_SECONDS,
    });
    (psm, spectrum, Some(measured_peaks), identified)
}

fn precursor_mz(form: &PeptideForm, charge: u8) -> f64 {
    let residues: f64 = form
        .sequence
        .chars()
        .filter_map(|r| {
            tuffold_core::chemistry::residue_mass(
                r,
                tuffold_core::chemistry::MassType::Monoisotopic,
            )
        })
        .sum();
    let deltas: f64 = form.modifications.iter().map(|m| m.mass_delta).sum();
    tuffold_core::chemistry::mz(
        residues + deltas + tuffold_core::chemistry::WATER_MONO,
        charge,
    )
}

/// Fragment peaks (each present with probability `quality`) plus noise,
/// sorted by m/z.
fn peaks(rng: &mut Rng, form: &PeptideForm, charge: u8, quality: f64) -> (Vec<f64>, Vec<f32>) {
    let mut peaks: Vec<(f64, f32)> = Vec::new();
    for ion in fragment_ions(form, Some(charge)).unwrap_or_default() {
        if rng.uniform() >= quality {
            continue;
        }
        let series = if ion.series == tuffold_core::fragments::IonSeries::Y {
            1.6
        } else {
            1.0
        };
        let height = (200.0 + 800.0 * rng.uniform()) * series / f64::from(ion.charge);
        peaks.push((ion.mz * (1.0 + 3e-6 * rng.normal()), height as f32));
    }
    let top = precursor_mz(form, charge) * f64::from(charge);
    for _ in 0..40 {
        let mz = 100.0 + (top - 100.0).max(400.0) * rng.uniform();
        peaks.push((mz, (40.0 + 260.0 * rng.uniform()) as f32));
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    peaks.into_iter().unzip()
}

/// Ions of `form` found in the peaks (20 ppm) and their summed intensity.
fn explained(form: &PeptideForm, charge: u8, mz: &[f64], intensity: &[f32]) -> (usize, f64) {
    let ions = fragment_ions(form, Some(charge)).unwrap_or_default();
    let matches = annotate(&ions, mz, intensity, Tolerance::Ppm(20.0));
    let total = matches.iter().map(|m| f64::from(intensity[m.peak])).sum();
    (matches.len(), total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trypsin_cuts_after_k_and_r_but_not_before_p() {
        // MAAAAAAK | RPAAAAAAAR (no cut in RP) | AAAAAAAK | AAK (too short)
        assert_eq!(
            digest("MAAAAAAKRPAAAAAAARAAAAAAAKAAK"),
            ["MAAAAAAK", "RPAAAAAAAR", "AAAAAAAK"]
        );
    }

    #[test]
    fn the_demo_is_deterministic_and_labelled_synthetic() {
        let first = generate(DEMO_SEED);
        let second = generate(DEMO_SEED);
        assert_eq!(first.batch, second.batch);
        assert_eq!(first.batch.psms.len(), 4 * SPECTRA_PER_RUN as usize);
        assert!(first.batch.software.name.contains("synthetic"));
        let decoys = first.batch.proteins.iter().filter(|p| p.is_decoy).count();
        assert_eq!(decoys, TARGETS);
        assert!(first.batch.psms.iter().all(|psm| !psm.evidence.is_empty()));
    }
}
