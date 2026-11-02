//! How well the decoy-alanine FLR estimate tracks the true false localization
//! rate on seeded synthetic phosphopeptide spectra (the generator of
//! `localize.rs`'s calibration test, over many seeds). Prints, asserts
//! nothing:
//!
//! `cargo test -p tuffold-validate --release --test localization_flr -- --ignored --nocapture`
//!
//! Environment: `SEEDS` (20), `SPECTRA` per seed (3000), `NOISE` = `union`
//! (each distinct ion of every candidate isoform and the Ala isoform at
//! `NOISE_RATE`, 4%), `ladder` (each isoform's ions drawn separately, so ions
//! shared by several wrong isoforms get several chances) or `uniform`
//! (`UNIFORM` random peaks per 100 Th, scored at 0.5 Da), and `DIAG`, the cut
//! for the breakdown by distance (0.9).

use tuffold_core::fragments::{Tolerance, fragment_ions};
use tuffold_core::peptide::{ModSite, Modification, PeptideForm};
use tuffold_validate::localize::{
    LocalizationInput, ScoredLocalization, SpectrumLocalization, localize_spectrum,
};

const TEMPLATES: [(&str, &[u32]); 5] = [
    ("ASGPTEYK", &[2, 5, 7]),
    ("SAGPTVSK", &[1, 5, 7]),
    ("GAPSTLDK", &[4, 5]),
    ("ALSTGPEK", &[3, 4]),
    ("SVAGTYPK", &[1, 5, 6]),
];
const CUTS: [f64; 4] = [0.75, 0.90, 0.95, 0.99];
/// Distances 1, 2, and 3 or more residues.
const STRATA: usize = 3;

fn env_or<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn phospho_at(sequence: &str, pos: u32) -> PeptideForm {
    let modification = Modification {
        site: ModSite::Residue(pos),
        mass_delta: 79.966_331,
        name: Some("Phospho".into()),
        accession: None,
    };
    PeptideForm::new(sequence, vec![modification]).expect("valid template")
}

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() as f64) / (u64::MAX as f64)
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Noise {
    Union,
    Ladder,
    Uniform,
}

/// One spectrum: the true isoform's ions kept at a random 5–70%, plus noise.
fn spectrum(
    rng: &mut Rng,
    sequence: &str,
    sites: &[u32],
    true_pos: u32,
    noise: Noise,
) -> Vec<(f64, f32)> {
    let rate: f64 = env_or("NOISE_RATE", 0.04);
    let survival = 0.05 + 0.65 * rng.unit();
    let mut peaks = Vec::new();
    for ion in fragment_ions(&phospho_at(sequence, true_pos), Some(2)).expect("fragments") {
        if rng.unit() < survival {
            peaks.push((ion.mz, 200.0 + 800.0 * rng.unit() as f32));
        }
    }
    let alanine = u32::try_from(sequence.find('A').expect("template has Ala")).expect("short") + 1;
    let isoforms: Vec<Vec<f64>> = sites
        .iter()
        .chain(std::iter::once(&alanine))
        .map(|&pos| {
            fragment_ions(&phospho_at(sequence, pos), Some(2))
                .expect("fragments")
                .iter()
                .map(|ion| ion.mz)
                .collect()
        })
        .collect();
    match noise {
        Noise::Ladder => {
            for mz in isoforms.iter().flatten() {
                if rng.unit() < rate {
                    peaks.push((*mz, 50.0 + 300.0 * rng.unit() as f32));
                }
            }
        }
        Noise::Union => {
            let mut distinct: Vec<f64> = isoforms.into_iter().flatten().collect();
            distinct.sort_by(f64::total_cmp);
            distinct.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
            for mz in distinct {
                if rng.unit() < rate {
                    peaks.push((mz, 50.0 + 300.0 * rng.unit() as f32));
                }
            }
        }
        Noise::Uniform => {
            let per_100: f64 = env_or("UNIFORM", 10.0);
            for _ in 0..(per_100 * 10.0) as usize {
                peaks.push((
                    100.0 + 1000.0 * rng.unit(),
                    50.0 + 300.0 * rng.unit() as f32,
                ));
            }
        }
    }
    if peaks.is_empty() {
        peaks.push((150.0, 100.0));
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    peaks
}

/// Distance stratum of `pos` from the nearest assigned target site.
fn stratum(scored: &ScoredLocalization, pos: u32) -> usize {
    let distance = scored
        .target_sites()
        .filter(|s| s.assigned)
        .map(|s| s.position.abs_diff(pos))
        .min()
        .unwrap_or(u32::MAX);
    distance.clamp(1, STRATA as u32) as usize - 1
}

#[derive(Default)]
struct Tally {
    targets_passing: [usize; 4],
    wrong_passing: [usize; 4],
    decoys_passing: [usize; 4],
    false_opportunities: usize,
    decoys_tested: usize,
    /// Unassigned target sites and decoys tested, by distance to the assigned site.
    opportunities_by_stratum: [usize; STRATA],
    decoys_by_stratum: [usize; STRATA],
    decoys_passing_by_stratum: [[usize; STRATA]; 4],
    /// By distance to the true site (0..8), at the `DIAG` cut: [tested, passing].
    wrong_by_distance: [[usize; 2]; 8],
    decoy_by_distance: [[usize; 2]; 8],
}

impl Tally {
    fn add(&mut self, scored: &ScoredLocalization, true_pos: u32, diag: f64) {
        self.false_opportunities += scored.false_target_sites as usize;
        for site in &scored.sites {
            let distance = (site.position.abs_diff(true_pos) as usize).min(7);
            let diag_pass = usize::from(site.probability >= diag);
            if site.is_decoy {
                self.decoys_tested += 1;
                self.decoy_by_distance[distance][0] += 1;
                self.decoy_by_distance[distance][1] += diag_pass;
                let s = stratum(scored, site.position);
                self.decoys_by_stratum[s] += 1;
                for (c, &cut) in CUTS.iter().enumerate() {
                    let pass = usize::from(site.probability >= cut);
                    self.decoys_passing[c] += pass;
                    self.decoys_passing_by_stratum[c][s] += pass;
                }
                continue;
            }
            let wrong = site.position != true_pos;
            if wrong {
                self.wrong_by_distance[distance][0] += 1;
                self.wrong_by_distance[distance][1] += diag_pass;
            }
            if !site.assigned {
                self.opportunities_by_stratum[stratum(scored, site.position)] += 1;
            }
            for (c, &cut) in CUTS.iter().enumerate() {
                if site.probability >= cut {
                    self.targets_passing[c] += 1;
                    self.wrong_passing[c] += usize::from(wrong);
                }
            }
        }
    }

    fn report(&self) {
        let pct = |num: f64, den: usize| 100.0 * num / den.max(1) as f64;
        println!(
            "false target opportunities {}, decoy sites tested {}",
            self.false_opportunities, self.decoys_tested
        );
        println!("cut    targets  wrong  true FLR  decoys  estimate  stratified");
        for (c, &cut) in CUTS.iter().enumerate() {
            // The estimate without the +1, to compare expectations.
            let estimate = self.false_opportunities as f64 * self.decoys_passing[c] as f64
                / self.decoys_tested.max(1) as f64;
            let pooled = self.decoys_passing[c] as f64 / self.decoys_tested.max(1) as f64;
            let stratified: f64 = (0..STRATA)
                .map(|s| {
                    let rate = if self.decoys_by_stratum[s] > 0 {
                        self.decoys_passing_by_stratum[c][s] as f64
                            / self.decoys_by_stratum[s] as f64
                    } else {
                        pooled
                    };
                    self.opportunities_by_stratum[s] as f64 * rate
                })
                .sum();
            println!(
                "{cut:.2}  {:>8} {:>6}  {:>7.3}%  {:>6}  {:>7.3}%  {:>9.3}%",
                self.targets_passing[c],
                self.wrong_passing[c],
                pct(self.wrong_passing[c] as f64, self.targets_passing[c]),
                self.decoys_passing[c],
                pct(estimate, self.targets_passing[c]),
                pct(stratified, self.targets_passing[c]),
            );
        }
        println!("distance to true site: wrong targets tested / passing, decoys tested / passing");
        for d in 1..8 {
            let (w, x) = (self.wrong_by_distance[d], self.decoy_by_distance[d]);
            if w[0] + x[0] > 0 {
                println!("  {d}  {:>7} {:>5}   {:>7} {:>5}", w[0], w[1], x[0], x[1]);
            }
        }
    }
}

#[test]
#[ignore = "prints a seeded sweep; run by hand"]
fn decoy_alanine_flr_sweep() {
    let seeds: u64 = env_or("SEEDS", 20);
    let spectra: i64 = env_or("SPECTRA", 3000);
    let diag: f64 = env_or("DIAG", 0.9);
    let noise = match std::env::var("NOISE").as_deref() {
        Ok("ladder") => Noise::Ladder,
        Ok("uniform") => Noise::Uniform,
        _ => Noise::Union,
    };
    let tolerance = if noise == Noise::Uniform {
        Tolerance::Da(0.5)
    } else {
        Tolerance::Ppm(20.0)
    };
    let mut tally = Tally::default();
    for seed in 0..seeds {
        let mut rng = Rng(0x5EED_CAFE_1234_5678 ^ (seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1));
        for id in 0..spectra {
            let (sequence, sites) = TEMPLATES[id as usize % TEMPLATES.len()];
            let true_pos = sites[rng.next_u64() as usize % sites.len()];
            let reported = phospho_at(sequence, sites[rng.next_u64() as usize % sites.len()]);
            let peaks = spectrum(&mut rng, sequence, sites, true_pos, noise);
            let mz: Vec<f64> = peaks.iter().map(|p| p.0).collect();
            let intensity: Vec<f32> = peaks.iter().map(|p| p.1).collect();
            let input = LocalizationInput {
                spectrum_id: id,
                form: &reported,
                charge: Some(2),
                mz: &mz,
                intensity: &intensity,
            };
            if let SpectrumLocalization::Scored(scored) = localize_spectrum(&input, tolerance) {
                tally.add(&scored, true_pos, diag);
            }
        }
    }
    println!("{seeds} seeds x {spectra} spectra");
    tally.report();
}
