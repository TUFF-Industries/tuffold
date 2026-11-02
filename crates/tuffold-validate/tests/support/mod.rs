//! A seeded synthetic search with known ground truth.
//!
//! Each spectrum is correct with probability 1/2 (its true peptide comes from
//! one of [`PRESENT`] proteins) or incorrect (a random target or decoy peptide
//! from all [`PROTEINS`]). Engines score independently: three features that
//! are each weakly higher for correct matches. `tools/crosscheck/synthetic.py`
//! ports this generator draw for draw to run mokapot on the same data, so the
//! order of random draws below must not change.

// Shared by several test binaries; each uses only part of it.
#![allow(dead_code)]

use tuffold_validate::{Input, PeptideRow, ProteinRow, PsmRow};

pub const PROTEINS: u64 = 2000;
pub const PRESENT: u64 = 200;
pub const PEPTIDES_PER_PROTEIN: u64 = 5;
const DECOY_OFFSET: i64 = 1_000_000;
const LETTERS: &[u8] = b"ACDEFGHILMNQSTVWY";

pub struct SplitMix(u64);

impl SplitMix {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn uniform(&mut self) -> f64 {
        (self.next() >> 11) as f64 * (1.0 / (1_u64 << 53) as f64)
    }

    pub fn normal(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * (1.0 - u1).ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    pub fn index(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

pub struct Synthetic {
    pub input: Input,
    /// The true peptide of each spectrum (index = spectrum id - 1).
    pub true_peptides: Vec<i64>,
}

fn sequence(body: u64, decoy: bool) -> String {
    let mut digits = Vec::new();
    let mut rest = body;
    loop {
        digits.push(LETTERS[(rest % LETTERS.len() as u64) as usize]);
        rest /= LETTERS.len() as u64;
        if rest == 0 {
            break;
        }
    }
    while digits.len() < 6 {
        digits.push(b'A');
    }
    if decoy {
        digits.reverse();
    }
    let mut text = String::from_utf8(digits).expect("ASCII");
    text.push('K');
    text
}

fn peptide_id(protein: u64, k: u64, decoy: bool) -> i64 {
    let id = i64::try_from(protein * PEPTIDES_PER_PROTEIN + k).expect("small") + 1;
    if decoy { DECOY_OFFSET + id } else { id }
}

pub fn generate(seed: u64, spectra: u64, engines: usize) -> Synthetic {
    generate_with(seed, spectra, engines, 0.0)
}

/// Like [`generate`], but a share `near_miss` of incorrect spectra match the
/// same wrong target peptide in every engine, scoring a little higher than
/// random matches (shift 0.5) — the correlated error decoys cannot model.
/// With `near_miss == 0` the draws are exactly [`generate`]'s.
pub fn generate_with(seed: u64, spectra: u64, engines: usize, near_miss: f64) -> Synthetic {
    generate_partial(seed, spectra, engines, near_miss, PEPTIDES_PER_PROTEIN)
}

/// Like [`generate_with`], but only the first `present_peptides` peptides of
/// each present protein are in the sample, so wrong matches can land on
/// absent peptides of present proteins (which sibling evidence lifts). With
/// all [`PEPTIDES_PER_PROTEIN`] the draws are exactly [`generate_with`]'s.
pub fn generate_partial(
    seed: u64,
    spectra: u64,
    engines: usize,
    near_miss: f64,
    present_peptides: u64,
) -> Synthetic {
    let mut rng = SplitMix::new(seed);
    let mut input = Input {
        engines: (0..engines).map(|e| format!("Engine {e}")).collect(),
        score_names: vec!["f1".into(), "f2".into(), "expect".into()],
        ..Input::default()
    };
    for protein in 0..PROTEINS {
        for decoy in [false, true] {
            let id = i64::try_from(protein).expect("small") + 1;
            input.proteins.push(ProteinRow {
                id: if decoy { DECOY_OFFSET + id } else { id },
                accession: format!("{}sp|P{protein:05}|SYN", if decoy { "rev_" } else { "" }),
                is_decoy: decoy,
            });
            for k in 0..PEPTIDES_PER_PROTEIN {
                let peptide = peptide_id(protein, k, decoy);
                input.peptides.push(PeptideRow {
                    id: peptide,
                    sequence: sequence(protein * PEPTIDES_PER_PROTEIN + k, decoy),
                    modifications: 0,
                    is_decoy: decoy,
                });
                let protein_id = if decoy { DECOY_OFFSET + id } else { id };
                input.peptide_proteins.push((peptide, protein_id));
            }
        }
    }
    let mut true_peptides = Vec::new();
    let mut psm_id = 0;
    for spectrum in 1..=spectra {
        let spectrum_id = i64::try_from(spectrum).expect("small");
        let correct = rng.uniform() < 0.5;
        let true_peptide = peptide_id(rng.index(PRESENT), rng.index(present_peptides), false);
        true_peptides.push(true_peptide);
        let shared = (near_miss > 0.0 && !correct && rng.uniform() < near_miss)
            .then(|| peptide_id(rng.index(PROTEINS), rng.index(PEPTIDES_PER_PROTEIN), false));
        for engine in 0..engines {
            let peptide = if correct {
                true_peptide
            } else if let Some(shared) = shared {
                shared
            } else {
                let decoy = rng.uniform() < 0.5;
                peptide_id(rng.index(PROTEINS), rng.index(PEPTIDES_PER_PROTEIN), decoy)
            };
            let shift = if correct {
                1.0
            } else if shared.is_some() {
                0.5
            } else {
                0.0
            };
            let f1 = rng.normal() + 1.6 * shift;
            let f2 = rng.normal() + 1.2 * shift;
            let expect = 10_f64.powf(-(rng.normal() + 1.0 * shift + 2.0));
            psm_id += 1;
            input.psms.push(PsmRow {
                id: psm_id,
                spectrum_id,
                peptide_id: peptide,
                engine,
                charge: 2,
                calculated_mz: None,
                experimental_mz: None,
                run_id: 1,
                retention_time: None,
                scores: vec![(0, f1), (1, f2), (2, expect)],
            });
        }
    }
    Synthetic {
        input,
        true_peptides,
    }
}

/// Renames every decoy protein so it names no target (as TPP's numbered
/// `reverse-5287` or a Scaffold file's `tr|…|…-DECOY` beside targets named
/// `A0QSM5`); peptide sequences, and so sequence pairing, are untouched.
pub fn with_unnamed_decoys(mut synthetic: Synthetic) -> Synthetic {
    for protein in &mut synthetic.input.proteins {
        if protein.is_decoy {
            protein.accession = format!("reverse-{}", protein.id);
        }
    }
    synthetic
}

/// Rotates every decoy peptide's body (same length and composition, still
/// ending in K), so no reversal of it names a target: decoys pair with
/// nothing by sequence, as shuffled decoys don't.
pub fn with_unpaired_decoy_sequences(mut synthetic: Synthetic) -> Synthetic {
    for peptide in synthetic.input.peptides.iter_mut().filter(|p| p.is_decoy) {
        let (body, end) = peptide.sequence.split_at(peptide.sequence.len() - 1);
        peptide.sequence = format!("{}{}{end}", &body[2..], &body[..2]);
    }
    synthetic
}

/// Keeps only one engine's PSMs.
pub fn only_engine(mut synthetic: Synthetic, engine: usize) -> Synthetic {
    synthetic.input.psms.retain(|psm| psm.engine == engine);
    synthetic
}

/// Each peptide's true elution time in [`with_retention_times`]: its
/// composition's part plus its own deviation (what a perfect predictor would
/// know).
pub fn elution_times(synthetic: &Synthetic, seed: u64) -> std::collections::HashMap<i64, f64> {
    synthetic
        .input
        .peptides
        .iter()
        .map(|peptide| {
            let hydrophobicity: f64 = peptide
                .sequence
                .bytes()
                .map(|residue| {
                    let index = LETTERS.iter().position(|&l| l == residue).unwrap_or(0);
                    ((index * 7) % LETTERS.len()) as f64 / (LETTERS.len() - 1) as f64
                })
                .sum::<f64>()
                / peptide.sequence.len() as f64;
            let mut own =
                SplitMix::new(seed ^ 0x5254_5254 ^ (peptide.id as u64).wrapping_mul(0x9E37_79B9));
            (
                peptide.id,
                300.0 + 3000.0 * hydrophobicity + 150.0 * own.normal(),
            )
        })
        .collect()
}

/// Seconds of the simulated gradient.
pub const GRADIENT: f64 = 4000.0;

/// Gives every spectrum a retention time, from a stream of its own so the
/// search's draws stay as they were. A peptide elutes at a time set mostly by
/// its composition (each letter's fixed hydrophobicity) plus a deviation of
/// its own that composition can't explain (σ 150 s); a correct spectrum is
/// taken within σ 20 s of its true peptide's time, an incorrect one at any
/// time of the gradient. The run is one (id 1).
pub fn with_retention_times(mut synthetic: Synthetic, seed: u64) -> Synthetic {
    let elution = elution_times(&synthetic, seed);
    let elution = |peptide: i64| elution[&peptide];
    // Correct spectra are matched to their true peptide (every engine).
    let truth = &synthetic.true_peptides;
    let correct: std::collections::HashSet<i64> = synthetic
        .input
        .psms
        .iter()
        .filter(|psm| {
            truth[usize::try_from(psm.spectrum_id - 1).expect("ids from 1")] == psm.peptide_id
        })
        .map(|psm| psm.spectrum_id)
        .collect();
    let mut rng = SplitMix::new(seed ^ 0x0052_5449_4D45);
    let times: Vec<f64> = (1..=truth.len())
        .map(|spectrum| {
            let spectrum_id = i64::try_from(spectrum).expect("small");
            if correct.contains(&spectrum_id) {
                elution(truth[spectrum - 1]) + 20.0 * rng.normal()
            } else {
                GRADIENT * rng.uniform()
            }
        })
        .collect();
    for psm in &mut synthetic.input.psms {
        psm.retention_time = Some(times[usize::try_from(psm.spectrum_id - 1).expect("ids from 1")]);
    }
    synthetic
}
