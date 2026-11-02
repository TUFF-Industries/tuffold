//! Peptide identifications: sequence, modifications, engine scores, and the
//! checks every one must pass.

use std::collections::HashMap;

use tuffold_core::chemistry::{PROTON, peptide_mass};
use tuffold_core::{ModSite, Modification, PeptideForm, Score};

use super::project::ModificationType;
use super::records::{Records, read_f32, read_i32};

pub(crate) const PEPTIDE_BLOCK: usize = 160;
const SEQUENCE: std::ops::Range<usize> = 4..34;
const SEQUENCE_END: u8 = 0x0a;
const TAIL: usize = 35;
const ANALYSIS: usize = 57;
const MOD_TYPES: [usize; 3] = [90, 94, 98];
const MOD_POSITIONS: [usize; 3] = [106, 110, 114];
const MORE_MOD_TYPES: usize = 102;
const MORE_MOD_POSITIONS: usize = 118;
const PROTEIN: usize = 122;
/// Unused integer and float slots.
const NONE: i32 = 0x7fff_ffff;
const ISOTOPE_SPACING: f64 = 1.003_355;
/// Longer "peptides" are a misread, not a peptide.
const MAX_PEPTIDE_LENGTH: usize = 1_000;

/// How a search engine's scores sit in a peptide identification. Rank 1 per
/// spectrum is the best `primary`; never Scaffold's own probability.
pub(crate) struct EngineScores {
    pub primary: usize,
    pub fields: &'static [(usize, &'static str)],
    /// Field names rest on value ranges only (docs/sf3-format.md).
    pub inferred: bool,
    /// Where the engine's expectation value comes from, when it has one.
    pub expectation: Option<Expectation>,
}

/// How an engine's expectation value sits in a peptide identification.
#[derive(Clone, Copy)]
pub(crate) enum Expectation {
    /// Mascot's, which follows from its ion score and identity threshold
    /// (offsets; not stored).
    Mascot { score: usize, identity: usize },
    /// Stored as −log10 of the value (MSFragger), under this name.
    NegLog10 { offset: usize, name: &'static str },
}

impl Expectation {
    fn score(self, block: &[u8]) -> Option<Score> {
        let (name, value) = match self {
            Self::Mascot { score, identity } => (
                "Mascot:expectation value",
                expectation(
                    read_f32(block, score)?.into(),
                    read_f32(block, identity)?.into(),
                ),
            ),
            Self::NegLog10 { offset, name } => {
                (name, 10f64.powf(-f64::from(read_f32(block, offset)?)))
            }
        };
        value.is_finite().then(|| Score {
            name: name.into(),
            accession: None,
            value,
        })
    }
}

/// Mascot's expectation value from its ion score `s` and identity threshold
/// `i`: the threshold is the score whose expectation is 0.05 (p < 0.05), and
/// expectation falls tenfold per 10 points, so E = 0.05 · 10^((i − s) / 10).
pub(crate) fn expectation(score: f64, identity: f64) -> f64 {
    0.05 * 10f64.powf((identity - score) / 10.0)
}

pub(crate) fn engine_scores(program: &str) -> Option<EngineScores> {
    let name = program.to_ascii_lowercase();
    Some(match name.as_str() {
        // Scaffold's "Mascot (Ion Score Only)" stores the same fields; only
        // its 66 (the ion score, not ion − identity) differs, and 66 isn't read.
        "mascot" | "mascot (ion score only)" => EngineScores {
            primary: 70,
            fields: &[
                (70, "Mascot:score"),
                (74, "Mascot:identity threshold"),
                (78, "Mascot:homology threshold"),
            ],
            inferred: false,
            expectation: Some(Expectation::Mascot {
                score: 70,
                identity: 74,
            }),
        },
        "maxquant" => EngineScores {
            primary: 66,
            fields: &[(66, "Andromeda:score"), (70, "MaxQuant:delta score")],
            inferred: false,
            expectation: None,
        },
        "spectrum mill" => EngineScores {
            primary: 70,
            fields: &[(70, "SpectrumMill:score"), (74, "SpectrumMill:%SPI")],
            inferred: false,
            expectation: None,
        },
        // Sequest and MSFragger are named by Scaffold's own score list
        // (docs/sf3-format.md): Sequest's 66 is Scaffold's discriminant
        // score (its verdict, left out), 70 deltaCn, 74 XCorr, 78 Sp.
        "sequest" => EngineScores {
            primary: 74,
            fields: &[
                (74, "SEQUEST:xcorr"),
                (70, "SEQUEST:deltacn"),
                (78, "SEQUEST:sp"),
            ],
            inferred: false,
            expectation: None,
        },
        // 70 hyperscore, 74 −log10 expectation value, 78 next-best
        // hyperscore; the values agree (Scaffold Q+S, PXD068782: 70 never
        // below 78, 74 rising with 70).
        "msfragger" => EngineScores {
            primary: 70,
            fields: &[(70, "MSFragger:hyperscore"), (78, "MSFragger:nextscore")],
            inferred: false,
            expectation: Some(Expectation::NegLog10 {
                offset: 74,
                name: "MSFragger:expectation value",
            }),
        },
        // Scaffold's score list says 66 and 74 −log expectation value, 70
        // hyperscore, but 70 ranks targets better *lower*, which no
        // hyperscore does, so the names stay unconfirmed. 66 separates
        // targets from decoys best and is ranked on.
        "x! tandem" | "x!tandem" => EngineScores {
            primary: 66,
            fields: &[
                (66, "X!Tandem score (sf3 offset 66)"),
                (70, "X!Tandem score (sf3 offset 70)"),
                (74, "X!Tandem score (sf3 offset 74)"),
            ],
            inferred: true,
            expectation: None,
        },
        _ => return None,
    })
}

/// One peptide identification, decoded but not yet checked against its
/// spectrum.
pub(crate) struct Candidate {
    pub decoded: Decoded,
    pub protein_record: i32,
    pub primary: f32,
    pub scores: Vec<Score>,
}

/// A candidate whose peptide is settled.
pub(crate) struct Resolved {
    pub peptide: PeptideForm,
    /// Monoisotopic neutral mass, when every residue is a standard one.
    pub mass: Option<f64>,
    pub protein_record: i32,
    pub primary: f32,
    pub scores: Vec<Score>,
}

/// A sequence and its stated modifications. Spectrum Mill writes a
/// modified residue in lower case without naming the modification;
/// those positions are `implied`.
pub(crate) struct Decoded {
    /// Upper-case.
    pub sequence: String,
    pub modifications: Vec<Modification>,
    pub implied: Vec<usize>,
}

/// Modifications a lower-case residue may stand for: the types the file
/// states on that residue elsewhere, by upper-case residue.
pub(crate) type ResidueOptions = HashMap<char, Vec<(String, f64)>>;

impl Decoded {
    fn stated_mass(&self) -> Option<f64> {
        peptide_mass(&self.sequence).map(|mass| {
            mass + self
                .modifications
                .iter()
                .map(|modification| modification.mass_delta)
                .sum::<f64>()
        })
    }

    /// Settle implied modifications: every lower-case occurrence of a
    /// residue takes the same one of its options, and the combination is
    /// the one, and only one, that fits the spectrum's [M+H]+ within
    /// `tolerance`. `None` when none or several fit.
    pub fn resolve(
        self,
        mh: f64,
        options: &ResidueOptions,
        tolerance: f64,
    ) -> Option<(PeptideForm, Option<f64>)> {
        let mass = self.stated_mass();
        if self.implied.is_empty() {
            let form = PeptideForm::new(&self.sequence, self.modifications).ok()?;
            return Some((form, mass));
        }
        let mass = mass?;
        let residues: Vec<char> = self.sequence.chars().collect();
        let mut letters: Vec<char> = self.implied.iter().map(|&index| residues[index]).collect();
        letters.sort_unstable();
        letters.dedup();
        let choices: Vec<&[(String, f64)]> = letters
            .iter()
            .map(|letter| options.get(letter).map_or(&[][..], Vec::as_slice))
            .collect();
        let combinations = choices
            .iter()
            .try_fold(1_usize, |product, choice| product.checked_mul(choice.len()))
            .unwrap_or(usize::MAX);
        if combinations == 0 || combinations > 64 {
            return None;
        }
        let mut fitting = Vec::new();
        for mut combination in 0..combinations {
            let picked: Vec<&(String, f64)> = choices
                .iter()
                .map(|choice| {
                    let pick = &choice[combination % choice.len()];
                    combination /= choice.len();
                    pick
                })
                .collect();
            let extra: f64 = self
                .implied
                .iter()
                .map(|&index| {
                    let letter = letters
                        .iter()
                        .position(|&l| l == residues[index])
                        .expect("listed");
                    picked[letter].1
                })
                .sum();
            if mass_error(mass + extra, mh).abs() <= tolerance {
                fitting.push((picked, mass + extra));
            }
        }
        let [(picked, total)] = fitting.as_slice() else {
            return None;
        };
        let mut modifications = self.modifications;
        for &index in &self.implied {
            let letter = letters
                .iter()
                .position(|&l| l == residues[index])
                .expect("listed");
            let (name, delta) = picked[letter];
            modifications.push(Modification {
                site: ModSite::Residue(u32::try_from(index + 1).ok()?),
                mass_delta: *delta,
                name: Some(name.clone()).filter(|name| !name.is_empty()),
                accession: None,
            });
        }
        let form = PeptideForm::new(&self.sequence, modifications).ok()?;
        Some((form, Some(*total)))
    }
}

/// Why a peptide identification could not be decoded. Any of these means
/// the layout is not what Tuffold knows, so the file is refused.
#[derive(Debug)]
pub(crate) struct Undecodable(pub String);

pub(crate) struct Peptides<'a> {
    pub tails: &'a Records,
    pub type_lists: &'a Records,
    pub position_lists: &'a Records,
    pub types: &'a HashMap<i32, ModificationType>,
}

impl Peptides<'_> {
    /// The sequence as stored: upper case, or lower case for a residue
    /// Spectrum Mill marks as modified.
    pub fn sequence(&self, block: &[u8]) -> Result<String, Undecodable> {
        let inline = &block[SEQUENCE];
        let end = inline
            .iter()
            .position(|&byte| byte == SEQUENCE_END)
            .unwrap_or(inline.len());
        let mut sequence = String::from_utf8_lossy(&inline[..end]).into_owned();
        let tail = read_i32(block, TAIL).unwrap_or(0);
        if tail != 0 {
            let rest = self
                .tails
                .string(tail)
                .ok_or_else(|| Undecodable(format!("sequence tail {tail} is missing")))?;
            sequence.push_str(&rest);
        }
        if sequence.is_empty()
            || sequence.len() > MAX_PEPTIDE_LENGTH
            || !sequence.bytes().all(|byte| byte.is_ascii_alphabetic())
        {
            return Err(Undecodable(format!(
                "sequence {sequence:?} is not residues"
            )));
        }
        Ok(sequence)
    }

    /// (0-based position, type id) pairs; a later one at the same position
    /// replaces an earlier one (Scaffold lists e.g. carbamidomethyl, then
    /// the NEM that replaced it).
    fn modification_ids(&self, block: &[u8]) -> Result<Vec<(i32, i32)>, Undecodable> {
        let mut pairs: Vec<(i32, i32)> = MOD_TYPES
            .iter()
            .zip(MOD_POSITIONS)
            .filter_map(|(&kind, position)| {
                let kind = read_i32(block, kind)?;
                (kind != NONE).then(|| (read_i32(block, position).unwrap_or(NONE), kind))
            })
            .collect();
        let (more_types, more_positions) = (
            read_i32(block, MORE_MOD_TYPES).unwrap_or(0),
            read_i32(block, MORE_MOD_POSITIONS).unwrap_or(0),
        );
        if more_types != 0 || more_positions != 0 {
            let types = self.type_lists.int_list(more_types);
            let positions = self.position_lists.int_list(more_positions);
            match (types, positions) {
                (Some(types), Some(positions)) if types.len() == positions.len() => {
                    pairs.extend(positions.into_iter().zip(types));
                }
                _ => {
                    return Err(Undecodable(format!(
                        "modification lists {more_types}/{more_positions} are missing or unequal"
                    )));
                }
            }
        }
        let mut by_position: Vec<(i32, i32)> = Vec::with_capacity(pairs.len());
        for (position, kind) in pairs {
            match by_position.iter_mut().find(|(known, _)| *known == position) {
                Some(slot) => slot.1 = kind,
                None => by_position.push((position, kind)),
            }
        }
        Ok(by_position)
    }

    pub fn decode(&self, block: &[u8]) -> Result<Decoded, Undecodable> {
        let stored = self.sequence(block)?;
        let sequence = stored.to_ascii_uppercase();
        let residues: Vec<char> = sequence.chars().collect();
        let mut modifications = Vec::new();
        let mut stated = Vec::new();
        for (position, kind) in self.modification_ids(block)? {
            let modification = self
                .types
                .get(&kind)
                .ok_or_else(|| Undecodable(format!("unknown modification type {kind}")))?;
            let index = usize::try_from(position)
                .ok()
                .filter(|&index| index < residues.len())
                .ok_or_else(|| {
                    Undecodable(format!(
                        "{} at {position} is outside {sequence}",
                        modification.name
                    ))
                })?;
            let site = if modification.residues.contains(residues[index]) {
                stated.push(index);
                ModSite::Residue(u32::try_from(index + 1).expect("short peptide"))
            } else if index == 0 && modification.residues.contains('n') {
                ModSite::NTerm
            } else if index + 1 == residues.len() && modification.residues.contains('c') {
                ModSite::CTerm
            } else {
                return Err(Undecodable(format!(
                    "{} cannot sit on {} of {sequence}",
                    modification.name, residues[index]
                )));
            };
            modifications.push(Modification {
                site,
                mass_delta: modification.mono_mass,
                name: Some(modification.name.clone()).filter(|name| !name.is_empty()),
                accession: None,
            });
        }
        let implied = stored
            .char_indices()
            .filter(|(index, residue)| residue.is_ascii_lowercase() && !stated.contains(index))
            .map(|(index, _)| index)
            .collect();
        Ok(Decoded {
            sequence,
            modifications,
            implied,
        })
    }

    /// The modification types the file states on each residue.
    pub fn residue_options(&self, records: &Records) -> ResidueOptions {
        let mut options = ResidueOptions::new();
        for id in records.ids() {
            let Some(block) = records.get(id) else {
                continue;
            };
            let (Ok(sequence), Ok(ids)) = (self.sequence(block), self.modification_ids(block))
            else {
                continue;
            };
            let residues: Vec<char> = sequence.to_ascii_uppercase().chars().collect();
            for (position, kind) in ids {
                let (Some(&residue), Some(modification)) = (
                    usize::try_from(position)
                        .ok()
                        .and_then(|index| residues.get(index)),
                    self.types.get(&kind),
                ) else {
                    continue;
                };
                if modification.residues.contains(residue) {
                    let known = options.entry(residue).or_default();
                    if !known.iter().any(|(name, _)| *name == modification.name) {
                        known.push((modification.name.clone(), modification.mono_mass));
                    }
                }
            }
        }
        options
    }

    pub fn candidate(&self, block: &[u8], engine: &EngineScores) -> Result<Candidate, Undecodable> {
        let decoded = self.decode(block)?;
        let primary = read_f32(block, engine.primary).unwrap_or(f32::NAN);
        let scores = engine
            .fields
            .iter()
            .filter_map(|&(offset, name)| {
                let value = read_f32(block, offset)?;
                value.is_finite().then(|| Score {
                    name: name.to_string(),
                    accession: None,
                    value: f64::from(value),
                })
            })
            .collect::<Vec<Score>>();
        let scores = scores
            .into_iter()
            .chain(engine.expectation.and_then(|kind| kind.score(block)))
            .collect();
        Ok(Candidate {
            decoded,
            protein_record: read_i32(block, PROTEIN).unwrap_or(0),
            primary,
            scores,
        })
    }
}

/// The analysis id a peptide identification belongs to.
pub(crate) fn analysis_of(block: &[u8]) -> i32 {
    read_i32(block, ANALYSIS).unwrap_or(0)
}

/// Neutral mass error in daltons against the spectrum's [M+H]+, allowing
/// the isotope errors engines search (−1..+3).
pub(crate) fn mass_error(mass: f64, mh: f64) -> f64 {
    let observed = mh - PROTON;
    (-1..=3)
        .map(|isotope| observed - mass - f64::from(isotope) * ISOTOPE_SPACING)
        .min_by(|a, b| a.abs().total_cmp(&b.abs()))
        .expect("five isotopes")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mascot_expectation_is_005_at_the_identity_threshold() {
        // At the threshold E = 0.05; 10 points above it, 0.005; 20 below, 5.
        assert!((expectation(35.0, 35.0) - 0.05).abs() < 1e-12);
        assert!((expectation(45.0, 35.0) - 0.005).abs() < 1e-12);
        assert!((expectation(15.0, 35.0) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn mass_errors_allow_isotope_offsets() {
        let mass = 1000.0;
        assert!(mass_error(mass, mass + PROTON).abs() < 1e-9);
        assert!(mass_error(mass, mass + PROTON + ISOTOPE_SPACING).abs() < 1e-9);
        assert!(
            (mass_error(mass, mass + PROTON + 57.02146) - 57.02146 + 3.0 * ISOTOPE_SPACING).abs()
                < 1e-9
        );
    }

    fn decoded(sequence: &str, implied: &[usize]) -> Decoded {
        Decoded {
            sequence: sequence.into(),
            modifications: Vec::new(),
            implied: implied.to_vec(),
        }
    }

    #[test]
    fn lower_case_residues_take_the_one_modification_that_fits() {
        let options: ResidueOptions = [(
            'C',
            vec![
                ("Carbamidomethyl".to_string(), 57.021_464),
                ("Nethylmaleimide".to_string(), 125.047_679),
            ],
        )]
        .into();
        let mass = peptide_mass("PEPCK").unwrap();
        let mh = mass + 125.047_679 + PROTON;
        let (form, total) = decoded("PEPCK", &[3]).resolve(mh, &options, 1.0).unwrap();
        assert_eq!(form.proforma(), "PEPC[+125.0477]K");
        assert!((total.unwrap() - mass - 125.047_679).abs() < 1e-9);
        // Nothing fits a mass 30 Da off; a residue without options can't be settled.
        assert!(
            decoded("PEPCK", &[3])
                .resolve(mh + 30.0, &options, 1.0)
                .is_none()
        );
        assert!(decoded("PEPMK", &[3]).resolve(mh, &options, 1.0).is_none());
        // Two options within the tolerance of each other are ambiguous.
        let close: ResidueOptions = [('C', vec![("A".into(), 57.0), ("B".into(), 57.5)])].into();
        let mh = mass + 57.2 + PROTON;
        assert!(decoded("PEPCK", &[3]).resolve(mh, &close, 1.0).is_none());
    }

    #[test]
    fn known_engines_have_scores_and_others_do_not() {
        assert_eq!(engine_scores("Mascot").unwrap().primary, 70);
        assert_eq!(
            engine_scores("Mascot (Ion Score Only)").unwrap().primary,
            70
        );
        assert!(engine_scores("X! Tandem").unwrap().inferred);
        // Scaffold's list for Sequest: Discriminant Score, deltaCn, XCorr,
        // Sp (slots 66, 70, 74, 78); the first is Scaffold's, not Sequest's.
        let sequest = engine_scores("Sequest").unwrap();
        assert!(!sequest.inferred);
        assert_eq!(sequest.primary, 74);
        assert_eq!(
            sequest.fields,
            &[
                (74, "SEQUEST:xcorr"),
                (70, "SEQUEST:deltacn"),
                (78, "SEQUEST:sp")
            ]
        );
        assert!(engine_scores("Byonic").is_none());
    }
}
