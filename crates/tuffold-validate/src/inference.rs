//! Protein inference: parsimony grouping and picked protein-group FDR.
//!
//! Only peptides passing the peptide FDR threshold take part, so the protein
//! list follows the threshold the user chose. A target peptide is evidence
//! only for target proteins, a decoy peptide only for decoys.
//!
//! 1. Proteins with identical peptide sets form one entity (indistinguishable
//!    proteins; the leader is a Swiss-Prot `sp|` entry when there is one,
//!    then the first accession).
//! 2. An entity whose peptides are a strict subset of another's is not a
//!    group; it is listed as a subset member.
//! 3. Shared peptides go to one group as razor peptides: the group with the
//!    most peptides, then the best peptide score, then the first leader.
//!    A group left with no unique or razor peptide is explained by the
//!    others and becomes a subset member too (Occam's razor).
//! 4. A group scores its best unique or razor peptide (Savitski et al. 2015),
//!    so shared evidence is never counted twice.
//! 5. Picked FDR: a decoy group is paired with the target group holding the
//!    accession it was derived from (`rev_X` with `X`); only the better of the
//!    pair is kept (ties to the decoy), then q-values follow target-decoy
//!    competition over the kept groups. Groups that lost their pair have no
//!    q-value and are never reported.

use std::collections::{BTreeMap, HashMap, HashSet};

use tuffold_core::decoy::target_accession;

use crate::fdr::confidence;
use crate::input::ProteinRow;
use crate::notes::{Note, NoteCode};
use crate::settings::{Level, Thresholds};

/// A peptide as protein inference sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct PeptideEvidence<'a> {
    pub id: i64,
    pub sequence: &'a str,
    pub is_decoy: bool,
    pub score: f64,
    pub q_value: f64,
    pub pep: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MemberRole {
    Leader,
    /// Same peptides as the leader; indistinguishable from it.
    Equivalent,
    /// Every peptide also belongs to this group (or is explained by it).
    Subset,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PeptideRole {
    /// Maps to this group only.
    Unique,
    /// Shared, and credited to this group.
    Razor,
    /// Shared, and credited to another group.
    Shared,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ProteinGroup {
    pub leader: i64,
    /// Every member, the leader first.
    pub members: Vec<(i64, MemberRole)>,
    pub peptides: Vec<(i64, PeptideRole)>,
    pub is_decoy: bool,
    pub score: f64,
    /// `None` when the group lost to its target-decoy partner.
    pub q_value: Option<f64>,
    /// Posterior error probability among the groups that competed (`None`
    /// for those that lost, like the q-value).
    pub pep: Option<f64>,
    /// Distinct sequences among unique and razor peptides.
    pub distinct_sequences: u32,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ProteinInference {
    pub groups: Vec<ProteinGroup>,
    pub notes: Vec<Note>,
}

/// Indistinguishable proteins and their (sorted, local) peptide indices.
struct Entity {
    peptides: Vec<usize>,
    proteins: Vec<i64>,
}

pub fn infer_proteins(
    peptides: &[PeptideEvidence<'_>],
    peptide_proteins: &[(i64, i64)],
    proteins: &[ProteinRow],
    thresholds: &Thresholds,
) -> ProteinInference {
    let mut passing: Vec<&PeptideEvidence<'_>> = peptides
        .iter()
        .filter(|peptide| {
            thresholds.passes(Level::Peptide, Some(peptide.q_value), Some(peptide.pep))
        })
        .collect();
    passing.sort_by_key(|peptide| peptide.id);
    let mut notes = Vec::new();
    if passing.is_empty() {
        notes.push(Note::new(
            NoteCode::NoPassingPeptides,
            format!(
                "No peptide passes {}, so no protein groups were formed.",
                thresholds.describe(Level::Peptide)
            ),
        ));
        return ProteinInference {
            groups: Vec::new(),
            notes,
        };
    }
    let catalog: HashMap<i64, &ProteinRow> = proteins.iter().map(|p| (p.id, p)).collect();
    let sets = protein_sets(&passing, peptide_proteins, &catalog);
    let covered: HashSet<usize> = sets.values().flatten().copied().collect();
    if covered.len() < passing.len() {
        notes.push(Note::new(
            NoteCode::PeptidesWithoutProteins,
            format!(
                "{} passing peptides map to no protein of their own kind (target or decoy) \
                 and are left out of protein groups.",
                passing.len() - covered.len()
            ),
        ));
    }
    let entities = entities(sets, &catalog);
    let groups = build_groups(&entities, &passing, &catalog);
    let (groups, unpaired) = picked_q_values(groups, &catalog);
    notes.extend(unpaired);
    ProteinInference { groups, notes }
}

/// Protein id -> sorted local indices of the passing peptides it explains.
fn protein_sets(
    passing: &[&PeptideEvidence<'_>],
    peptide_proteins: &[(i64, i64)],
    catalog: &HashMap<i64, &ProteinRow>,
) -> BTreeMap<i64, Vec<usize>> {
    let local: HashMap<i64, usize> = passing
        .iter()
        .enumerate()
        .map(|(index, peptide)| (peptide.id, index))
        .collect();
    let mut sets: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
    for (peptide, protein) in peptide_proteins {
        let (Some(&index), Some(row)) = (local.get(peptide), catalog.get(protein)) else {
            continue;
        };
        if row.is_decoy == passing[index].is_decoy {
            sets.entry(*protein).or_default().push(index);
        }
    }
    for set in sets.values_mut() {
        set.sort_unstable();
        set.dedup();
    }
    sets
}

fn entities(sets: BTreeMap<i64, Vec<usize>>, catalog: &HashMap<i64, &ProteinRow>) -> Vec<Entity> {
    let mut by_set: BTreeMap<Vec<usize>, Vec<i64>> = BTreeMap::new();
    for (protein, set) in sets {
        by_set.entry(set).or_default().push(protein);
    }
    by_set
        .into_iter()
        .map(|(peptides, mut proteins)| {
            proteins.sort_by(|a, b| {
                let key = |id: &i64| {
                    let accession = catalog[id].accession.as_str();
                    (!accession.starts_with("sp|"), accession.to_string())
                };
                key(a).cmp(&key(b))
            });
            Entity { peptides, proteins }
        })
        .collect()
}

/// Groups in razor order, with members, peptide roles, and scores;
/// q-values are filled later.
fn build_groups(
    entities: &[Entity],
    passing: &[&PeptideEvidence<'_>],
    catalog: &HashMap<i64, &ProteinRow>,
) -> Vec<ProteinGroup> {
    let mut containing: HashMap<usize, Vec<usize>> = HashMap::new();
    for (index, entity) in entities.iter().enumerate() {
        for &peptide in &entity.peptides {
            containing.entry(peptide).or_default().push(index);
        }
    }
    let candidates = razor_order(entities, &containing, passing, catalog);
    let mut owner: HashMap<usize, usize> = HashMap::new();
    for &candidate in &candidates {
        for &peptide in &entities[candidate].peptides {
            owner.entry(peptide).or_insert(candidate);
        }
    }
    // A candidate credited with nothing is explained by the others.
    let kept: Vec<usize> = candidates
        .into_iter()
        .filter(|&c| entities[c].peptides.iter().any(|p| owner[p] == c))
        .collect();
    let kept_set: HashSet<usize> = kept.iter().copied().collect();
    let groups = kept
        .iter()
        .map(|&group| {
            let sharers = |peptide: usize| {
                containing[&peptide]
                    .iter()
                    .filter(|index| kept_set.contains(index))
                    .count()
            };
            (
                group,
                make_group(&entities[group], group, &owner, sharers, passing),
            )
        })
        .collect();
    attach_subsets(entities, &owner, groups)
}

/// Entities that are not strict subsets of another, ordered for razor
/// credit: most peptides, then best peptide score, then leader accession.
fn razor_order(
    entities: &[Entity],
    containing: &HashMap<usize, Vec<usize>>,
    passing: &[&PeptideEvidence<'_>],
    catalog: &HashMap<i64, &ProteinRow>,
) -> Vec<usize> {
    let is_strict_subset = |index: usize| {
        let entity = &entities[index];
        containing[&entity.peptides[0]].iter().any(|&other| {
            other != index
                && entities[other].peptides.len() > entity.peptides.len()
                && contains_all(&entities[other].peptides, &entity.peptides)
        })
    };
    let best_score = |entity: &Entity| best_of(entity.peptides.iter().copied(), passing);
    let mut candidates: Vec<usize> = (0..entities.len())
        .filter(|&index| !is_strict_subset(index))
        .collect();
    candidates.sort_by(|&a, &b| {
        let (ea, eb) = (&entities[a], &entities[b]);
        eb.peptides
            .len()
            .cmp(&ea.peptides.len())
            .then(best_score(eb).total_cmp(&best_score(ea)))
            .then_with(|| {
                catalog[&ea.proteins[0]]
                    .accession
                    .cmp(&catalog[&eb.proteins[0]].accession)
            })
    });
    candidates
}

/// One group from its entity. `sharers` counts the kept groups containing
/// a peptide.
fn make_group(
    entity: &Entity,
    index: usize,
    owner: &HashMap<usize, usize>,
    sharers: impl Fn(usize) -> usize,
    passing: &[&PeptideEvidence<'_>],
) -> ProteinGroup {
    let peptides = entity
        .peptides
        .iter()
        .map(|&p| {
            let role = if sharers(p) == 1 {
                PeptideRole::Unique
            } else if owner[&p] == index {
                PeptideRole::Razor
            } else {
                PeptideRole::Shared
            };
            (passing[p].id, role)
        })
        .collect();
    let credited: Vec<usize> = entity
        .peptides
        .iter()
        .copied()
        .filter(|p| owner[p] == index)
        .collect();
    let sequences: HashSet<&str> = credited.iter().map(|&p| passing[p].sequence).collect();
    let members = entity
        .proteins
        .iter()
        .enumerate()
        .map(|(i, &protein)| {
            let role = if i == 0 {
                MemberRole::Leader
            } else {
                MemberRole::Equivalent
            };
            (protein, role)
        })
        .collect();
    ProteinGroup {
        leader: entity.proteins[0],
        members,
        peptides,
        is_decoy: passing[entity.peptides[0]].is_decoy,
        score: best_of(credited.iter().copied(), passing),
        q_value: None,
        pep: None,
        distinct_sequences: u32::try_from(sequences.len()).unwrap_or(u32::MAX),
    }
}

fn best_of(peptides: impl Iterator<Item = usize>, passing: &[&PeptideEvidence<'_>]) -> f64 {
    peptides
        .map(|p| passing[p].score)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Entities that are not groups join, as subset members, the groups that
/// contain all their peptides, or else the groups credited with them.
fn attach_subsets(
    entities: &[Entity],
    owner: &HashMap<usize, usize>,
    mut groups: Vec<(usize, ProteinGroup)>,
) -> Vec<ProteinGroup> {
    let position: HashMap<usize, usize> = groups
        .iter()
        .enumerate()
        .map(|(position, (entity, _))| (*entity, position))
        .collect();
    for (index, entity) in entities.iter().enumerate() {
        if position.contains_key(&index) {
            continue;
        }
        let mut owners: Vec<usize> = entity.peptides.iter().map(|p| owner[p]).collect();
        owners.sort_unstable();
        owners.dedup();
        let supersets: Vec<usize> = owners
            .iter()
            .copied()
            .filter(|o| contains_all(&entities[*o].peptides, &entity.peptides))
            .collect();
        let hosts = if supersets.is_empty() {
            owners
        } else {
            supersets
        };
        for host in hosts {
            let members = &mut groups[position[&host]].1.members;
            members.extend(entity.proteins.iter().map(|&p| (p, MemberRole::Subset)));
        }
    }
    groups.into_iter().map(|(_, group)| group).collect()
}

/// Whether sorted `outer` contains every element of sorted `inner`.
fn contains_all(outer: &[usize], inner: &[usize]) -> bool {
    let mut outer = outer.iter().peekable();
    inner.iter().all(|wanted| {
        while outer.next_if(|&&value| value < *wanted).is_some() {}
        outer.next_if_eq(&wanted).is_some()
    })
}

/// Fills q-values by picked competition; also returns a note when most
/// decoy groups found no target partner (their accessions carry no
/// recognized decoy label), so the estimate fell back to unpaired counting.
fn picked_q_values(
    mut groups: Vec<ProteinGroup>,
    catalog: &HashMap<i64, &ProteinRow>,
) -> (Vec<ProteinGroup>, Option<Note>) {
    let named = |group: &ProteinGroup| -> Vec<String> {
        group
            .members
            .iter()
            .filter(|(_, role)| *role != MemberRole::Subset)
            .map(|(protein, _)| catalog[protein].accession.clone())
            .collect()
    };
    let mut target_of: HashMap<String, usize> = HashMap::new();
    for (index, group) in groups.iter().enumerate() {
        if !group.is_decoy {
            for accession in named(group) {
                target_of.entry(accession).or_insert(index);
            }
        }
    }
    // Target group -> its best decoy competitor.
    let mut rival: HashMap<usize, usize> = HashMap::new();
    let mut partnered = HashSet::new();
    for (index, group) in groups.iter().enumerate().filter(|(_, g)| g.is_decoy) {
        let partner = named(group).iter().find_map(|accession| {
            target_accession(accession).and_then(|target| target_of.get(target).copied())
        });
        if let Some(partner) = partner {
            partnered.insert(index);
            let current = rival.entry(partner).or_insert(index);
            if groups[index].score > groups[*current].score {
                *current = index;
            }
        }
    }
    let winning_decoys: HashSet<usize> = rival
        .iter()
        .filter(|&(&target, &decoy)| groups[decoy].score >= groups[target].score)
        .map(|(_, &decoy)| decoy)
        .collect();
    let kept: Vec<usize> = (0..groups.len())
        .filter(|&index| {
            if groups[index].is_decoy {
                !partnered.contains(&index) || winning_decoys.contains(&index)
            } else {
                !rival
                    .get(&index)
                    .is_some_and(|decoy| winning_decoys.contains(decoy))
            }
        })
        .collect();
    let scores: Vec<f64> = kept.iter().map(|&i| groups[i].score).collect();
    let decoys: Vec<bool> = kept.iter().map(|&i| groups[i].is_decoy).collect();
    for (index, (q, pep)) in kept.iter().zip(confidence(&scores, &decoys)) {
        groups[*index].q_value = Some(q);
        groups[*index].pep = Some(pep);
    }
    let decoy_groups = groups.iter().filter(|group| group.is_decoy).count();
    // A labelled decoy whose target was never identified is unpaired by
    // nature; only decoys without a recognizable label are a problem.
    let labelled = groups
        .iter()
        .filter(|group| {
            group.is_decoy
                && named(group)
                    .iter()
                    .any(|accession| target_accession(accession).is_some())
        })
        .count();
    let note = (labelled * 2 < decoy_groups).then(|| {
        Note::new(
            NoteCode::UnpairedDecoys,
            format!(
                "Only {labelled} of {decoy_groups} decoy protein groups carry a label (such as \
                 rev_ or -DECOY) naming the target they were made from, so most decoys could \
                 not compete with a target and the protein FDR estimate is less reliable."
            ),
        )
    });
    (groups, note)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn protein(id: i64, accession: &str) -> ProteinRow {
        ProteinRow {
            id,
            accession: accession.to_string(),
            is_decoy: tuffold_core::decoy::is_decoy_accession(accession),
        }
    }

    fn peptide(id: i64, sequence: &'static str, score: f64) -> PeptideEvidence<'static> {
        PeptideEvidence {
            id,
            sequence,
            is_decoy: false,
            score,
            q_value: 0.001,
            pep: 0.0005,
        }
    }

    fn roles(group: &ProteinGroup) -> Vec<(i64, PeptideRole)> {
        let mut roles = group.peptides.clone();
        roles.sort_unstable();
        roles
    }

    #[test]
    fn contains_all_walks_sorted_sets() {
        assert!(contains_all(&[1, 2, 3, 5], &[2, 5]));
        assert!(!contains_all(&[1, 2, 3], &[2, 4]));
        assert!(contains_all(&[1], &[]));
    }

    #[test]
    fn parsimony_groups_equivalents_folds_subsets_and_credits_razor_peptides() {
        // Peptides 1..=5. Proteins:
        //   A {1,2,3}   B {1,2,3}  (indistinguishable from A; sp| leads)
        //   C {2,3}     (strict subset of A/B)
        //   D {3,4}     (shares 3 with A/B; owns 4)
        //   E {5}       (alone)
        //   F {4}       (subset of D)
        let proteins = [
            protein(1, "tr|A"),
            protein(2, "sp|B"),
            protein(3, "C"),
            protein(4, "D"),
            protein(5, "E"),
            protein(6, "F"),
        ];
        let map = [
            (1, 1),
            (2, 1),
            (3, 1),
            (1, 2),
            (2, 2),
            (3, 2),
            (2, 3),
            (3, 3),
            (3, 4),
            (4, 4),
            (5, 5),
            (4, 6),
        ];
        let peptides = [
            peptide(1, "AAAK", 10.0),
            peptide(2, "CCCK", 9.0),
            peptide(3, "DDDK", 8.0),
            peptide(4, "EEEK", 7.0),
            peptide(5, "FFFK", 6.0),
        ];
        let inference = infer_proteins(&peptides, &map, &proteins, &Thresholds::default());
        assert!(inference.notes.is_empty());
        let groups = &inference.groups;
        assert_eq!(groups.len(), 3, "{groups:#?}");
        // Razor order: {1,2,3} (3 peptides), then {3,4}, then {5}.
        assert_eq!(groups[0].leader, 2, "sp|B leads tr|A");
        assert_eq!(
            groups[0].members,
            [
                (2, MemberRole::Leader),
                (1, MemberRole::Equivalent),
                (3, MemberRole::Subset)
            ]
        );
        assert_eq!(
            roles(&groups[0]),
            [
                (1, PeptideRole::Unique),
                (2, PeptideRole::Unique),
                (3, PeptideRole::Razor)
            ]
        );
        assert_eq!((groups[0].score, groups[0].distinct_sequences), (10.0, 3));
        assert_eq!(groups[1].leader, 4);
        assert_eq!(
            groups[1].members,
            [(4, MemberRole::Leader), (6, MemberRole::Subset)]
        );
        assert_eq!(
            roles(&groups[1]),
            [(3, PeptideRole::Shared), (4, PeptideRole::Unique)]
        );
        // D is scored by peptide 4 only (7), not by the shared 3 (8).
        assert_eq!((groups[1].score, groups[1].distinct_sequences), (7.0, 1));
        assert_eq!(groups[2].leader, 5);
    }

    #[test]
    fn a_group_explained_by_others_becomes_a_subset() {
        // A {1,2}, B {2,3}, C {1,3}: no strict subsets. Razor order by best
        // score: A (peptide 1 = 10), then C (10 too; accession order A < C),
        // then B. A takes 1 and 2, C takes 3; B is left with nothing.
        let proteins = [protein(1, "A"), protein(2, "B"), protein(3, "C")];
        let map = [(1, 1), (2, 1), (2, 2), (3, 2), (1, 3), (3, 3)];
        let peptides = [
            peptide(1, "AK", 10.0),
            peptide(2, "CK", 5.0),
            peptide(3, "DK", 4.0),
        ];
        let groups = infer_proteins(&peptides, &map, &proteins, &Thresholds::default()).groups;
        let leaders: Vec<i64> = groups.iter().map(|g| g.leader).collect();
        assert_eq!(leaders, [1, 3]);
        // B's peptides are credited to A (2) and C (3); neither has both, so
        // it joins both as a subset.
        assert!(groups[0].members.contains(&(2, MemberRole::Subset)));
        assert!(groups[1].members.contains(&(2, MemberRole::Subset)));
        // With B gone, peptide 2 belongs to A alone and 3 to C alone; peptide
        // 1 is shared by A and C and credited to A.
        assert_eq!(
            roles(&groups[0]),
            [(1, PeptideRole::Razor), (2, PeptideRole::Unique)]
        );
        assert_eq!(
            roles(&groups[1]),
            [(1, PeptideRole::Shared), (3, PeptideRole::Unique)]
        );
    }

    #[test]
    fn peptides_below_the_threshold_are_not_evidence() {
        let proteins = [protein(1, "A")];
        let mut weak = peptide(1, "AK", 1.0);
        weak.q_value = 0.2;
        let inference = infer_proteins(
            std::slice::from_ref(&weak),
            &[(1, 1)],
            &proteins,
            &Thresholds::default(),
        );
        assert!(inference.groups.is_empty());
        assert_eq!(inference.notes[0].code, NoteCode::NoPassingPeptides);
        assert!(inference.notes[0].message.contains("1% FDR"));
        // By probability the same peptide (PEP 0.0005, 99.95%) is evidence.
        let by_probability = Thresholds::probability_default();
        let inference = infer_proteins(&[weak], &[(1, 1)], &proteins, &by_probability);
        assert_eq!(inference.groups.len(), 1);
    }

    #[test]
    fn decoys_without_a_recognized_label_are_reported_unpaired() {
        // Decoys flagged only by mzIdentML isDecoy, accessions unlabelled.
        let proteins = [
            protein(1, "X"),
            ProteinRow {
                id: 2,
                accession: "X_shuffled".into(),
                is_decoy: true,
            },
        ];
        let decoy = PeptideEvidence {
            is_decoy: true,
            ..peptide(2, "KX", 1.0)
        };
        let inference = infer_proteins(
            &[peptide(1, "XK", 9.0), decoy],
            &[(1, 1), (2, 2)],
            &proteins,
            &Thresholds::default(),
        );
        assert_eq!(inference.notes.len(), 1);
        assert_eq!(inference.notes[0].code, NoteCode::UnpairedDecoys);
        assert!(inference.notes[0].message.starts_with("Only 0 of 1 decoy"));
    }

    #[test]
    fn decoys_compete_with_their_own_target_only() {
        // Target groups: X (score 9), Y (score 3), Z (score 8).
        // Decoy groups: rev_Y (score 5) beats Y; rev_X (score 2) loses to X;
        // rev_W (score 4) has no target partner and stays.
        let proteins = [
            protein(1, "X"),
            protein(2, "Y"),
            protein(3, "Z"),
            protein(4, "rev_Y"),
            protein(5, "rev_X"),
            protein(6, "rev_W"),
        ];
        let map = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6)];
        let decoy = |id, sequence, score| PeptideEvidence {
            is_decoy: true,
            ..peptide(id, sequence, score)
        };
        let peptides = [
            peptide(1, "XK", 9.0),
            peptide(2, "YK", 3.0),
            peptide(3, "ZK", 8.0),
            decoy(4, "KY", 5.0),
            decoy(5, "KX", 2.0),
            decoy(6, "KW", 4.0),
        ];
        let inference = infer_proteins(&peptides, &map, &proteins, &Thresholds::default());
        assert!(inference.notes.is_empty(), "2 of 3 decoys paired");
        let groups = inference.groups;
        let q = |leader: i64| groups.iter().find(|g| g.leader == leader).unwrap().q_value;
        // Kept, best first: X 9 (T), Z 8 (T), rev_Y 5 (D), rev_W 4 (D).
        // FDR: 1/1, 1/2, 2/2, 3/2 -> monotone from the bottom: 1/2, 1/2, 1, 1.
        assert_eq!(q(1), Some(0.5));
        assert_eq!(q(3), Some(0.5));
        assert_eq!(q(4), Some(1.0));
        assert_eq!(q(6), Some(1.0));
        assert_eq!(q(2), None, "Y lost to rev_Y");
        assert_eq!(q(5), None, "rev_X lost to X");
        // PEPs over the same kept groups, best first T T D D: the isotonic
        // decoy fractions are 0, 0, 1, 1, floored at the q-values 1/2, 1/2,
        // 1, 1, so PEPs 1/2, 1/2, 1, 1. Losers: unknown.
        let pep = |leader: i64| groups.iter().find(|g| g.leader == leader).unwrap().pep;
        let peps: Vec<Option<f64>> = [1, 3, 4, 6, 2, 5].into_iter().map(pep).collect();
        assert_eq!(
            peps,
            [Some(0.5), Some(0.5), Some(1.0), Some(1.0), None, None]
        );
    }
}
