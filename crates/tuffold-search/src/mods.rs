//! Sage's peptides as Tuffold peptide forms. Sage stores one f32 mass per
//! residue and one per terminus, and never stacks modifications in one slot
//! (fixed ones skip occupied residues, variable ones fill only empty slots),
//! so each mass is exactly one configured rule. Matching it back writes the
//! rule's exact mass and Unimod name; a mass that matches no rule, or
//! several different ones, keeps Sage's value unnamed rather than a guess.

use sage_core::peptide::Peptide;
use tuffold_core::{ModSite, Modification, PeptideForm};

use crate::settings::{ModRule, SearchSettings, Site};

/// f32 masses of a few hundred Da agree with f64 to ~3e-5. Narrower than
/// TMTpro and iTRAQ8plex's 1.8e-3 difference, so those stay apart.
const MASS_MATCH: f64 = 5e-4;

pub(crate) fn peptide_form(peptide: &Peptide, settings: &SearchSettings) -> PeptideForm {
    let sequence = String::from_utf8_lossy(&peptide.sequence).into_owned();
    let rules: Vec<&ModRule> = settings.fixed.iter().chain(&settings.variable).collect();
    let at = |accept: &dyn Fn(Site) -> bool| -> Vec<&ModRule> {
        rules
            .iter()
            .copied()
            .filter(|rule| accept(rule.site))
            .collect()
    };
    let mut modifications: Vec<Modification> = Vec::new();
    if let Some(delta) = peptide.nterm {
        modifications.push(matched(ModSite::NTerm, delta, &at(&Site::is_n_terminal)));
    }
    for (ix, (residue, &delta)) in sequence.chars().zip(&peptide.modifications).enumerate() {
        if delta != 0.0 {
            let position = u32::try_from(ix + 1).unwrap_or(u32::MAX);
            let here = at(&|site| site == Site::Residue(residue));
            modifications.push(matched(ModSite::Residue(position), delta, &here));
        }
    }
    if let Some(delta) = peptide.cterm {
        modifications.push(matched(ModSite::CTerm, delta, &at(&Site::is_c_terminal)));
    }
    // Every site above lies on the sequence, so this cannot fail; if it ever
    // did, keep the residues and say nothing rather than invent positions.
    PeptideForm::new(&sequence, modifications)
        .unwrap_or_else(|_| PeptideForm::unmodified(&sequence))
}

/// The one configured rule `delta` is; unnamed when none match, or several
/// that are different modifications. Rules for one modification at
/// different sites of a kind (`[:Acetyl` and `^:Acetyl` both match an
/// N-terminal +42.0106) agree on everything written, so either names it.
fn matched(site: ModSite, delta: f32, candidates: &[&ModRule]) -> Modification {
    let delta = f64::from(delta);
    let hits: Vec<&ModRule> = candidates
        .iter()
        .copied()
        .filter(|rule| (rule.mass - delta).abs() <= MASS_MATCH)
        .collect();
    let same = |a: &ModRule, b: &ModRule| {
        a.name.is_some() && a.name == b.name && a.accession == b.accession && a.mass == b.mass
    };
    match hits.split_first() {
        Some((rule, rest)) if rest.iter().all(|other| same(rule, other)) => Modification {
            site,
            mass_delta: rule.mass,
            name: rule.name.clone(),
            accession: rule.accession.clone(),
        },
        _ => Modification {
            site,
            mass_delta: delta,
            name: None,
            accession: None,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rule(text: &str) -> ModRule {
        ModRule::parse(text).unwrap()
    }

    #[test]
    fn a_mass_is_the_one_rule_it_matches() {
        let tmt = rule("K:TMT6plex");
        let other = rule("K:42.010565");
        let candidates = [&tmt, &other];
        let found = matched(ModSite::Residue(3), 229.162_932_f64 as f32, &candidates);
        assert_eq!(
            (found.mass_delta, found.name.as_deref()),
            (229.162_932, Some("TMT6plex"))
        );
        // Matches nothing: Sage's mass, unnamed.
        let odd = matched(ModSite::NTerm, 1.5, &candidates);
        assert_eq!((odd.name, odd.mass_delta), (None, 1.5));
    }

    #[test]
    fn close_rules_stay_apart_and_ties_stay_unnamed() {
        let tmtpro = rule("K:TMTpro");
        let itraq = rule("K:iTRAQ8plex");
        let candidates = [&tmtpro, &itraq];
        let found = matched(ModSite::Residue(1), 304.207_146_f64 as f32, &candidates);
        assert_eq!(found.name.as_deref(), Some("TMTpro"));
        // Two different rules Sage's f32 can't tell apart (TMTpro, and an
        // unnamed mass 5e-5 from it): which one is unknowable, so no name.
        let twin = rule("K:304.2072");
        let tied = matched(
            ModSite::Residue(1),
            304.207_146_f64 as f32,
            &[&tmtpro, &twin],
        );
        assert_eq!(tied.name, None);
    }

    #[test]
    fn one_modification_at_two_terminal_rules_is_named() {
        let protein = rule("[:Acetyl");
        let peptide = rule("^:Acetyl");
        let found = matched(ModSite::NTerm, 42.010_565_f64 as f32, &[&protein, &peptide]);
        assert_eq!(
            (found.name.as_deref(), found.mass_delta),
            (Some("Acetyl"), 42.010_565)
        );
    }
}
