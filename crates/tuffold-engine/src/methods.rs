//! A methods paragraph for publication, written from what the project
//! actually did: the engines and versions imported, how each was ranked,
//! the thresholds, and what passed. Every number in it comes from the
//! project, and nothing is claimed that wasn't done.

use tuffold_predict::FEATURE_NAMES;
use tuffold_store::{Project, StoreError, ValidationSummary};
use tuffold_validate::{
    DEEP_RT_ERROR, Level, REPLICATES, REPLICATES_IMPUTED, RT_ERROR, SIBLINGS, SIBLINGS_IMPUTED,
    ThresholdMode,
};

#[derive(Debug, thiserror::Error)]
pub enum MethodsError {
    #[error(transparent)]
    Store(#[from] StoreError),
    #[error("this project has no current validation results; validate it first")]
    NotValidated,
}

/// Paragraphs of methods text, then the references they cite.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Methods {
    pub paragraphs: Vec<String>,
    pub references: Vec<&'static str>,
}

impl Methods {
    pub fn text(&self) -> String {
        let mut text = self.paragraphs.join("\n\n");
        text.push_str("\n\nReferences\n");
        for (n, reference) in self.references.iter().enumerate() {
            text.push_str(&format!("{}. {reference}\n", n + 1));
        }
        text
    }
}

const ELIAS: &str = "Elias JE, Gygi SP. Target-decoy search strategy for increased confidence \
                     in large-scale protein identifications by mass spectrometry. Nat Methods. \
                     2007;4(3):207-214.";
const PERCOLATOR: &str = "Käll L, Canterbury JD, Weston J, Noble WS, MacCoss MJ. \
                          Semi-supervised learning for peptide identification from shotgun \
                          proteomics datasets. Nat Methods. 2007;4(11):923-925.";
const IPROPHET: &str = "Shteynberg D, Deutsch EW, Lam H, Eng JK, Sun Z, Tasman N, Mendoza L, \
                        Moritz RL, Aebersold R, Nesvizhskii AI. iProphet: multi-level \
                        integrative analysis of shotgun proteomic data improves peptide and \
                        protein identification rates and error estimates. Mol Cell Proteomics. \
                        2011;10(12):M111.007690.";
const ALPHAPEPTDEEP: &str = "Zeng WF, Zhou XX, Willems S, Ammar C, Wahle M, Bludau I, Voytik E, \
                             Strauss MT, Mann M. AlphaPeptDeep: a modular deep learning \
                             framework to predict peptide properties for proteomics. Nat Commun. \
                             2022;13(1):7238.";
const PICKED: &str = "Savitski MM, Wilhelm M, Hahne H, Kuster B, Bantscheff M. A scalable \
                      approach for protein false discovery rate estimation in large proteomic \
                      data sets. Mol Cell Proteomics. 2015;14(9):2394-2404.";
const SAGE: &str = "Lazear MR. Sage: an open-source tool for fast proteomics searching and \
                    quantification at scale. J Proteome Res. 2023;22(11):3652-3659.";
const MAXLFQ: &str = "Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate \
                      proteome-wide label-free quantification by delayed normalization and \
                      maximal peptide ratio extraction, termed MaxLFQ. Mol Cell Proteomics. \
                      2014;13(9):2513-2526.";
const TUKEY: &str = "Tukey JW. Exploratory Data Analysis. Reading, MA: Addison-Wesley; 1977.";
const NSAF: &str = "Zybailov B, Mosley AL, Sardiu ME, Coleman MK, Florens L, Washburn MP. \
                    Statistical analysis of membrane proteome expression changes in \
                    Saccharomyces cerevisiae. J Proteome Res. 2006;5(9):2339-2347.";
const BH: &str = "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical \
                  and powerful approach to multiple testing. J R Stat Soc Series B. \
                  1995;57(1):289-300.";
const ASCORE: &str = "Beausoleil SA, Villén J, Gerber SA, Rush J, Gygi SP. A probability-based \
                      approach for high-throughput protein phosphorylation analysis and site \
                      localization. Nat Biotechnol. 2006;24(10):1285-1292.";
const PTMRS: &str = "Taus T, Köcher T, Pichler P, Paschke C, Schmidt A, Henrich C, Mechtler K. \
                     Universal and confident phosphorylation site localization using phosphoRS. \
                     J Proteome Res. 2011;10(12):5354-5362.";
const DECOY_FLR: &str = "Ramsbottom KA, Prakash A, Riverol Y, Camacho OM, Martin MJ, Vizcaíno JA, \
                         Deutsch EW, Jones AR. Method for independent estimation of the false \
                         localization rate for phosphoproteomics. J Proteome Res. \
                         2022;21(7):1603-1615.";

/// Cite `reference`, returning its number in the list.
fn cite(references: &mut Vec<&'static str>, reference: &'static str) -> usize {
    match references.iter().position(|r| *r == reference) {
        Some(i) => i + 1,
        None => {
            references.push(reference);
            references.len()
        }
    }
}

fn list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [init @ .., last] => format!("{} and {last}", init.join(", ")),
    }
}

/// The methods of a validated project.
pub fn methods(project: &Project) -> Result<Methods, MethodsError> {
    let summary = project
        .validation_summary()?
        .filter(|s| s.validated_at.is_some())
        .ok_or(MethodsError::NotValidated)?;
    let mut references = Vec::new();
    let mut paragraphs = vec![identification(project, &mut references)?];
    paragraphs.push(validation(&summary, &mut references));
    paragraphs.push(results(&summary));
    paragraphs.push(quantitation(project, &mut references)?);
    Ok(Methods {
        paragraphs,
        references,
    })
}

fn identification(
    project: &Project,
    references: &mut Vec<&'static str>,
) -> Result<String, StoreError> {
    let mut software: Vec<String> = Vec::new();
    let mut searched = false;
    let mut synthetic = false;
    for import in project.imports()? {
        searched |= import.software == "Sage";
        synthetic |= import.format == "synthetic";
        let named = match &import.software_version {
            Some(version) => format!("{} {version}", import.software),
            None => import.software.clone(),
        };
        if !software.contains(&named) {
            software.push(named);
        }
    }
    let project_summary = project.summary()?;
    let mut text = String::new();
    if synthetic {
        text.push_str(
            "THESE DATA ARE SYNTHETIC: Tuffold's generated demo project; nothing was measured. ",
        );
    }
    text.push_str(&format!(
        "Peptide-spectrum matches (PSMs) for {} spectra from {} MS runs ({} biosamples) were \
         validated with Tuffold {}, an open-source (Apache-2.0) proteomics validation tool.",
        project_summary.spectra,
        project_summary.ms_runs,
        project_summary.biosamples,
        env!("CARGO_PKG_VERSION")
    ));
    if searched {
        let n = cite(references, SAGE);
        text.push_str(&format!(
            " Spectra were searched with Tuffold's built-in search, which uses the search core \
             of Sage [{n}] and generates reversed decoy sequences."
        ));
    }
    text.push_str(&format!(" Search results came from {}.", list(&software)));
    Ok(text)
}

fn validation(summary: &ValidationSummary, references: &mut Vec<&'static str>) -> String {
    let decoys = cite(references, ELIAS);
    let mut text = format!(
        "False discovery rates were estimated by target-decoy competition [{decoys}] at the \
         spectrum, peptide, and protein-group level, with the +1 correction, and \
         posterior error probabilities by isotonic regression of the decoy fraction."
    );
    let rescored: Vec<String> = summary
        .engines
        .iter()
        .filter(|e| e.method == "rescored")
        .map(|e| e.name.clone())
        .collect();
    if summary.settings.rescore && !rescored.is_empty() {
        let n = cite(references, PERCOLATOR);
        text.push_str(&format!(
            " Rank-1 PSMs of {} were rescored semi-supervised against decoys, \
             cross-validated over three folds [{n}], by logistic regression or a linear \
             discriminant, whichever accepted more training PSMs.",
            list(&rescored)
        ));
        let uses = |feature: &str| {
            summary
                .engines
                .iter()
                .any(|e| e.method == "rescored" && e.features.iter().any(|f| f == feature))
        };
        let timed = uses(RT_ERROR);
        let (spectra, deep_time) = (uses(FEATURE_NAMES[0]), uses(DEEP_RT_ERROR));
        if spectra || deep_time {
            let n = cite(references, ALPHAPEPTDEEP);
            let what = match (spectra, deep_time) {
                (true, true) => {
                    "the agreement of each PSM's spectrum with the fragment intensities \
                     AlphaPeptDeep's pretrained model predicts for its peptide (spectral angle and \
                     correlations; collision energy and instrument calibrated on a sample of \
                     rank-1 PSMs), and its retention-time error against AlphaPeptDeep's \
                     retention-time model fine-tuned per run on confident target peptides, \
                     cross-fitted on the same folds"
                }
                (true, false) => {
                    "the agreement of each PSM's spectrum with the fragment intensities \
                     AlphaPeptDeep's pretrained model predicts for its peptide (spectral angle and \
                     correlations; collision energy and instrument calibrated on a sample of \
                     rank-1 PSMs)"
                }
                _ => {
                    "each PSM's retention-time error against AlphaPeptDeep's retention-time \
                     model fine-tuned per run on confident target peptides, cross-fitted on the \
                     same folds"
                }
            };
            text.push_str(&format!(" Features included {what} [{n}]."));
        }
        if timed {
            text.push_str(
                " Features included each PSM's retention-time error against a per-run model \
                 of retention time from amino-acid composition (ridge regression with a \
                 running-median calibration, trained on confident target peptides and \
                 cross-fitted on the same folds).",
            );
        }
        if uses(REPLICATES) {
            let n = cite(references, IPROPHET);
            text.push_str(&format!(
                " A second rescoring added the number of other spectra of each PSM's peptide \
                 scoring above the first rescoring's 0.5% cut [{n}]; a decoy PSM took the larger \
                 of its own count and that of the target peptide its sequence was made from, so \
                 that decoys model the errors this evidence adds."
            ));
        }
        if uses(REPLICATES_IMPUTED) {
            let n = cite(references, IPROPHET);
            text.push_str(&format!(
                " A second rescoring added the number of other spectra of each PSM's peptide \
                 scoring above the first rescoring's 0.5% cut [{n}]; decoy sequences named no \
                 target peptide, so each decoy peptide took a count drawn from those of \
                 low-scoring target PSMs (its PSMs keeping their own count when larger), so that \
                 decoys model the errors this evidence adds."
            ));
        }
    }
    if let Some(scoring) = summary
        .peptide_scoring
        .as_ref()
        .filter(|scoring| scoring.method == "pooled")
    {
        let n = cite(references, IPROPHET);
        let siblings = if scoring.features.iter().any(|feature| feature == SIBLINGS) {
            ", and the number of other identified peptides of their proteins (a decoy \
             peptide counting those of its paired target protein, so that decoys model the \
             errors this evidence adds),"
        } else if scoring
            .features
            .iter()
            .any(|feature| feature == SIBLINGS_IMPUTED)
        {
            ", and the number of other identified peptides of their proteins (decoy \
             proteins named no target, so each decoy peptide took a count drawn from those \
             of low-scoring targets, so that decoys model the errors this evidence adds),"
        } else {
            ""
        };
        text.push_str(&format!(
            " Peptides were scored from all their spectra, as in iProphet [{n}]: their best \
             and second-best spectrum scores and the numbers of other spectra and charges \
             matching them{siblings} were combined by the same cross-validated, \
             semi-supervised rescoring, against decoy peptides; modified forms were not \
             pooled."
        ));
    }
    let unranked: Vec<String> = summary
        .engines
        .iter()
        .filter(|e| e.method == "unscored")
        .map(|e| e.name.clone())
        .collect();
    if !unranked.is_empty() {
        text.push_str(&format!(
            " {} could not be ranked (no decoys, or no complete score) and did not contribute.",
            list(&unranked)
        ));
    }
    if summary.engines.len() > 1 {
        text.push_str(
            " Where several engines matched a spectrum, each engine's calibrated vote was \
             compared and the strongest decided the peptide.",
        );
    }
    let picked = cite(references, PICKED);
    let thresholds = &summary.settings.thresholds;
    let mode = match thresholds.mode {
        ThresholdMode::Fdr => "",
        ThresholdMode::Probability => " (probability = 1 − posterior error probability)",
    };
    text.push_str(&format!(
        " Proteins were grouped by parsimony (indistinguishable proteins together; shared \
         peptides assigned as razor peptides to the group with the most peptides), and \
         protein-group FDR was estimated by picked target-decoy competition [{picked}]. \
         Thresholds were {} for spectra, {} for peptides, and {} for protein groups{mode}, \
         with at least {} peptide{} per protein group.",
        thresholds.describe(Level::Psm),
        thresholds.describe(Level::Peptide),
        thresholds.describe(Level::Protein),
        thresholds.min_peptides,
        if thresholds.min_peptides == 1 {
            ""
        } else {
            "s"
        }
    ));
    if summary.localization.scored_psms > 0 {
        let ptmrs = cite(references, PTMRS);
        let ascore = cite(references, ASCORE);
        let rate = summary
            .localization
            .estimated_flr
            .map_or_else(String::new, |r| {
                let flr = cite(references, DECOY_FLR);
                format!(
                    " (estimated false localization rate {} from decoy amino acids [{flr}])",
                    tuffold_validate::format_rate(r)
                )
            });
        text.push_str(&format!(
            " Variable modifications with more than one candidate residue on a peptide were \
             localized by scoring each candidate arrangement against its site-determining b and \
             y fragment ions (20 ppm, or 0.5 Da for low-resolution spectra; up to 256 \
             arrangements per spectrum): pairwise binomial scores over peak depths 1 to 10 in \
             100 m/z windows, weighted as in the Ascore algorithm [{ascore}], gave arrangement \
             and site probabilities in the manner of phosphoRS [{ptmrs}], and each site's Ascore \
             is reported. {} of {} ambiguous passing PSMs ({} of {} sites) were localized at >= \
             {:.0}% site probability{rate}.",
            summary.localization.localized_psms,
            summary.localization.scored_psms,
            summary.localization.localized_sites,
            summary.localization.total_target_sites,
            summary.localization.probability_threshold * 100.0,
        ));
    }
    text
}

fn results(summary: &ValidationSummary) -> String {
    let level = |name: &str, counts: &tuffold_store::LevelCounts| {
        let fdr = counts.estimated_fdr().map_or_else(String::new, |f| {
            format!(", decoy-estimated FDR {}", tuffold_validate::format_rate(f))
        });
        format!("{} {name} ({} decoys{fdr})", counts.targets, counts.decoys)
    };
    format!(
        "At these thresholds, {}, {}, and {} passed.",
        level("spectra", &summary.spectra),
        level("peptides", &summary.peptides),
        level("protein groups", &summary.protein_groups)
    )
}

fn quantitation(
    project: &Project,
    references: &mut Vec<&'static str>,
) -> Result<String, StoreError> {
    let nsaf = cite(references, NSAF);
    let mut text = format!(
        "Spectral counts were taken from the passing spectra of each protein group's unique \
         and razor peptides, normalized to the mean total per biosample, and reported with \
         NSAF [{nsaf}]."
    );
    if project.intensity_list()?.has_ms1() {
        let maxlfq = cite(references, MAXLFQ);
        let tukey = cite(references, TUKEY);
        text.push_str(&format!(
            " MS1 precursor intensities were measured by tracing each identified precursor's \
             monoisotopic and M+1 peaks (10 ppm) through the MS1 scans and integrating its \
             elution peak. Peptide intensities were normalized by the half-sample mode of \
             pairwise peptide ratios between biosamples, fitted by least squares as in \
             MaxLFQ [{maxlfq}]. Protein intensities were summarized from the precursor areas \
             by Tukey's median polish [{tukey}] over the biosamples connected by shared \
             peptides, without imputing missing values."
        ));
    }
    let bh = cite(references, BH);
    text.push_str(&format!(
        " Two categories were compared with Welch's t-test on log2 values (for spectral \
         counts without replicates, Fisher's exact test on pooled counts), and several at once \
         with one-way ANOVA, with Benjamini-Hochberg q-values [{bh}]."
    ));
    Ok(text)
}
