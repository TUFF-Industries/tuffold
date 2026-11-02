//! Tuffold validation: from engine PSMs to FDR-controlled spectra,
//! peptides, and protein groups.
//!
//! Pure functions over plain rows ([`Input`]); the store loads and saves them
//! and the engine crate runs the workflow. Every level uses target-decoy
//! competition with the +1 correction ([`fdr`]), so the decoys a search
//! reported are the only model of error — nothing is assumed about score
//! distributions.
//!
//! - [`validate`] ranks each engine's PSMs ([`rescore`]), combines engines
//!   per spectrum, scores peptides, and infers proteins.
//! - [`infer_proteins`] reruns protein inference alone, for a new peptide
//!   threshold.

pub mod fdr;
mod features;
pub mod inference;
mod input;
mod lda;
pub mod localize;
mod logistic;
mod notes;
mod pooling;
mod replicates;
pub mod rescore;
mod retention;
mod run;
mod settings;
mod siblings;
pub mod transfer;

pub use inference::{
    MemberRole, PeptideEvidence, PeptideRole, ProteinGroup, ProteinInference, infer_proteins,
};
pub use input::{Input, PeptideRow, PredictedFeatures, ProteinRow, PsmRow};
pub use localize::{
    DEFAULT_SITE_PROBABILITY, LocalizationInput, LocalizationSummary, MAX_ISOFORMS,
    ScoredLocalization, SiteScore, SpectrumLocalization, ambiguous_isoform_count, estimate_flr,
    format_sites, localize_spectrum, summarize_localizations,
};
pub use notes::{Note, NoteCode};
pub use pooling::{SIBLINGS, SIBLINGS_IMPUTED};
pub use replicates::{REPLICATES, REPLICATES_IMPUTED};
pub use rescore::{BestScoreReason, Method, RescoreModel, UnscoredReason};
pub use retention::{DEEP_RT_ERROR, RT_ERROR, RetentionPredictor, TrainedRetention};
pub use run::{
    EngineScore, EngineSummary, PeptideResult, PsmResult, SpectrumMatch, Validation,
    ValidationError, engine_scores, validate, validate_with,
};
pub use settings::{Level, Settings, ThresholdMode, ThresholdRequest, Thresholds, format_rate};
