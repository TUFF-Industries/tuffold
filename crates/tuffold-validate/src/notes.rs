//! Findings about a validation run that the user should know.

/// Stable vocabulary: front ends and tests branch on these codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NoteCode {
    /// A level (or the whole project) has no decoys, so no FDR is estimated.
    NoDecoys,
    /// An engine's PSMs could not be ranked (no decoys or no shared score).
    EngineUnscored,
    /// Rescoring was asked for but an engine fell back to its best score.
    RescoringFallback,
    /// Peptides passing the threshold with no protein to group under.
    PeptidesWithoutProteins,
    /// No peptide passed the peptide threshold, so no protein groups exist.
    NoPassingPeptides,
    /// Most decoy protein groups had no target partner for picked FDR.
    UnpairedDecoys,
    /// An engine's rank-1 PSMs hold so few decoys that its results look
    /// already filtered, so Tuffold's FDR checks little.
    FewDecoys,
    /// Sibling peptides were not used: decoys could not be paired with their
    /// targets, so they could not count the errors that evidence adds
    /// (validations before decoys' counts were imputed).
    SiblingsUnpaired,
    /// Decoys could not be paired with their targets by name, so their
    /// sibling counts were imputed from low-scoring targets'.
    SiblingsImputed,
    /// What the deep-learning models did, or why they weren't used.
    DeepLearning,
    /// Replicate spectra were not used: decoys could not be paired with
    /// their target peptides (validations before decoys' counts were
    /// imputed).
    ReplicatesUnpaired,
    /// Decoys could not be paired with their target peptides, so their
    /// replicate counts were imputed from low-scoring targets'.
    ReplicatesImputed,
    /// Which fragment tolerance PTM site localization used, and the matches
    /// that decided it.
    Localization,
}

impl NoteCode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NoDecoys => "no-decoys",
            Self::EngineUnscored => "engine-unscored",
            Self::RescoringFallback => "rescoring-fallback",
            Self::PeptidesWithoutProteins => "peptides-without-proteins",
            Self::NoPassingPeptides => "no-passing-peptides",
            Self::UnpairedDecoys => "unpaired-decoys",
            Self::FewDecoys => "few-decoys",
            Self::SiblingsUnpaired => "siblings-unpaired",
            Self::SiblingsImputed => "siblings-imputed",
            Self::DeepLearning => "deep-learning",
            Self::ReplicatesUnpaired => "replicates-unpaired",
            Self::ReplicatesImputed => "replicates-imputed",
            Self::Localization => "localization",
        }
    }

    pub fn parse(code: &str) -> Option<Self> {
        [
            Self::NoDecoys,
            Self::EngineUnscored,
            Self::RescoringFallback,
            Self::PeptidesWithoutProteins,
            Self::NoPassingPeptides,
            Self::UnpairedDecoys,
            Self::FewDecoys,
            Self::SiblingsUnpaired,
            Self::SiblingsImputed,
            Self::DeepLearning,
            Self::ReplicatesUnpaired,
            Self::ReplicatesImputed,
            Self::Localization,
        ]
        .into_iter()
        .find(|known| known.as_str() == code)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Note {
    pub code: NoteCode,
    pub message: String,
}

impl Note {
    pub fn new(code: NoteCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_round_trip() {
        for code in [
            NoteCode::NoDecoys,
            NoteCode::EngineUnscored,
            NoteCode::RescoringFallback,
            NoteCode::PeptidesWithoutProteins,
            NoteCode::NoPassingPeptides,
            NoteCode::UnpairedDecoys,
            NoteCode::FewDecoys,
            NoteCode::SiblingsUnpaired,
            NoteCode::SiblingsImputed,
            NoteCode::DeepLearning,
            NoteCode::ReplicatesUnpaired,
            NoteCode::ReplicatesImputed,
            NoteCode::Localization,
        ] {
            assert_eq!(NoteCode::parse(code.as_str()), Some(code));
        }
        assert_eq!(NoteCode::parse("unknown"), None);
    }
}
