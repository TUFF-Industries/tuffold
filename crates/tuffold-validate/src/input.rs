//! What validation reads from a project.

/// The project's top-ranked PSMs with every score, and the peptide-protein
/// map. Rows refer to one another by store ids.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Input {
    /// Search engine names; [`PsmRow::engine`] indexes this.
    pub engines: Vec<String>,
    /// Score names; the first element of each [`PsmRow::scores`] pair
    /// indexes this.
    pub score_names: Vec<String>,
    /// Rank-1 PSMs only: lower ranks never compete for a spectrum.
    pub psms: Vec<PsmRow>,
    pub peptides: Vec<PeptideRow>,
    pub proteins: Vec<ProteinRow>,
    /// Distinct `(peptide id, protein id)` evidence pairs.
    pub peptide_proteins: Vec<(i64, i64)>,
    /// Features computed from each PSM's spectrum before validation (e.g.
    /// how well predicted fragment intensities match the measured ones).
    pub predicted: Option<PredictedFeatures>,
}

/// Per-PSM features from outside validation, higher always better. A PSM
/// without them (no peaks, a peptide the predictor can't take) gets each
/// column's lowest value among the engine's PSMs.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PredictedFeatures {
    pub names: Vec<String>,
    /// PSM id → one value per name.
    pub values: std::collections::HashMap<i64, Vec<f64>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PsmRow {
    pub id: i64,
    pub spectrum_id: i64,
    pub peptide_id: i64,
    pub engine: usize,
    pub charge: i32,
    pub calculated_mz: Option<f64>,
    pub experimental_mz: Option<f64>,
    /// The MS run the spectrum belongs to (a store id).
    pub run_id: i64,
    /// The spectrum's retention time in seconds, when known.
    pub retention_time: Option<f64>,
    /// `(score name index, value)`.
    pub scores: Vec<(usize, f64)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PeptideRow {
    pub id: i64,
    /// Unmodified residues.
    pub sequence: String,
    pub modifications: u32,
    /// A peptide is a decoy when every protein it maps to is a decoy.
    pub is_decoy: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinRow {
    pub id: i64,
    pub accession: String,
    pub is_decoy: bool,
}
