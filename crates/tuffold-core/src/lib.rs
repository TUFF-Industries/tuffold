//! Tuffold data model (v2).
//!
//! - [`batch`]: what importers produce — runs, spectra, peptide forms, PSMs
//!   with every engine score, proteins, and evidence.
//! - [`peptide`]: modified peptide forms and their ProForma-style keys.
//! - [`chemistry`]: residue and terminal masses.
//! - [`decoy`]: decoy recognition by accession convention.
//! - [`fragments`]: b/y fragment ions and peak annotation.
//! - [`coverage`]: where peptides fall on a protein sequence.
//! - [`digest`]: tryptic peptides a protein could show (for emPAI).
//! - [`reporter`]: isobaric labels (TMT, TMTpro) and their reporter ions.

pub mod batch;
pub mod centroid;
pub mod chemistry;
pub mod coverage;
pub mod decoy;
pub mod digest;
pub mod fragments;
pub mod peptide;
pub mod reporter;
pub mod unimod;
pub mod uniprot;

pub use batch::{
    BioSample, DecoyIndex, Evidence, ImportBatch, MsRun, PeakSpectrum, Protein, Psm, Score,
    Software, SourceFormat, Spectrum, SpectrumKey,
};
pub use coverage::SequenceCoverage;
pub use peptide::{InvalidModification, ModSite, Modification, PeptideForm};
