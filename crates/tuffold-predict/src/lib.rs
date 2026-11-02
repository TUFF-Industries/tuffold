//! Deep-learned predictions for rescoring: fragment intensities and
//! retention times from AlphaPeptDeep's pretrained models (Zeng et al.,
//! Nat Commun 2022; Apache-2.0), run in pure Rust on candle.
//!
//! The models are not shipped: the user installs AlphaPeptDeep's released
//! model archive once (`tuffold models --install`), whose `generic` MS2 and
//! RT weights are converted from PyTorch's format to safetensors. The layers are
//! written out here to match the PyTorch originals; they were checked
//! against a PyTorch implementation (fragment predictions correlate with
//! measured spectra at a median 0.93 on PXD006447, the LSTM matches
//! PyTorch's to 1e-6).

mod encode;
mod features;
mod ms2;
mod nn;
mod rt;

use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use candle_core::Device;

pub use encode::{Encoded, encode};
pub use features::{FEATURE_NAMES, MIN_FIT, Ms2Scorer, Predictions, PsmSpectrum, Tolerance};
pub use ms2::Instrument;
pub use rt::DeepRetention;

/// AlphaPeptDeep's model archive (`pretrained_models_v3.zip`).
pub const MODELS_URL: &str = "https://github.com/MannLabs/alphapeptdeep/releases/download/pre-trained-models/pretrained_models_v3.zip";
/// Converted from the archive's `generic/ms2.pth` and `generic/rt.pth` at
/// install, so everyday runs never parse a pickle.
const MS2_FILE: &str = "ms2.safetensors";
const RT_FILE: &str = "rt.safetensors";
/// Largest model file an archive may hold (the real ones are 16 and 3 MB).
const MAX_MODEL_FILE: u64 = 256 << 20;

#[derive(Debug, thiserror::Error)]
pub enum PredictError {
    #[error("{file}: {problem}")]
    Weights { file: String, problem: String },
    #[error("{0}")]
    Archive(String),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// The installed models.
pub struct Models {
    ms2: ms2::Ms2Model,
    rt: nn::Weights,
    device: Device,
}

impl Models {
    /// Where models are installed: `TUFFOLD_MODELS`, else the user's data
    /// folder (`~/Library/Application Support/Tuffold/models` on macOS,
    /// `$XDG_DATA_HOME/tuffold/models` or `~/.local/share/tuffold/models`
    /// elsewhere).
    pub fn default_dir() -> Option<PathBuf> {
        if let Some(dir) = std::env::var_os("TUFFOLD_MODELS") {
            return Some(PathBuf::from(dir));
        }
        let home = std::env::var_os("HOME").map(PathBuf::from);
        if cfg!(target_os = "macos") {
            return home.map(|home| home.join("Library/Application Support/Tuffold/models"));
        }
        std::env::var_os("XDG_DATA_HOME")
            .map(PathBuf::from)
            .or_else(|| home.map(|home| home.join(".local/share")))
            .map(|data| data.join("tuffold/models"))
    }

    /// Copy the generic MS2 and RT models out of AlphaPeptDeep's archive
    /// into `dir`, and check they load.
    pub fn install(archive: &Path, dir: &Path) -> Result<(), PredictError> {
        let file = std::fs::File::open(archive)?;
        let mut zip = zip::ZipArchive::new(file)
            .map_err(|error| PredictError::Archive(format!("{}: {error}", archive.display())))?;
        std::fs::create_dir_all(dir)?;
        for (wanted, converted) in [("ms2.pth", MS2_FILE), ("rt.pth", RT_FILE)] {
            let name = (0..zip.len())
                .filter_map(|i| zip.name_for_index(i).map(str::to_string))
                .find(|name| name.ends_with(&format!("generic/{wanted}")))
                .ok_or_else(|| {
                    PredictError::Archive(format!(
                        "{} has no generic/{wanted}; is it AlphaPeptDeep's model archive?",
                        archive.display()
                    ))
                })?;
            let mut entry = zip
                .by_name(&name)
                .map_err(|error| PredictError::Archive(format!("{name}: {error}")))?;
            let mut bytes = Vec::new();
            (&mut entry)
                .take(MAX_MODEL_FILE + 1)
                .read_to_end(&mut bytes)?;
            if bytes.len() as u64 > MAX_MODEL_FILE {
                return Err(PredictError::Archive(format!(
                    "{name} is too large to be a model"
                )));
            }
            let unpacked = dir.join(wanted);
            std::fs::write(&unpacked, bytes)?;
            let weights = nn::Weights::from_pth(&unpacked);
            std::fs::remove_file(&unpacked)?;
            weights?.save(&dir.join(converted))?;
        }
        Self::load(dir)?
            .map(|_| ())
            .ok_or_else(|| PredictError::Archive("the installed models don't load".into()))
    }

    /// `None` when no models are installed in `dir`.
    pub fn load(dir: &Path) -> Result<Option<Self>, PredictError> {
        let (ms2, rt) = (dir.join(MS2_FILE), dir.join(RT_FILE));
        if !ms2.exists() || !rt.exists() {
            return Ok(None);
        }
        let device = Device::Cpu;
        Ok(Some(Self {
            ms2: ms2::Ms2Model::load(&ms2, &device)?,
            rt: nn::Weights::load(&rt)?,
            device,
        }))
    }

    /// Fragment-intensity scoring at the settings that fit `sample` best;
    /// `None` when the model doesn't fit these spectra.
    pub fn calibrate(&self, sample: &[PsmSpectrum<'_>]) -> Option<Ms2Scorer<'_>> {
        features::calibrate(&self.ms2, sample)
    }

    /// Predicted b⁺, b²⁺, y⁺, y²⁺ intensities of equal-length peptides, one
    /// row per cleavage site.
    pub fn fragment_intensities(
        &self,
        peptides: &[&Encoded],
        charge: u8,
        nce: f64,
        instrument: Instrument,
    ) -> Result<Vec<Vec<[f32; 4]>>, PredictError> {
        self.ms2
            .predict(peptides, charge, nce, instrument)
            .map_err(|error| PredictError::Weights {
                file: MS2_FILE.into(),
                problem: error.to_string(),
            })
    }

    /// The pretrained model's gradient-relative retention times of
    /// equal-length peptides (before any fine-tuning).
    pub fn pretrained_retention(&self, peptides: &[&Encoded]) -> Result<Vec<f32>, PredictError> {
        rt::pretrained(&self.rt, peptides, &self.device).map_err(|error| PredictError::Weights {
            file: RT_FILE.into(),
            problem: error.to_string(),
        })
    }

    /// A retention-time predictor for these peptides, fine-tuned per run by
    /// validation.
    pub fn retention(&self, peptides: HashMap<i64, Encoded>) -> DeepRetention {
        DeepRetention::new(self.rt.clone(), peptides, self.device.clone())
    }
}
