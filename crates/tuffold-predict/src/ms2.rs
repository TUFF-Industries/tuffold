//! AlphaPeptDeep's MS2 model (generic, `ms2.pth`): a 4-layer BERT encoder
//! over residues, modifications and position, with the precursor charge,
//! collision energy and instrument appended to every position; a linear
//! decoder gives each cleavage site's b⁺, b²⁺, y⁺ and y²⁺ intensities
//! (relative, the strongest about 1). The neutral-loss branch is not run.

use std::path::Path;

use candle_core::{Device, Result, Tensor};

use crate::PredictError;
use crate::encode::Encoded;
use crate::nn::{
    Weights, batch_tensors, decoder, gelu, layer_norm, linear, modification_embedding,
};

const HEADS: usize = 8;
/// AlphaPeptDeep's BERT layer-norm epsilon.
const LAYER_NORM_EPS: f64 = 1e-8;
const CHARGE_FACTOR: f64 = 0.1;
const NCE_FACTOR: f64 = 0.01;
/// The instrument one-hot's width.
const INSTRUMENT_SLOTS: usize = 8;
/// b⁺, b²⁺, y⁺, y²⁺ (the model's own order).
pub const CHANNELS: usize = 4;

/// Instruments the model was trained on, by its index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Instrument {
    QExactive = 0,
    Lumos = 1,
    TimsTof = 2,
    SciexTof = 3,
}

pub(crate) struct Ms2Model {
    weights: Weights,
    layers: usize,
    device: Device,
}

impl Ms2Model {
    pub fn load(path: &Path, device: &Device) -> std::result::Result<Self, PredictError> {
        Self::new(Weights::load(path)?, path, device)
    }

    pub fn new(
        weights: Weights,
        path: &Path,
        device: &Device,
    ) -> std::result::Result<Self, PredictError> {
        let layers = weights
            .names()
            .filter_map(|name| {
                name.strip_prefix("hidden_nn.bert.layer.")?
                    .split('.')
                    .next()?
                    .parse::<usize>()
                    .ok()
            })
            .max()
            .map_or(0, |last| last + 1);
        if layers == 0 {
            return Err(PredictError::Weights {
                file: path.display().to_string(),
                problem: "no BERT layers".into(),
            });
        }
        Ok(Self {
            weights,
            layers,
            device: device.clone(),
        })
    }

    /// Intensities for equal-length peptides at one charge: per peptide, one
    /// row per cleavage site (`b₁/y₍ₙ₋₁₎` first) of [`CHANNELS`] values.
    pub fn predict(
        &self,
        batch: &[&Encoded],
        charge: u8,
        nce: f64,
        instrument: Instrument,
    ) -> Result<Vec<Vec<[f32; CHANNELS]>>> {
        let w = &self.weights;
        let size = batch.len();
        let (residues, modifications) = batch_tensors(batch, &self.device)?;
        let length = residues.dim(1)?;
        let residue_x = w
            .get("input_nn.aa_emb.weight")?
            .index_select(&residues.flatten_all()?, 0)?
            .reshape((size, length, ()))?;
        let modification_x =
            modification_embedding(&modifications, w.get("input_nn.mod_nn.nn.weight")?)?;
        let x = Tensor::cat(&[residue_x, modification_x], 2)?
            .broadcast_add(&w.get("input_nn.pos_encoder.pe")?.narrow(1, 0, length)?)?;
        let mut meta = vec![0f32; INSTRUMENT_SLOTS + 1];
        meta[instrument as usize] = 1.0;
        meta[INSTRUMENT_SLOTS] = (nce * NCE_FACTOR) as f32;
        let meta = linear(
            &Tensor::from_vec(meta, (1, INSTRUMENT_SLOTS + 1), &self.device)?,
            w.get("meta_nn.nn.weight")?,
            Some(w.get("meta_nn.nn.bias")?),
        )?;
        let charge = Tensor::full(
            (f64::from(charge) * CHARGE_FACTOR) as f32,
            (1, 1),
            &self.device,
        )?;
        let meta = Tensor::cat(&[meta, charge], 1)?
            .unsqueeze(1)?
            .broadcast_as((size, length, 8))?
            .contiguous()?;
        let input = Tensor::cat(&[x, meta], 2)?;
        let mut hidden = input.clone();
        for layer in 0..self.layers {
            hidden = bert_layer(&hidden, w, &format!("hidden_nn.bert.layer.{layer}."))?;
        }
        let hidden = (hidden + (input * 0.2)?)?;
        let out = decoder(&hidden, w, "output_nn.")?
            .narrow(1, 3, length - 3)?
            .to_dtype(candle_core::DType::F32)?
            .to_vec3::<f32>()?;
        Ok(out
            .into_iter()
            .map(|rows| {
                rows.into_iter()
                    .map(|row| [row[0], row[1], row[2], row[3]])
                    .collect()
            })
            .collect())
    }
}

/// One HuggingFace `BertLayer` (post-norm, exact GELU), no attention mask.
fn bert_layer(x: &Tensor, w: &Weights, prefix: &str) -> Result<Tensor> {
    let (size, length, hidden) = x.dims3()?;
    let head = hidden / HEADS;
    let project = |name: &str| -> Result<Tensor> {
        linear(
            x,
            w.get(&format!("{prefix}attention.self.{name}.weight"))?,
            Some(w.get(&format!("{prefix}attention.self.{name}.bias"))?),
        )?
        .reshape((size, length, HEADS, head))?
        .transpose(1, 2)?
        .contiguous()
    };
    let (query, key, value) = (project("query")?, project("key")?, project("value")?);
    let scores = (query.matmul(&key.t()?.contiguous()?)? / (head as f64).sqrt())?;
    let attention = candle_nn::ops::softmax_last_dim(&scores)?;
    let context = attention
        .matmul(&value)?
        .transpose(1, 2)?
        .reshape((size, length, hidden))?;
    let attended = linear(
        &context,
        w.get(&format!("{prefix}attention.output.dense.weight"))?,
        Some(w.get(&format!("{prefix}attention.output.dense.bias"))?),
    )?;
    let attended = layer_norm(
        &(attended + x)?,
        w.get(&format!("{prefix}attention.output.LayerNorm.weight"))?,
        w.get(&format!("{prefix}attention.output.LayerNorm.bias"))?,
        LAYER_NORM_EPS,
    )?;
    let intermediate = linear(
        &attended,
        w.get(&format!("{prefix}intermediate.dense.weight"))?,
        Some(w.get(&format!("{prefix}intermediate.dense.bias"))?),
    )?;
    let intermediate = gelu(&intermediate)?;
    let output = linear(
        &intermediate,
        w.get(&format!("{prefix}output.dense.weight"))?,
        Some(w.get(&format!("{prefix}output.dense.bias"))?),
    )?;
    layer_norm(
        &(output + attended)?,
        w.get(&format!("{prefix}output.LayerNorm.weight"))?,
        w.get(&format!("{prefix}output.LayerNorm.bias"))?,
        LAYER_NORM_EPS,
    )
}
