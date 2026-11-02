//! AlphaPeptDeep's retention-time model (generic, `rt.pth`): residues
//! (one-hot) and modifications through three convolutions, a 2-layer
//! bidirectional LSTM and an attention-weighted sum, then a linear decoder.
//!
//! The pretrained model predicts a gradient-relative time; each run's
//! gradient differs, so its attention pooling and decoder are fine-tuned per
//! run (and per rescoring fold) on that run's confident peptides, over the
//! frozen convolutions and LSTM, whose outputs are computed once per
//! peptide. On PXD006447 the median error falls from 2.5 min (pretrained,
//! calibrated) to 1.5 min, against 4.5–5 min for the classical model.
//! Fine-tuning every layer reached 1.2 min but took minutes per fold (the
//! unrolled LSTM's backward pass); the head alone takes seconds.

use std::collections::HashMap;

use candle_core::{D, DType, Device, Result, Tensor, Var};
use candle_nn::optim::{AdamW, Optimizer, ParamsAdamW};
use tuffold_validate::{RetentionPredictor, TrainedRetention};

use crate::encode::Encoded;
use crate::nn::{Lookup, Weights, batch_tensors, decoder, linear, modification_embedding};

const RESIDUE_SLOTS: usize = 27;
const LSTM_LAYERS: usize = 2;
const EPOCHS: usize = 30;
const BATCH: usize = 64;
const PREDICT_BATCH: usize = 1024;
const LEARNING_RATE: f64 = 3e-3;
/// The parameters fine-tuned: attention pooling and decoder.
const HEAD: [&str; 2] = ["rt_encoder.attn_sum.", "rt_decoder."];
/// Fixed so fine-tuning is deterministic.
const SEED: u64 = 0x5eed_0f7f_f01d;

/// Parameters, fixed or being trained.
struct Params {
    tensors: HashMap<String, Tensor>,
}

impl Lookup for Params {
    fn at(&self, name: &str) -> Result<&Tensor> {
        self.tensors
            .get(name)
            .ok_or_else(|| candle_core::Error::Msg(format!("weight {name} is missing")))
    }
}

/// Fine-tunes the pretrained model on a run's peptides.
pub struct DeepRetention {
    base: Weights,
    peptides: HashMap<i64, Encoded>,
    device: Device,
}

impl DeepRetention {
    pub(crate) fn new(base: Weights, peptides: HashMap<i64, Encoded>, device: Device) -> Self {
        Self {
            base,
            peptides,
            device,
        }
    }

    /// The fine-tuned head, and its predictions for the training peptides
    /// (their encoder outputs are at hand).
    fn fine_tune(&self, training: &[(i64, &Encoded, f32)]) -> Result<(Params, HashMap<i64, f64>)> {
        let mut vars = HashMap::new();
        let mut fixed = HashMap::new();
        for (name, tensor) in self.base.tensors() {
            if HEAD.iter().any(|head| name.starts_with(head)) {
                vars.insert(name.clone(), Var::from_tensor(tensor)?);
            } else {
                fixed.insert(name.clone(), tensor.clone());
            }
        }
        let mut optimizer = AdamW::new(
            vars.values().cloned().collect(),
            ParamsAdamW {
                lr: LEARNING_RATE,
                weight_decay: 0.0,
                ..ParamsAdamW::default()
            },
        )?;
        // The frozen encoder's outputs, once per batch of equal lengths.
        let order: Vec<usize> = (0..training.len()).collect();
        let mut encoded_batches = Vec::new();
        let mut batch_ids = Vec::new();
        for batch in batches(&order, |index| training[index].1.len(), BATCH) {
            let peptides: Vec<&Encoded> = batch.iter().map(|&i| training[i].1).collect();
            let targets: Vec<f32> = batch.iter().map(|&i| training[i].2).collect();
            encoded_batches.push((
                encoder(&self.base, &peptides, &self.device)?,
                Tensor::from_vec(targets, batch.len(), &self.device)?,
            ));
            batch_ids.push(batch.iter().map(|&i| training[i].0).collect::<Vec<i64>>());
        }
        let mut schedule: Vec<usize> = (0..encoded_batches.len()).collect();
        let mut random = SEED;
        for _ in 0..EPOCHS {
            shuffle(&mut schedule, &mut random);
            for &index in &schedule {
                let (hidden, targets) = &encoded_batches[index];
                let params = Params {
                    tensors: vars
                        .iter()
                        .map(|(name, var)| (name.clone(), var.as_tensor().clone()))
                        .collect(),
                };
                let loss = (head(&params, hidden)? - targets)?.abs()?.mean_all()?;
                optimizer.backward_step(&loss)?;
            }
        }
        let params = Params {
            tensors: vars
                .into_iter()
                .map(|(name, var)| (name, var.as_tensor().detach()))
                .chain(fixed)
                .collect(),
        };
        let mut known = HashMap::with_capacity(training.len());
        for ((hidden, _), ids) in encoded_batches.iter().zip(&batch_ids) {
            let values = head(&params, hidden)?
                .to_dtype(DType::F32)?
                .to_vec1::<f32>()?;
            known.extend(ids.iter().zip(values).map(|(&id, v)| (id, f64::from(v))));
        }
        Ok((params, known))
    }
}

impl RetentionPredictor for DeepRetention {
    fn train(&self, training: &[(i64, f64)]) -> Option<Box<dyn TrainedRetention + '_>> {
        let scale = training.iter().map(|(_, time)| *time).fold(0.0, f64::max);
        let examples: Vec<(i64, &Encoded, f32)> = training
            .iter()
            .filter_map(|&(id, time)| Some((id, self.peptides.get(&id)?, (time / scale) as f32)))
            .collect();
        if scale <= 0.0 || examples.is_empty() {
            return None;
        }
        let (params, known) = self.fine_tune(&examples).ok()?;
        Some(Box::new(Trained {
            params,
            known,
            owner: self,
        }))
    }
}

struct Trained<'a> {
    params: Params,
    known: HashMap<i64, f64>,
    owner: &'a DeepRetention,
}

impl TrainedRetention for Trained<'_> {
    fn predict(&self, peptides: &[i64]) -> Vec<Option<f64>> {
        // Rows repeat peptides (one per PSM): predict each once.
        let mut unique: Vec<i64> = peptides
            .iter()
            .copied()
            .filter(|id| !self.known.contains_key(id) && self.owner.peptides.contains_key(id))
            .collect();
        unique.sort_unstable();
        unique.dedup();
        let encoded = |i: usize| &self.owner.peptides[&unique[i]];
        let order: Vec<usize> = (0..unique.len()).collect();
        let mut predicted = self.known.clone();
        for batch in batches(&order, |i| encoded(i).len(), PREDICT_BATCH) {
            let inputs: Vec<&Encoded> = batch.iter().map(|&i| encoded(i)).collect();
            let Ok(values) = encoder(&self.owner.base, &inputs, &self.owner.device)
                .and_then(|hidden| head(&self.params, &hidden))
                .and_then(|t| t.to_dtype(DType::F32)?.to_vec1::<f32>())
            else {
                continue;
            };
            for (&i, value) in batch.iter().zip(values) {
                predicted.insert(unique[i], f64::from(value));
            }
        }
        peptides
            .iter()
            .map(|id| predicted.get(id).copied())
            .collect()
    }
}

pub(crate) fn pretrained(
    weights: &Weights,
    batch: &[&Encoded],
    device: &Device,
) -> Result<Vec<f32>> {
    head(weights, &encoder(weights, batch, device)?)?
        .to_dtype(DType::F32)?
        .to_vec1()
}

/// Convolutions and LSTM on equal-length peptides: `[B, L, 256]`.
fn encoder(p: &impl Lookup, batch: &[&Encoded], device: &Device) -> Result<Tensor> {
    let (residues, modifications) = batch_tensors(batch, device)?;
    let one_hot = candle_nn::encoding::one_hot(residues, RESIDUE_SLOTS, 1f32, 0f32)?;
    let modification_x =
        modification_embedding(&modifications, p.at("rt_encoder.mod_nn.nn.weight")?)?;
    let x = Tensor::cat(&[one_hot, modification_x], 2)?
        .transpose(1, 2)?
        .contiguous()?;
    let mut parts = vec![x.clone()];
    for (name, padding) in [("short", 1), ("medium", 2), ("long", 3)] {
        let bias = p.at(&format!("rt_encoder.input_cnn.cnn_{name}.bias"))?;
        let conv = x
            .conv1d(
                p.at(&format!("rt_encoder.input_cnn.cnn_{name}.weight"))?,
                padding,
                1,
                1,
                1,
            )?
            .broadcast_add(&bias.reshape((1, (), 1))?)?;
        parts.push(conv);
    }
    let x = Tensor::cat(&parts, 1)?.transpose(1, 2)?.contiguous()?;
    lstm(p, &x)
}

/// Attention-weighted sum over positions, then the decoder: one prediction
/// per peptide.
fn head(p: &impl Lookup, x: &Tensor) -> Result<Tensor> {
    let attention = candle_nn::ops::softmax(
        &linear(x, p.at("rt_encoder.attn_sum.attn.0.weight")?, None)?,
        1,
    )?;
    let pooled = x.broadcast_mul(&attention)?.sum(1)?;
    decoder(&pooled, p, "rt_decoder.")?.squeeze(D::Minus1)
}

/// Two bidirectional LSTM layers from learned initial states (PyTorch's gate
/// order: input, forget, cell, output).
fn lstm(p: &impl Lookup, x: &Tensor) -> Result<Tensor> {
    let prefix = "rt_encoder.hidden_nn.rnn.";
    let (size, length, _) = x.dims3()?;
    let h0 = p.at("rt_encoder.hidden_nn.rnn_h0")?;
    let c0 = p.at("rt_encoder.hidden_nn.rnn_c0")?;
    let mut input = x.clone();
    for layer in 0..LSTM_LAYERS {
        let mut directions = Vec::with_capacity(2);
        for (direction, suffix) in ["", "_reverse"].into_iter().enumerate() {
            let state = layer * 2 + direction;
            let projected = linear(
                &input,
                p.at(&format!("{prefix}weight_ih_l{layer}{suffix}"))?,
                Some(p.at(&format!("{prefix}bias_ih_l{layer}{suffix}"))?),
            )?;
            let recurrent = p.at(&format!("{prefix}weight_hh_l{layer}{suffix}"))?;
            let recurrent_bias = p.at(&format!("{prefix}bias_hh_l{layer}{suffix}"))?;
            let hidden_size = recurrent.dim(1)?;
            let mut h = h0
                .get(state)?
                .broadcast_as((size, hidden_size))?
                .contiguous()?;
            let mut c = c0
                .get(state)?
                .broadcast_as((size, hidden_size))?
                .contiguous()?;
            let mut outputs = vec![None; length];
            let steps: Vec<usize> = if direction == 0 {
                (0..length).collect()
            } else {
                (0..length).rev().collect()
            };
            for t in steps {
                let gates = (projected.narrow(1, t, 1)?.squeeze(1)?
                    + linear(&h, recurrent, Some(recurrent_bias))?)?;
                let chunks = gates.chunk(4, 1)?;
                let (i, f, g, o) = (
                    sigmoid(&chunks[0])?,
                    sigmoid(&chunks[1])?,
                    chunks[2].tanh()?,
                    sigmoid(&chunks[3])?,
                );
                c = ((f * &c)? + (i * g)?)?;
                h = (o * c.tanh()?)?;
                outputs[t] = Some(h.clone());
            }
            let outputs: Vec<Tensor> = outputs.into_iter().flatten().collect();
            directions.push(Tensor::stack(&outputs, 1)?);
        }
        input = Tensor::cat(&directions, 2)?;
    }
    Ok(input)
}

/// Written out so it differentiates everywhere.
fn sigmoid(x: &Tensor) -> Result<Tensor> {
    (x.neg()?.exp()? + 1.0)?.recip()
}

/// Indices grouped by peptide length (the models take equal lengths per
/// batch), in order within each length, cut into batches.
fn batches(order: &[usize], length: impl Fn(usize) -> usize, size: usize) -> Vec<Vec<usize>> {
    let mut by_length: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for &index in order {
        by_length.entry(length(index)).or_default().push(index);
    }
    let mut out: Vec<Vec<usize>> = by_length
        .into_values()
        .flat_map(|group| {
            group
                .chunks(size)
                .map(<[usize]>::to_vec)
                .collect::<Vec<_>>()
        })
        .collect();
    // Interleave lengths across an epoch rather than training on one length
    // at a time (the order list was shuffled; keep its first index's order).
    out.sort_by_key(|batch| order.iter().position(|&i| i == batch[0]));
    out
}

fn shuffle(order: &mut [usize], state: &mut u64) {
    for i in (1..order.len()).rev() {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        let j = (*state % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }
}
