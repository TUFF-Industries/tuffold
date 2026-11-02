//! Layers shared by the models, written out on candle tensors so they match
//! the PyTorch originals exactly (checked against them; see `lib.rs`).

use std::collections::HashMap;
use std::path::Path;

use candle_core::{Device, Result, Tensor};

use crate::PredictError;

/// A model's weights by PyTorch parameter name.
#[derive(Clone)]
pub(crate) struct Weights(HashMap<String, Tensor>);

impl Weights {
    /// From a PyTorch `.pth` file (installing only: it parses a pickle).
    pub fn from_pth(path: &Path) -> std::result::Result<Self, PredictError> {
        let tensors =
            candle_core::pickle::read_all(path).map_err(|error| PredictError::Weights {
                file: path.display().to_string(),
                problem: error.to_string(),
            })?;
        Ok(Self(tensors.into_iter().collect()))
    }

    pub fn load(path: &Path) -> std::result::Result<Self, PredictError> {
        candle_core::safetensors::load(path, &Device::Cpu)
            .map(Self)
            .map_err(|error| PredictError::Weights {
                file: path.display().to_string(),
                problem: error.to_string(),
            })
    }

    pub fn save(&self, path: &Path) -> std::result::Result<(), PredictError> {
        candle_core::safetensors::save(&self.0, path).map_err(|error| PredictError::Weights {
            file: path.display().to_string(),
            problem: error.to_string(),
        })
    }

    pub fn get(&self, name: &str) -> Result<&Tensor> {
        self.0
            .get(name)
            .ok_or_else(|| candle_core::Error::Msg(format!("weight {name} is missing")))
    }

    pub fn names(&self) -> impl Iterator<Item = &String> {
        self.0.keys()
    }

    pub fn tensors(&self) -> impl Iterator<Item = (&String, &Tensor)> {
        self.0.iter()
    }
}

/// `x · Wᵀ + b` over the last dimension, for any leading dimensions.
pub(crate) fn linear(x: &Tensor, weight: &Tensor, bias: Option<&Tensor>) -> Result<Tensor> {
    let dims = x.dims().to_vec();
    let inputs = *dims.last().expect("a dimension");
    let rows = x.elem_count() / inputs;
    let out = x.reshape((rows, inputs))?.matmul(&weight.t()?)?;
    let out = match bias {
        Some(bias) => out.broadcast_add(bias)?,
        None => out,
    };
    let mut shape = dims;
    *shape.last_mut().expect("a dimension") = weight.dim(0)?;
    out.reshape(shape)
}

pub(crate) fn layer_norm(x: &Tensor, weight: &Tensor, bias: &Tensor, eps: f64) -> Result<Tensor> {
    candle_nn::ops::layer_norm(&x.contiguous()?, weight, bias, eps as f32)
}

/// Exact (erf) GELU, computed on all cores: candle's own is single-threaded
/// and was most of the MS2 model's time. erf from Abramowitz & Stegun
/// 7.1.26 (absolute error below 1.5e-7).
pub(crate) fn gelu(x: &Tensor) -> Result<Tensor> {
    let shape = x.shape().clone();
    let mut values = x.flatten_all()?.to_vec1::<f32>()?;
    let threads = std::thread::available_parallelism().map_or(1, usize::from);
    let chunk = values.len().div_ceil(threads).max(4096);
    std::thread::scope(|scope| {
        for part in values.chunks_mut(chunk) {
            scope.spawn(move || {
                for v in part {
                    *v = 0.5 * *v * (1.0 + erf(*v * std::f32::consts::FRAC_1_SQRT_2));
                }
            });
        }
    });
    Tensor::from_vec(values, shape, x.device())
}

fn erf(x: f32) -> f32 {
    let t = 1.0 / (1.0 + 0.327_591_1 * x.abs());
    let poly = t
        * (0.254_829_6
            + t * (-0.284_496_74 + t * (1.421_413_7 + t * (-1.453_152_1 + t * 1.061_405_4))));
    let y = 1.0 - poly * (-x * x).exp();
    y.copysign(x)
}

/// AlphaPeptDeep's `Decoder_Linear`: linear, PReLU (one slope), linear.
pub(crate) fn decoder(x: &Tensor, weights: &impl Lookup, prefix: &str) -> Result<Tensor> {
    let hidden = linear(
        x,
        weights.at(&format!("{prefix}nn.0.weight"))?,
        Some(weights.at(&format!("{prefix}nn.0.bias"))?),
    )?;
    let positive = hidden.relu()?;
    let negative = (&hidden - &positive)?;
    let hidden =
        (positive + negative.broadcast_mul(weights.at(&format!("{prefix}nn.1.weight"))?)?)?;
    linear(
        &hidden,
        weights.at(&format!("{prefix}nn.2.weight"))?,
        Some(weights.at(&format!("{prefix}nn.2.bias"))?),
    )
}

/// `Mod_Embedding_FixFirstK`: the first six elements (C H N O P S) as they
/// are, the rest through one bias-free linear layer.
pub(crate) fn modification_embedding(modifications: &Tensor, weight: &Tensor) -> Result<Tensor> {
    let first = modifications.narrow(2, 0, 6)?;
    let rest = modifications.narrow(2, 6, modifications.dim(2)? - 6)?;
    Tensor::cat(&[first, linear(&rest, weight, None)?], 2)
}

/// Read access to parameters, whether fixed or being trained.
pub(crate) trait Lookup {
    fn at(&self, name: &str) -> Result<&Tensor>;
}

impl Lookup for Weights {
    fn at(&self, name: &str) -> Result<&Tensor> {
        self.get(name)
    }
}

/// A batch of equal-length peptides as tensors: residue indices `[B, L]`
/// and modification compositions `[B, L, 109]`.
pub(crate) fn batch_tensors(
    batch: &[&crate::encode::Encoded],
    device: &Device,
) -> Result<(Tensor, Tensor)> {
    let length = batch[0].len();
    let residues: Vec<u32> = batch
        .iter()
        .flat_map(|p| p.residues.iter().copied())
        .collect();
    let modifications: Vec<f32> = batch
        .iter()
        .flat_map(|p| p.modifications.iter().copied())
        .collect();
    Ok((
        Tensor::from_vec(residues, (batch.len(), length), device)?,
        Tensor::from_vec(
            modifications,
            (batch.len(), length, crate::encode::MOD_FEATURES),
            device,
        )?,
    ))
}
