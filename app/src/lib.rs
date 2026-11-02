//! Tuffold desktop application: GPUI views over the headless engine crates.
//!
//! The binary (`main.rs`) only boots the platform; everything here is also
//! driven headlessly by the snapshot tests.

pub mod browse;
pub mod project;
pub mod quant;
pub mod sorted;
pub mod theme;
mod views;
pub mod workspace;
