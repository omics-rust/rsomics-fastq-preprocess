#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Coherent single-pass FASTQ trimming and filtering operations.

mod input;
mod output;
/// Streaming SE/PE execution and typed reports.
pub mod pipeline;
/// Pure trimming and filtering transforms.
pub mod transform;

pub use pipeline::{
    FilterMetrics, IoSpec, Mode, Operation, PipelineConfig, PreprocessReport, TrimMetrics, execute,
};
pub use transform::filter::{FilterConfig, FilterOutcome, PhredEncoding};
pub use transform::trim::{FixedTrim, PolyTailConfig, TrimConfig};
