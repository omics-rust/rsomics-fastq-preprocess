//! Pure per-record preprocessing transforms.

/// Whole-read fastp-style filters.
pub mod filter;
/// Fastp-compatible rolling quality cuts.
pub mod quality;
/// Fixed and poly-tail trimming.
pub mod trim;
