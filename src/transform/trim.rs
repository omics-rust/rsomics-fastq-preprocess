use std::num::NonZeroUsize;

use rsomics_common::{Result, RsomicsError};
use rsomics_seqio::OwnedRecord;

/// Fixed 5-prime and 3-prime trimming.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct FixedTrim {
    /// Bases removed from the 5-prime end.
    pub front: usize,
    /// Bases removed from the 3-prime end.
    pub tail: usize,
}

/// Poly-G or dominant poly-X trimming policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PolyTailConfig {
    /// Minimum accepted tail span.
    pub min_len: usize,
    /// Absolute mismatch cap.
    pub max_mismatches: usize,
    /// One mismatch is allowed per this many scanned bases.
    pub mismatch_per_bases: NonZeroUsize,
}

impl Default for PolyTailConfig {
    fn default() -> Self {
        Self {
            min_len: 10,
            max_mismatches: 5,
            mismatch_per_bases: NonZeroUsize::new(8).expect("8 is nonzero"),
        }
    }
}

/// Complete trimming policy for both mates.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TrimConfig {
    /// Fixed trimming for R1 or single-end reads.
    pub fixed_r1: FixedTrim,
    /// Fixed trimming for R2.
    pub fixed_r2: FixedTrim,
    /// Optional poly-G trimming.
    pub poly_g: Option<PolyTailConfig>,
    /// Optional dominant poly-X trimming.
    pub poly_x: Option<PolyTailConfig>,
}

impl TrimConfig {
    pub(crate) fn validate(&self) -> Result<()> {
        validate_poly_tails(self.poly_g, self.poly_x)
    }
}

/// Per-record trimming counters.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RecordTrimMetrics {
    /// Bases removed by fixed trimming.
    pub fixed_bases: u64,
    /// Bases removed by poly-G trimming.
    pub poly_g_bases: u64,
    /// Whether a poly-G tail matched.
    pub poly_g_matched: bool,
    /// Bases removed by dominant poly-X trimming.
    pub poly_x_bases: u64,
    /// Whether a poly-X tail matched.
    pub poly_x_matched: bool,
}

/// Applies fixed, poly-G, then poly-X trimming.
///
/// Returns an error when the record is not a complete FASTQ record or a
/// poly-tail minimum is zero.
pub fn trim_record(
    record: &mut OwnedRecord,
    fixed: FixedTrim,
    poly_g: Option<PolyTailConfig>,
    poly_x: Option<PolyTailConfig>,
) -> Result<RecordTrimMetrics> {
    let quality_len = record
        .qual
        .as_ref()
        .ok_or_else(|| RsomicsError::InvalidInput("FASTQ record requires quality scores".into()))?
        .len();
    if record.seq.len() != quality_len {
        return Err(RsomicsError::InvalidInput(format!(
            "FASTQ sequence/quality length mismatch: {} vs {quality_len}",
            record.seq.len()
        )));
    }
    validate_poly_tails(poly_g, poly_x)?;
    Ok(trim_seqio_record(record, fixed, poly_g, poly_x))
}

pub(crate) fn trim_seqio_record(
    record: &mut OwnedRecord,
    fixed: FixedTrim,
    poly_g: Option<PolyTailConfig>,
    poly_x: Option<PolyTailConfig>,
) -> RecordTrimMetrics {
    let mut metrics = RecordTrimMetrics::default();
    apply_fixed(record, fixed, &mut metrics);
    if let Some(config) = poly_g {
        apply_poly_g(record, config, &mut metrics);
    }
    if let Some(config) = poly_x {
        apply_poly_x(record, config, &mut metrics);
    }
    metrics
}

fn validate_poly_tails(
    poly_g: Option<PolyTailConfig>,
    poly_x: Option<PolyTailConfig>,
) -> Result<()> {
    if poly_g.is_some_and(|config| config.min_len == 0) {
        return Err(RsomicsError::ConfigError(
            "poly-G minimum tail length must be positive".into(),
        ));
    }
    if poly_x.is_some_and(|config| config.min_len == 0) {
        return Err(RsomicsError::ConfigError(
            "poly-X minimum tail length must be positive".into(),
        ));
    }
    Ok(())
}

fn apply_fixed(record: &mut OwnedRecord, config: FixedTrim, metrics: &mut RecordTrimMetrics) {
    let original = record.seq.len();
    let start = config.front.min(original);
    let end = original.saturating_sub(config.tail).max(start);
    replace_range(record, start, end);
    metrics.fixed_bases = metric_count(original - record.seq.len());
}

fn apply_poly_g(record: &mut OwnedRecord, config: PolyTailConfig, metrics: &mut RecordTrimMetrics) {
    if let Some(cut) = find_poly_g_3p(&record.seq, config) {
        let removed = record.seq.len() - cut;
        truncate(record, cut);
        metrics.poly_g_matched = true;
        metrics.poly_g_bases = metric_count(removed);
    }
}

fn apply_poly_x(record: &mut OwnedRecord, config: PolyTailConfig, metrics: &mut RecordTrimMetrics) {
    if let Some(cut) = find_poly_x_3p(&record.seq, config) {
        let removed = record.seq.len() - cut;
        truncate(record, cut);
        metrics.poly_x_matched = true;
        metrics.poly_x_bases = metric_count(removed);
    }
}

fn replace_range(record: &mut OwnedRecord, start: usize, end: usize) {
    if start > 0 {
        record.seq.drain(..start);
        if let Some(quality) = record.qual.as_mut() {
            quality.drain(..start);
        }
    }
    truncate(record, end - start);
}

fn truncate(record: &mut OwnedRecord, length: usize) {
    record.seq.truncate(length);
    if let Some(quality) = record.qual.as_mut() {
        quality.truncate(length);
    }
}

/// Exact translation of fastp v1.3.6 `PolyX::trimPolyG`.
fn find_poly_g_3p(seq: &[u8], config: PolyTailConfig) -> Option<usize> {
    if seq.is_empty() {
        return None;
    }
    let min_len = config.min_len.max(1);
    let mut mismatches = 0usize;
    let mut scanned_index = 0usize;
    let mut first_g_position = seq.len() - 1;
    while scanned_index < seq.len() {
        let position = seq.len() - scanned_index - 1;
        if seq[position] == b'G' {
            first_g_position = position;
        } else {
            mismatches += 1;
        }
        let allowed = (scanned_index + 1) / config.mismatch_per_bases.get();
        if mismatches > config.max_mismatches
            || (mismatches > allowed && scanned_index >= min_len - 1)
        {
            break;
        }
        scanned_index += 1;
    }
    (scanned_index >= min_len).then_some(first_g_position)
}

/// Exact safe translation of fastp v1.3.6 `PolyX::trimPolyX`.
fn find_poly_x_3p(seq: &[u8], config: PolyTailConfig) -> Option<usize> {
    if seq.is_empty() {
        return None;
    }
    let min_len = config.min_len.max(1);
    // fastp tie order is A, T, C, G and only uppercase bases are recognized.
    let bases = *b"ATCG";
    let mut counts = [0usize; 4];
    let mut position = 0usize;
    while position < seq.len() {
        let base = seq[seq.len() - position - 1];
        if let Some(index) = bases.iter().position(|candidate| *candidate == base) {
            counts[index] += 1;
        } else if base == b'N' {
            for count in &mut counts {
                *count += 1;
            }
        }
        let compared = position + 1;
        let allowed = config
            .max_mismatches
            .min(compared / config.mismatch_per_bases.get());
        let should_break = counts
            .iter()
            .all(|count| compared.saturating_sub(*count) > allowed);
        if should_break
            && (position >= config.mismatch_per_bases.get() || position + 1 >= min_len - 1)
        {
            break;
        }
        position += 1;
    }
    if position + 1 < min_len {
        return None;
    }
    let dominant = counts
        .iter()
        .enumerate()
        .max_by_key(|(index, count)| (**count, std::cmp::Reverse(*index)))
        .map(|(index, _)| bases[index])
        .expect("four poly-X counters are always present");
    position = position.min(seq.len() - 1);
    while position > 0 && seq[seq.len() - position - 1] != dominant {
        position -= 1;
    }
    if seq[seq.len() - position - 1] != dominant {
        return None;
    }
    Some(seq.len() - position - 1)
}

fn metric_count(value: usize) -> u64 {
    u64::try_from(value).expect("usize fits in u64 on every supported 64-bit target")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(seq: &[u8]) -> OwnedRecord {
        OwnedRecord {
            id: b"read".to_vec(),
            seq: seq.to_vec(),
            qual: Some(vec![b'I'; seq.len()]),
        }
    }

    #[test]
    fn fixed_trim_keeps_quality_in_lockstep() {
        let mut value = record(b"NNACGTAA");
        let metrics = trim_record(&mut value, FixedTrim { front: 2, tail: 2 }, None, None).unwrap();
        assert_eq!(value.seq, b"ACGT");
        assert_eq!(value.qual.as_deref(), Some(b"IIII".as_slice()));
        assert_eq!(metrics.fixed_bases, 4);
    }

    #[test]
    fn poly_g_exact_cut_includes_allowed_internal_mismatch() {
        let mut value = record(b"ACGTACGTGGGGGAGGGG");
        let config = PolyTailConfig {
            min_len: 9,
            max_mismatches: 1,
            mismatch_per_bases: NonZeroUsize::new(8).expect("nonzero"),
        };
        let metrics = trim_record(&mut value, FixedTrim::default(), Some(config), None).unwrap();
        assert_eq!(value.seq, b"ACGTACGT");
        assert_eq!(metrics.poly_g_bases, 10);
    }

    #[test]
    fn lowercase_poly_g_is_not_fastp_poly_g() {
        let mut value = record(b"ACGTACGTgggggggggg");
        let metrics = trim_record(
            &mut value,
            FixedTrim::default(),
            Some(PolyTailConfig::default()),
            None,
        )
        .unwrap();
        assert!(!metrics.poly_g_matched);
        assert_eq!(value.seq, b"ACGTACGTgggggggggg");
    }

    #[test]
    fn zero_poly_x_minimum_is_rejected() {
        let mut value = record(b"AAAAAAAAAAAA");
        let config = PolyTailConfig {
            min_len: 0,
            ..PolyTailConfig::default()
        };
        assert!(trim_record(&mut value, FixedTrim::default(), None, Some(config)).is_err());
    }

    #[test]
    fn malformed_owned_record_is_rejected() {
        let mut missing = OwnedRecord {
            id: b"read".to_vec(),
            seq: b"ACGT".to_vec(),
            qual: None,
        };
        assert!(trim_record(&mut missing, FixedTrim::default(), None, None).is_err());

        let mut mismatched = record(b"ACGT");
        mismatched.qual.as_mut().unwrap().pop();
        assert!(trim_record(&mut mismatched, FixedTrim::default(), None, None).is_err());
    }
}
