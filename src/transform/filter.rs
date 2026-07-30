use rsomics_common::{Result, RsomicsError};

/// Explicit FASTQ quality-score encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhredEncoding {
    /// Modern Sanger and Illumina 1.8+ encoding.
    Phred33,
    /// Legacy Illumina 1.3-1.7 encoding.
    Phred64,
}

impl PhredEncoding {
    /// Raw ASCII offset.
    #[must_use]
    pub const fn offset(self) -> u8 {
        match self {
            Self::Phred33 => 33,
            Self::Phred64 => 64,
        }
    }
}

/// Ordered filter outcome matching fastp's first-failure attribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterOutcome {
    /// The read passes every enabled filter.
    Pass,
    /// Too many bases are below the quality threshold.
    LowQuality,
    /// Mean quality is below the configured threshold.
    LowAverageQuality,
    /// Too many uppercase `N` bases are present.
    TooManyN,
    /// The read is shorter than the accepted minimum.
    TooShort,
    /// The read is longer than the accepted maximum.
    TooLong,
    /// Adjacent-base change fraction is below the configured threshold.
    LowComplexity,
}

/// Whole-read filtering policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FilterConfig {
    /// Explicit quality-score encoding.
    pub phred: PhredEncoding,
    /// Enable per-base quality, average-quality, and N-count filters.
    pub quality_enabled: bool,
    /// Minimum Phred score for a qualified base.
    pub qualified_quality_phred: u8,
    /// Maximum percentage of unqualified bases; equality passes.
    pub unqualified_percent_limit: u8,
    /// Minimum integer mean Phred score; zero disables this check.
    pub average_quality: u8,
    /// Maximum number of uppercase `N` bases; equality passes.
    pub n_base_limit: usize,
    /// Enable minimum and maximum length filtering.
    pub length_enabled: bool,
    /// Minimum accepted read length.
    pub length_required: usize,
    /// Maximum accepted read length; zero means unlimited.
    pub length_limit: usize,
    /// Optional minimum adjacent-base change percentage.
    pub complexity_threshold: Option<u8>,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            phred: PhredEncoding::Phred33,
            quality_enabled: true,
            qualified_quality_phred: 15,
            unqualified_percent_limit: 40,
            average_quality: 0,
            n_base_limit: 5,
            length_enabled: true,
            length_required: 15,
            length_limit: 0,
            complexity_threshold: None,
        }
    }
}

impl FilterConfig {
    /// Checks one FASTQ record using fastp v1.3.6 filter precedence.
    pub fn check(self, seq: &[u8], qual: &[u8]) -> Result<FilterOutcome> {
        if seq.is_empty()
            && (self.quality_enabled
                || (self.length_enabled && self.length_required > 0)
                || self.complexity_threshold.is_some())
        {
            return Ok(FilterOutcome::TooShort);
        }
        let offset = self.phred.offset();
        if let Some(invalid) = qual.iter().copied().find(|quality| *quality < offset) {
            return Err(RsomicsError::InvalidInput(format!(
                "quality byte {invalid} is below the explicit Phred+{offset} offset"
            )));
        }

        if self.quality_enabled {
            let qualified_threshold = u16::from(self.qualified_quality_phred) + u16::from(offset);
            let mut low_quality = 0usize;
            let mut n_count = 0usize;
            let mut total_quality = 0u64;
            for (&base, &quality) in seq.iter().zip(qual) {
                if u16::from(quality) < qualified_threshold {
                    low_quality += 1;
                }
                if base == b'N' {
                    n_count += 1;
                }
                total_quality += u64::from(quality - offset);
            }
            if metric_count(low_quality) * 100
                > u64::from(self.unqualified_percent_limit) * metric_count(seq.len())
            {
                return Ok(FilterOutcome::LowQuality);
            }
            if self.average_quality > 0
                && total_quality / metric_count(seq.len()) < u64::from(self.average_quality)
            {
                return Ok(FilterOutcome::LowAverageQuality);
            }
            if n_count > self.n_base_limit {
                return Ok(FilterOutcome::TooManyN);
            }
        }

        if self.length_enabled {
            if seq.len() < self.length_required {
                return Ok(FilterOutcome::TooShort);
            }
            if self.length_limit > 0 && seq.len() > self.length_limit {
                return Ok(FilterOutcome::TooLong);
            }
        }

        if let Some(threshold) = self.complexity_threshold
            && !passes_complexity(seq, threshold)
        {
            return Ok(FilterOutcome::LowComplexity);
        }

        Ok(FilterOutcome::Pass)
    }
}

/// Returns whether the adjacent-base change percentage meets the threshold.
///
/// fastp rejects reads of zero or one base whenever this filter is enabled.
#[must_use]
pub fn passes_complexity(seq: &[u8], threshold_percent: u8) -> bool {
    if seq.len() <= 1 {
        return false;
    }
    let changes = seq.windows(2).filter(|pair| pair[0] != pair[1]).count();
    metric_count(changes) * 100 >= u64::from(threshold_percent) * metric_count(seq.len() - 1)
}

fn metric_count(value: usize) -> u64 {
    u64::try_from(value).expect("usize fits in u64 on every supported 64-bit target")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fastp_quality_boundaries_are_strict() {
        let config = FilterConfig {
            length_enabled: false,
            ..FilterConfig::default()
        };
        assert_eq!(
            config
                .check(b"ACGTACGTACGTACGTACGT", b"!!!!!!!!IIIIIIIIIIII")
                .unwrap(),
            FilterOutcome::Pass
        );
        assert_eq!(
            config
                .check(b"ACGTACGTACGTACGTACGT", b"!!!!!!!!!IIIIIIIIIII")
                .unwrap(),
            FilterOutcome::LowQuality
        );
    }

    #[test]
    fn average_quality_precedes_n_count() {
        let config = FilterConfig {
            average_quality: 40,
            length_enabled: false,
            ..FilterConfig::default()
        };
        assert_eq!(
            config.check(b"NNNNNNACGT", b"!!!!!!!!!!").unwrap(),
            FilterOutcome::LowQuality
        );
        let config = FilterConfig {
            qualified_quality_phred: 0,
            average_quality: 40,
            length_enabled: false,
            ..FilterConfig::default()
        };
        assert_eq!(
            config.check(b"NNNNNNACGT", b"!!!!!!!!!!").unwrap(),
            FilterOutcome::LowAverageQuality
        );
    }

    #[test]
    fn complexity_rejects_short_reads_and_uses_case_sensitive_changes() {
        assert!(!passes_complexity(b"", 30));
        assert!(!passes_complexity(b"A", 30));
        assert!(!passes_complexity(b"AAAAAAAAAA", 30));
        assert!(passes_complexity(b"AaAaAaAaAa", 30));
        assert!(passes_complexity(b"ACACACACAC", 30));
    }

    #[test]
    fn explicit_phred64_rejects_phred33_bytes() {
        let config = FilterConfig {
            phred: PhredEncoding::Phred64,
            ..FilterConfig::default()
        };
        assert!(config.check(b"ACGT", b"!!!!").is_err());
    }

    #[test]
    fn exact_fastp_n_filter_counts_only_uppercase_n() {
        let config = FilterConfig {
            qualified_quality_phred: 0,
            n_base_limit: 0,
            length_enabled: false,
            ..FilterConfig::default()
        };
        assert_eq!(config.check(b"nnnn", b"IIII").unwrap(), FilterOutcome::Pass);
        assert_eq!(
            config.check(b"Nnnn", b"IIII").unwrap(),
            FilterOutcome::TooManyN
        );
    }

    #[test]
    fn phred64_accepts_inclusive_ascii_boundaries() {
        let config = FilterConfig {
            phred: PhredEncoding::Phred64,
            quality_enabled: false,
            length_enabled: false,
            ..FilterConfig::default()
        };
        assert_eq!(config.check(b"AC", b"@~").unwrap(), FilterOutcome::Pass);
        assert!(config.check(b"A", b"?").is_err());
    }
}
