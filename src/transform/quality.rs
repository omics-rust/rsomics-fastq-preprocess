use std::ops::Range;

use rsomics_common::{Result, RsomicsError};

use super::filter::PhredEncoding;
use super::trim::FixedTrim;

/// A fastp-compatible rolling quality window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QualityWindow {
    /// Window length in 1..=1000 bases.
    pub window_size: usize,
    /// Minimum mean Phred score in 1..=30.
    pub mean_quality: u8,
}

/// Optional quality cuts, applied front then right or tail.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct QualityCut {
    /// Remove low-quality leading windows.
    pub front: Option<QualityWindow>,
    /// Remove low-quality trailing windows; ignored when right is enabled.
    pub tail: Option<QualityWindow>,
    /// Stop at the first low-quality window while retaining its good prefix.
    pub right: Option<QualityWindow>,
}

impl QualityCut {
    pub(crate) fn enabled(self) -> bool {
        self.front.is_some() || self.tail.is_some() || self.right.is_some()
    }

    pub(crate) fn validate(self) -> Result<()> {
        for window in [self.front, self.tail, self.right].into_iter().flatten() {
            if !(1..=1000).contains(&window.window_size) || !(1..=30).contains(&window.mean_quality)
            {
                return Err(RsomicsError::ConfigError(
                    "quality-cut windows require size 1..=1000 and mean quality 1..=30".into(),
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn range(
        self,
        seq: &[u8],
        quality: &[u8],
        fixed: FixedTrim,
        phred: PhredEncoding,
    ) -> Option<Range<usize>> {
        let length = seq.len();
        let mut start = fixed.front;
        let mut end = length.saturating_sub(fixed.tail);
        if let Some(window) = self.front {
            if start >= end || end - start <= window.window_size {
                return None;
            }
            start = first_window(quality, start..end, window, phred, false).0;
            if start > 0 {
                start += window.window_size - 1;
            }
            while start < length && seq[start] == b'N' {
                start += 1;
            }
        }
        if let Some(window) = self.right {
            if start >= end || end - start <= window.window_size {
                return None;
            }
            let (mut index, found) = first_window(quality, start..end, window, phred, true);
            if found {
                let threshold = phred.offset() + window.mean_quality;
                while index < length - 1 && quality[index] >= threshold {
                    index += 1;
                }
                end = index;
            }
        } else if let Some(window) = self.tail {
            let width = window.window_size;
            if start >= end || end - start <= width {
                return None;
            }
            let mut total: usize = quality[end - (width - 1)..end]
                .iter()
                .map(|q| usize::from(*q))
                .sum();
            let threshold = width * usize::from(phred.offset() + window.mean_quality);
            let mut last = start + width - 1;
            for index in (start + width..end).rev() {
                total += usize::from(quality[index - width + 1]);
                if index < end - 1 {
                    total -= usize::from(quality[index + 1]);
                }
                if total >= threshold {
                    last = index;
                    break;
                }
            }
            if last < length - 1 {
                last -= width - 1;
            }
            while seq[last] == b'N' {
                last = last.checked_sub(1)?;
            }
            end = last + 1;
        }
        (start < end && start < length.saturating_sub(1)).then_some(start..end)
    }
}

fn first_window(
    quality: &[u8],
    range: Range<usize>,
    window: QualityWindow,
    phred: PhredEncoding,
    below: bool,
) -> (usize, bool) {
    let width = window.window_size;
    let threshold = width * usize::from(phred.offset() + window.mean_quality);
    let mut total: usize = quality[range.start..range.start + width - 1]
        .iter()
        .map(|q| usize::from(*q))
        .sum();
    for index in range.start..range.end - width {
        total += usize::from(quality[index + width - 1]);
        if index > range.start {
            total -= usize::from(quality[index - 1]);
        }
        if (total < threshold) == below {
            return (index, true);
        }
    }
    (range.end - width, false)
}
