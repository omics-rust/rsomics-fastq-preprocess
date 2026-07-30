use std::io::BufRead;

use rsomics_common::{Context, Result, RsomicsError};
use rsomics_seqio::{Format, OwnedRecord, PathReader, Reader, Record};

pub(crate) trait RecordSource {
    fn format(&self) -> Format;
    fn read_record(&mut self) -> Result<Option<Record<'_>>>;
}

impl RecordSource for PathReader {
    fn format(&self) -> Format {
        self.format()
    }

    fn read_record(&mut self) -> Result<Option<Record<'_>>> {
        self.read_record()
    }
}

impl<R: BufRead> RecordSource for Reader<R> {
    fn format(&self) -> Format {
        self.format()
    }

    fn read_record(&mut self) -> Result<Option<Record<'_>>> {
        self.read_record()
    }
}

pub(crate) fn ensure_fastq(source: &impl RecordSource, label: &str) -> Result<()> {
    if source.format() != Format::Fastq {
        return Err(RsomicsError::InvalidInput(format!(
            "{label} contains FASTA; FASTQ input is required"
        )));
    }
    Ok(())
}

pub(crate) fn next_owned(
    source: &mut impl RecordSource,
    label: &str,
    completed: u64,
) -> Result<Option<OwnedRecord>> {
    source
        .read_record()
        .rs_with_context(|| format!("reading {label} after record {completed}"))
        .map(|record| record.map(Record::to_owned))
}

pub(crate) fn same_mate_id(left: &[u8], right: &[u8]) -> bool {
    let Some((left_id, left_role)) = mate_identity(left) else {
        return false;
    };
    let Some((right_id, right_role)) = mate_identity(right) else {
        return false;
    };
    left_id == right_id
        && matches!(
            (left_role, right_role),
            (None, None) | (Some(MateRole::R1), Some(MateRole::R2))
        )
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum MateRole {
    R1,
    R2,
}

fn mate_identity(id: &[u8]) -> Option<(&[u8], Option<MateRole>)> {
    let mut fields = id.splitn(2, |byte| byte.is_ascii_whitespace());
    let token = fields.next().unwrap_or(id);
    let comment = fields.next().and_then(|value| {
        value
            .split(|byte| byte.is_ascii_whitespace())
            .find(|field| !field.is_empty())
    });
    let (normalized, suffix_role) = if let Some(value) = token.strip_suffix(b"/1") {
        (value, Some(MateRole::R1))
    } else if let Some(value) = token.strip_suffix(b"/2") {
        (value, Some(MateRole::R2))
    } else {
        (token, None)
    };
    let comment_role = match comment {
        Some(value) if value.starts_with(b"1:") => Some(MateRole::R1),
        Some(value) if value.starts_with(b"2:") => Some(MateRole::R2),
        _ => None,
    };
    if suffix_role.is_some() && comment_role.is_some() && suffix_role != comment_role {
        return None;
    }
    Some((normalized, suffix_role.or(comment_role)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mate_ids_allow_standard_suffixes_and_comments() {
        assert!(same_mate_id(b"read/1", b"read/2"));
        assert!(same_mate_id(b"read 1:N:0:1", b"read 2:N:0:1"));
        assert!(same_mate_id(b"read comment", b"read other-comment"));
        assert!(!same_mate_id(b"read-a/1", b"read-b/2"));
        assert!(!same_mate_id(b"read/2", b"read/1"));
        assert!(!same_mate_id(b"read/1", b"read/1"));
        assert!(!same_mate_id(b"read 2:N:0:1", b"read 1:N:0:1"));
        assert!(!same_mate_id(b"read/1 2:N:0:1", b"read/2 2:N:0:1"));
    }
}
