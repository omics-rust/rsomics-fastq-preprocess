use std::path::Path;

use rayon::prelude::*;
use rsomics_common::{Context, Result, RsomicsError};
use rsomics_seqio::{OwnedRecord, open_path, open_reader};
use serde::Serialize;

use crate::input::{RecordSource, ensure_fastq, next_owned, same_mate_id};
use crate::output::{OutputSink, TransactionOutput, finish_pair, validate_output_path};
use crate::transform::filter::{FilterConfig, FilterOutcome};
use crate::transform::trim::{RecordTrimMetrics, TrimConfig, trim_record};

const CHUNK_RECORDS: usize = 4_096;

/// Selected public operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Operation {
    /// Compose trimming and filtering in one pass.
    Run,
    /// Apply trimming and post-trim length filtering.
    Trim,
    /// Apply whole-read filters without trimming.
    Filter,
}

/// Input layout inferred from CLI arguments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Mode {
    /// Single-end FASTQ.
    #[serde(rename = "SE")]
    SingleEnd,
    /// Paired-end FASTQ.
    #[serde(rename = "PE")]
    PairedEnd,
}

/// User-visible input and output paths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IoSpec {
    /// R1 input or `-` for stdin in single-end mode.
    pub input_r1: String,
    /// Optional R2 input.
    pub input_r2: Option<String>,
    /// R1 output or `-` for stdout in single-end mode.
    pub output_r1: String,
    /// Optional R2 output.
    pub output_r2: Option<String>,
    /// Gzip compression level for `.gz` outputs.
    pub gzip_level: u32,
}

/// Complete pipeline policy.
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Operation name reported in JSON.
    pub operation: Operation,
    /// Optional trimming stage.
    pub trim: Option<TrimConfig>,
    /// Optional filtering stage.
    pub filter: Option<FilterConfig>,
}

/// Aggregate trimming counters.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct TrimMetrics {
    /// Reads with a poly-G match.
    pub poly_g_trimmed_reads: u64,
    /// Bases removed by poly-G trimming.
    pub poly_g_trimmed_bases: u64,
    /// Reads with a dominant poly-X match.
    pub poly_x_trimmed_reads: u64,
    /// Bases removed by poly-X trimming.
    pub poly_x_trimmed_bases: u64,
    /// Bases removed by fixed trimming.
    pub fixed_trimmed_bases: u64,
}

/// Aggregate whole-read filter counters.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct FilterMetrics {
    /// Reads attributed to excessive unqualified-base percentage.
    pub reads_failed_quality: u64,
    /// Reads attributed to low integer mean quality.
    pub reads_failed_average_quality: u64,
    /// Reads attributed to excessive N bases.
    pub reads_failed_n_bases: u64,
    /// Reads shorter than the configured minimum.
    pub reads_failed_too_short: u64,
    /// Reads longer than the configured maximum.
    pub reads_failed_too_long: u64,
    /// Reads attributed to low complexity.
    pub reads_failed_complexity: u64,
}

/// Stable JSON report for all three subcommands.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct PreprocessReport {
    /// Executed operation.
    pub operation: Operation,
    /// Single-end or paired-end layout.
    pub mode: Mode,
    /// R1 input path.
    pub input_r1: String,
    /// R2 input path.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_r2: Option<String>,
    /// R1 output path.
    pub output_r1: String,
    /// R2 output path.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output_r2: Option<String>,
    /// Input records; paired reads count as two.
    pub reads_in: u64,
    /// Emitted records; paired reads count as two.
    pub reads_out: u64,
    /// Input pairs.
    pub pairs_in: u64,
    /// Emitted pairs.
    pub pairs_out: u64,
    /// Input sequence bases.
    pub bases_in: u64,
    /// Emitted sequence bases.
    pub bases_out: u64,
    /// Trimming counters.
    pub trimming: TrimMetrics,
    /// Filtering counters.
    pub filtering: FilterMetrics,
}

impl PreprocessReport {
    fn new(config: &PipelineConfig, io: &IoSpec, mode: Mode) -> Self {
        Self {
            operation: config.operation,
            mode,
            input_r1: io.input_r1.clone(),
            input_r2: io.input_r2.clone(),
            output_r1: io.output_r1.clone(),
            output_r2: io.output_r2.clone(),
            reads_in: 0,
            reads_out: 0,
            pairs_in: 0,
            pairs_out: 0,
            bases_in: 0,
            bases_out: 0,
            trimming: TrimMetrics::default(),
            filtering: FilterMetrics::default(),
        }
    }

    fn merge(&mut self, delta: &Self) -> Result<()> {
        macro_rules! add {
            ($left:expr, $right:expr, $name:literal) => {
                $left = $left.checked_add($right).ok_or_else(|| {
                    RsomicsError::InvalidInput(concat!($name, " exceeds u64 capacity").into())
                })?
            };
        }
        add!(self.reads_in, delta.reads_in, "read count");
        add!(self.reads_out, delta.reads_out, "output read count");
        add!(self.pairs_in, delta.pairs_in, "pair count");
        add!(self.pairs_out, delta.pairs_out, "output pair count");
        add!(self.bases_in, delta.bases_in, "input base count");
        add!(self.bases_out, delta.bases_out, "output base count");
        add!(
            self.trimming.poly_g_trimmed_reads,
            delta.trimming.poly_g_trimmed_reads,
            "poly-G-trimmed read count"
        );
        add!(
            self.trimming.poly_g_trimmed_bases,
            delta.trimming.poly_g_trimmed_bases,
            "poly-G-trimmed base count"
        );
        add!(
            self.trimming.poly_x_trimmed_reads,
            delta.trimming.poly_x_trimmed_reads,
            "poly-X-trimmed read count"
        );
        add!(
            self.trimming.poly_x_trimmed_bases,
            delta.trimming.poly_x_trimmed_bases,
            "poly-X-trimmed base count"
        );
        add!(
            self.trimming.fixed_trimmed_bases,
            delta.trimming.fixed_trimmed_bases,
            "fixed-trimmed base count"
        );
        add!(
            self.filtering.reads_failed_quality,
            delta.filtering.reads_failed_quality,
            "quality failure count"
        );
        add!(
            self.filtering.reads_failed_average_quality,
            delta.filtering.reads_failed_average_quality,
            "average-quality failure count"
        );
        add!(
            self.filtering.reads_failed_n_bases,
            delta.filtering.reads_failed_n_bases,
            "N-base failure count"
        );
        add!(
            self.filtering.reads_failed_too_short,
            delta.filtering.reads_failed_too_short,
            "too-short failure count"
        );
        add!(
            self.filtering.reads_failed_too_long,
            delta.filtering.reads_failed_too_long,
            "too-long failure count"
        );
        add!(
            self.filtering.reads_failed_complexity,
            delta.filtering.reads_failed_complexity,
            "complexity failure count"
        );
        Ok(())
    }
}

/// Executes a single-end or paired-end pipeline.
pub fn execute(io: &IoSpec, config: &PipelineConfig) -> Result<PreprocessReport> {
    let mode = validate_io(io)?;
    match mode {
        Mode::SingleEnd => execute_single_end(io, config),
        Mode::PairedEnd => execute_paired_end(io, config),
    }
}

fn validate_io(io: &IoSpec) -> Result<Mode> {
    if io.gzip_level > 9 {
        return Err(RsomicsError::ConfigError(
            "gzip level must be in 0..=9".into(),
        ));
    }
    validate_output_path(&io.output_r1)?;
    match (&io.input_r2, &io.output_r2) {
        (None, None) => Ok(Mode::SingleEnd),
        (Some(input_r2), Some(output_r2)) => {
            if io.input_r1 == "-" || input_r2 == "-" {
                return Err(RsomicsError::ConfigError(
                    "paired-end mode requires two file inputs; interleaved stdin is not in this slice"
                        .into(),
                ));
            }
            if io.output_r1 == "-" || output_r2 == "-" {
                return Err(RsomicsError::ConfigError(
                    "paired-end mode requires two file outputs; interleaved stdout is not in this slice"
                        .into(),
                ));
            }
            if input_paths_alias(Path::new(&io.input_r1), Path::new(input_r2))? {
                return Err(RsomicsError::ConfigError(
                    "paired-end inputs must be distinct files".into(),
                ));
            }
            validate_output_path(output_r2)?;
            if normalized_output(&io.output_r1)? == normalized_output(output_r2)? {
                return Err(RsomicsError::ConfigError(
                    "paired-end outputs must be distinct paths".into(),
                ));
            }
            Ok(Mode::PairedEnd)
        }
        _ => Err(RsomicsError::ConfigError(
            "--in2 and --out2 must be supplied together for paired-end mode".into(),
        )),
    }
}

fn input_paths_alias(left: &Path, right: &Path) -> Result<bool> {
    if left == right {
        return Ok(true);
    }
    match same_file::is_same_file(left, right) {
        Ok(value) => Ok(value),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(error) => Err(RsomicsError::Io(std::io::Error::new(
            error.kind(),
            format!(
                "comparing paired inputs {} and {}: {error}",
                left.display(),
                right.display()
            ),
        ))),
    }
}

fn normalized_output(path: &str) -> Result<std::path::PathBuf> {
    let path = Path::new(path);
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent = parent
        .canonicalize()
        .rs_with_context(|| format!("resolving output parent {}", parent.display()))?;
    let name = path.file_name().ok_or_else(|| {
        RsomicsError::ConfigError(format!("output {} has no file name", path.display()))
    })?;
    Ok(parent.join(name))
}

fn execute_single_end(io: &IoSpec, config: &PipelineConfig) -> Result<PreprocessReport> {
    let mut output = OutputSink::create(&io.output_r1, io.gzip_level)?;
    let mut report = PreprocessReport::new(config, io, Mode::SingleEnd);
    if io.input_r1 == "-" {
        let stdin = std::io::stdin();
        let mut reader = open_reader(stdin.lock()).rs_context("opening stdin")?;
        ensure_fastq(&reader, "stdin")?;
        process_single_reader(&mut reader, &mut output, config, &mut report, "stdin")?;
    } else {
        let mut reader = open_path(Path::new(&io.input_r1))
            .rs_with_context(|| format!("opening {}", io.input_r1))?;
        ensure_fastq(&reader, &io.input_r1)?;
        process_single_reader(&mut reader, &mut output, config, &mut report, &io.input_r1)?;
    }
    output.finish()?;
    Ok(report)
}

fn process_single_reader(
    reader: &mut impl RecordSource,
    output: &mut OutputSink,
    config: &PipelineConfig,
    report: &mut PreprocessReport,
    label: &str,
) -> Result<()> {
    let mut completed = 0u64;
    loop {
        let mut chunk = Vec::with_capacity(CHUNK_RECORDS);
        while chunk.len() < CHUNK_RECORDS {
            let Some(record) = next_owned(reader, label, completed)? else {
                break;
            };
            completed = completed.checked_add(1).ok_or_else(|| {
                RsomicsError::InvalidInput("read count exceeds u64 capacity".into())
            })?;
            chunk.push(record);
        }
        if chunk.is_empty() {
            break;
        }
        let processed: Result<Vec<_>> = if rayon::current_num_threads() == 1 {
            chunk
                .into_iter()
                .map(|record| process_single(record, config))
                .collect()
        } else {
            chunk
                .into_par_iter()
                .map(|record| process_single(record, config))
                .collect()
        };
        for item in processed? {
            report.merge(&item.delta)?;
            if let Some(record) = item.output {
                output.write(&record)?;
            }
        }
    }
    Ok(())
}

fn execute_paired_end(io: &IoSpec, config: &PipelineConfig) -> Result<PreprocessReport> {
    let input_r2 = io.input_r2.as_deref().expect("validated paired input");
    let output_r2 = io.output_r2.as_deref().expect("validated paired output");
    let mut first_reader = open_path(Path::new(&io.input_r1))
        .rs_with_context(|| format!("opening {}", io.input_r1))?;
    let mut second_reader =
        open_path(Path::new(input_r2)).rs_with_context(|| format!("opening {input_r2}"))?;
    ensure_fastq(&first_reader, &io.input_r1)?;
    ensure_fastq(&second_reader, input_r2)?;
    let mut first_output = TransactionOutput::create(Path::new(&io.output_r1), io.gzip_level)?;
    let mut second_output = TransactionOutput::create(Path::new(output_r2), io.gzip_level)?;
    let mut report = PreprocessReport::new(config, io, Mode::PairedEnd);
    let mut completed = 0u64;

    loop {
        let mut chunk = Vec::with_capacity(CHUNK_RECORDS);
        while chunk.len() < CHUNK_RECORDS {
            let left = next_owned(&mut first_reader, &io.input_r1, completed)?;
            let right = next_owned(&mut second_reader, input_r2, completed)?;
            match (left, right) {
                (Some(left), Some(right)) => {
                    if !same_mate_id(&left.id, &right.id) {
                        return Err(RsomicsError::InvalidInput(format!(
                            "paired FASTQ identifiers diverge at pair {}: {:?} vs {:?}",
                            completed + 1,
                            String::from_utf8_lossy(&left.id),
                            String::from_utf8_lossy(&right.id)
                        )));
                    }
                    completed = completed.checked_add(1).ok_or_else(|| {
                        RsomicsError::InvalidInput("pair count exceeds u64 capacity".into())
                    })?;
                    chunk.push((left, right));
                }
                (None, None) => break,
                _ => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "paired FASTQ record counts diverge after {completed} complete pairs"
                    )));
                }
            }
        }
        if chunk.is_empty() {
            break;
        }
        let processed: Result<Vec<_>> = if rayon::current_num_threads() == 1 {
            chunk
                .into_iter()
                .map(|pair| process_pair(pair, config))
                .collect()
        } else {
            chunk
                .into_par_iter()
                .map(|pair| process_pair(pair, config))
                .collect()
        };
        for item in processed? {
            report.merge(&item.delta)?;
            if let Some((left, right)) = item.output {
                first_output.write(&left)?;
                second_output.write(&right)?;
            }
        }
    }
    finish_pair(first_output, second_output)?;
    Ok(report)
}

struct ProcessedSingle {
    delta: PreprocessReport,
    output: Option<OwnedRecord>,
}

fn process_single(mut record: OwnedRecord, config: &PipelineConfig) -> Result<ProcessedSingle> {
    let input_bases = metric_count(record.seq.len());
    let mut delta = blank_delta(config, Mode::SingleEnd);
    delta.reads_in = 1;
    delta.bases_in = input_bases;

    if let Some(trim) = config.trim.as_ref() {
        let metrics = trim_record(&mut record, trim.fixed_r1, trim.poly_g, trim.poly_x);
        add_trim_record(&mut delta.trimming, metrics);
    }

    let outcome = check_record(&record, config.filter)?;
    if outcome == FilterOutcome::Pass {
        normalize_output_quality(&mut record, config.filter);
        delta.reads_out = 1;
        delta.bases_out = metric_count(record.seq.len());
        Ok(ProcessedSingle {
            delta,
            output: Some(record),
        })
    } else {
        add_filter_outcome(&mut delta.filtering, outcome);
        Ok(ProcessedSingle {
            delta,
            output: None,
        })
    }
}

struct ProcessedPair {
    delta: PreprocessReport,
    output: Option<(OwnedRecord, OwnedRecord)>,
}

fn process_pair(
    (mut left, mut right): (OwnedRecord, OwnedRecord),
    config: &PipelineConfig,
) -> Result<ProcessedPair> {
    let mut delta = blank_delta(config, Mode::PairedEnd);
    delta.reads_in = 2;
    delta.pairs_in = 1;
    delta.bases_in = metric_count(left.seq.len())
        .checked_add(metric_count(right.seq.len()))
        .expect("two live record buffers fit in the supported address space");

    if let Some(trim) = config.trim.as_ref() {
        let left_metrics = trim_record(&mut left, trim.fixed_r1, trim.poly_g, trim.poly_x);
        let right_metrics = trim_record(&mut right, trim.fixed_r2, trim.poly_g, trim.poly_x);
        add_trim_record(&mut delta.trimming, left_metrics);
        add_trim_record(&mut delta.trimming, right_metrics);
    }

    let left_outcome = check_record(&left, config.filter)?;
    let right_outcome = check_record(&right, config.filter)?;
    if left_outcome == FilterOutcome::Pass && right_outcome == FilterOutcome::Pass {
        normalize_output_quality(&mut left, config.filter);
        normalize_output_quality(&mut right, config.filter);
        delta.reads_out = 2;
        delta.pairs_out = 1;
        delta.bases_out = metric_count(left.seq.len())
            .checked_add(metric_count(right.seq.len()))
            .expect("two live record buffers fit in the supported address space");
        Ok(ProcessedPair {
            delta,
            output: Some((left, right)),
        })
    } else {
        add_filter_outcome(&mut delta.filtering, left_outcome);
        add_filter_outcome(&mut delta.filtering, right_outcome);
        Ok(ProcessedPair {
            delta,
            output: None,
        })
    }
}

fn normalize_output_quality(record: &mut OwnedRecord, filter: Option<FilterConfig>) {
    if filter.is_some_and(|config| config.phred == crate::PhredEncoding::Phred64) {
        for quality in record
            .qual
            .as_mut()
            .expect("validated FASTQ records always carry quality")
        {
            *quality = quality
                .checked_sub(31)
                .expect("Phred+64 validation precedes output normalization");
        }
    }
}

fn check_record(record: &OwnedRecord, filter: Option<FilterConfig>) -> Result<FilterOutcome> {
    let Some(filter) = filter else {
        return Ok(FilterOutcome::Pass);
    };
    filter.check_seqio_record(
        &record.seq,
        record
            .qual
            .as_deref()
            .expect("validated FASTQ records always carry quality"),
    )
}

fn blank_delta(config: &PipelineConfig, mode: Mode) -> PreprocessReport {
    PreprocessReport {
        operation: config.operation,
        mode,
        input_r1: String::new(),
        input_r2: None,
        output_r1: String::new(),
        output_r2: None,
        reads_in: 0,
        reads_out: 0,
        pairs_in: 0,
        pairs_out: 0,
        bases_in: 0,
        bases_out: 0,
        trimming: TrimMetrics::default(),
        filtering: FilterMetrics::default(),
    }
}

fn add_trim_record(total: &mut TrimMetrics, record: RecordTrimMetrics) {
    total.fixed_trimmed_bases += record.fixed_bases;
    total.poly_g_trimmed_reads += u64::from(record.poly_g_matched);
    total.poly_g_trimmed_bases += record.poly_g_bases;
    total.poly_x_trimmed_reads += u64::from(record.poly_x_matched);
    total.poly_x_trimmed_bases += record.poly_x_bases;
}

fn add_filter_outcome(total: &mut FilterMetrics, outcome: FilterOutcome) {
    match outcome {
        FilterOutcome::Pass => {}
        FilterOutcome::LowQuality => total.reads_failed_quality += 1,
        FilterOutcome::LowAverageQuality => total.reads_failed_average_quality += 1,
        FilterOutcome::TooManyN => total.reads_failed_n_bases += 1,
        FilterOutcome::TooShort => total.reads_failed_too_short += 1,
        FilterOutcome::TooLong => total.reads_failed_too_long += 1,
        FilterOutcome::LowComplexity => total.reads_failed_complexity += 1,
    }
}

fn metric_count(value: usize) -> u64 {
    u64::try_from(value).expect("usize fits in u64 on every supported 64-bit target")
}
