use std::num::NonZeroUsize;

use clap::{Args, Parser, Subcommand};
use rayon::{ThreadPool, ThreadPoolBuilder};
use rsomics_common::{OutputArgs, Result, RsomicsError, ToolMeta};

use rsomics_fastq_preprocess::{
    FilterConfig, FixedTrim, IoSpec, PhredEncoding, PipelineConfig, PolyTailConfig,
    PreprocessReport, TrimConfig, execute,
};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-fastq-preprocess",
    version,
    about = "Single-pass FASTQ trimming and filtering",
    subcommand_required = true,
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,

    #[command(flatten)]
    pub output: OutputArgs,

    #[command(flatten)]
    pub threads: ThreadArgs,
}

#[derive(Debug, Clone, Args)]
#[command(next_help_heading = "Global options")]
pub struct ThreadArgs {
    /// Number of worker threads.
    #[arg(short = 't', long, global = true)]
    threads: Option<NonZeroUsize>,
}

impl ThreadArgs {
    pub fn build(&self) -> Result<ThreadPool> {
        let mut builder = ThreadPoolBuilder::new();
        if let Some(threads) = self.threads {
            builder = builder.num_threads(threads.get());
        }
        builder.build().map_err(|error| {
            RsomicsError::ConfigError(format!("creating worker thread pool failed: {error}"))
        })
    }
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Trim and filter in one FASTQ traversal.
    Run(RunArgs),
    /// Apply fixed and poly-tail trimming.
    Trim(TrimArgs),
    /// Apply whole-read quality, N, length, and optional complexity filters.
    Filter(FilterOnlyArgs),
}

#[derive(Debug, Args)]
struct RunArgs {
    #[command(flatten)]
    io: IoArgs,
    #[command(flatten)]
    trim: TrimOptions,
    #[command(flatten)]
    filter: FilterOptions,
}

#[derive(Debug, Args)]
struct TrimArgs {
    #[command(flatten)]
    io: IoArgs,
    #[command(flatten)]
    trim: TrimOptions,
    #[command(flatten)]
    output: TrimOutputOptions,
}

#[derive(Debug, Args)]
struct FilterOnlyArgs {
    #[command(flatten)]
    io: IoArgs,
    #[command(flatten)]
    filter: FilterOptions,
}

#[derive(Debug, Clone, Args)]
#[command(next_help_heading = "Input/output")]
struct IoArgs {
    /// R1 FASTQ input; `-` reads stdin in single-end mode.
    #[arg(short = 'i', long = "in1", default_value = "-")]
    input_r1: String,

    /// R2 FASTQ input; requires `--out2` and file-based paired mode.
    #[arg(short = 'I', long = "in2")]
    input_r2: Option<String>,

    /// R1 FASTQ output; `-` writes stdout in single-end mode.
    #[arg(short = 'o', long = "out1", default_value = "-")]
    output_r1: String,

    /// R2 FASTQ output; requires `--in2`.
    #[arg(short = 'O', long = "out2")]
    output_r2: Option<String>,

    /// Gzip level for output paths ending in `.gz`.
    #[arg(long, default_value_t = 4, value_parser = clap::value_parser!(u32).range(0..=9))]
    gzip_level: u32,
}

impl IoArgs {
    fn build(&self) -> IoSpec {
        IoSpec {
            input_r1: self.input_r1.clone(),
            input_r2: self.input_r2.clone(),
            output_r1: self.output_r1.clone(),
            output_r2: self.output_r2.clone(),
            gzip_level: self.gzip_level,
        }
    }
}

#[derive(Debug, Clone, Args)]
#[command(next_help_heading = "Trimming")]
struct TrimOptions {
    /// Trim this many bases from the R1 5-prime end.
    #[arg(short = 'f', long, default_value_t = 0)]
    trim_front1: usize,

    /// Trim this many bases from the R1 3-prime end.
    #[arg(long, default_value_t = 0)]
    trim_tail1: usize,

    /// Trim this many bases from the R2 5-prime end.
    #[arg(short = 'F', long)]
    trim_front2: Option<usize>,

    /// Trim this many bases from the R2 3-prime end.
    #[arg(long)]
    trim_tail2: Option<usize>,

    /// Trim poly-G tails.
    #[arg(short = 'g', long)]
    trim_poly_g: bool,

    /// Trim dominant poly-A/C/G/T tails.
    #[arg(short = 'x', long)]
    trim_poly_x: bool,

    /// Minimum poly-G tail span.
    #[arg(long, default_value_t = 10)]
    poly_g_min_len: usize,

    /// Minimum poly-X tail span.
    #[arg(long, default_value_t = 10)]
    poly_x_min_len: usize,

    /// Absolute mismatch cap in poly-tail scans.
    #[arg(long, default_value_t = 5)]
    poly_max_mismatches: usize,

    /// Allow one poly-tail mismatch per this many scanned bases.
    #[arg(long, default_value_t = 8)]
    poly_mismatch_per_bases: usize,
}

impl TrimOptions {
    fn build(&self) -> Result<TrimConfig> {
        if self.trim_poly_g && self.poly_g_min_len == 0 {
            return Err(RsomicsError::ConfigError(
                "--poly-g-min-len must be positive".into(),
            ));
        }
        if self.trim_poly_x && self.poly_x_min_len == 0 {
            return Err(RsomicsError::ConfigError(
                "--poly-x-min-len must be positive".into(),
            ));
        }
        let mismatch_per_bases =
            NonZeroUsize::new(self.poly_mismatch_per_bases).ok_or_else(|| {
                RsomicsError::ConfigError("--poly-mismatch-per-bases must be positive".into())
            })?;
        let poly = |min_len| PolyTailConfig {
            min_len,
            max_mismatches: self.poly_max_mismatches,
            mismatch_per_bases,
        };
        Ok(TrimConfig {
            fixed_r1: FixedTrim {
                front: self.trim_front1,
                tail: self.trim_tail1,
            },
            fixed_r2: FixedTrim {
                front: self.trim_front2.unwrap_or(self.trim_front1),
                tail: self.trim_tail2.unwrap_or(self.trim_tail1),
            },
            poly_g: self.trim_poly_g.then(|| poly(self.poly_g_min_len)),
            poly_x: self.trim_poly_x.then(|| poly(self.poly_x_min_len)),
        })
    }
}

#[derive(Debug, Clone, Args)]
#[command(next_help_heading = "Output filtering")]
struct TrimOutputOptions {
    /// Quality encoding; Phred+64 input is always emitted as Phred+33.
    #[arg(long, default_value_t = 33, value_parser = clap::value_parser!(u8).range(33..=64))]
    phred_offset: u8,

    /// Enable post-trim minimum-length filtering with this threshold.
    #[arg(short = 'l', long)]
    length_required: Option<usize>,

    /// Enable post-trim maximum-length filtering; zero means unlimited.
    #[arg(long)]
    length_limit: Option<usize>,
}

impl TrimOutputOptions {
    fn build(&self) -> Result<Option<FilterConfig>> {
        let phred = parse_phred(self.phred_offset)?;
        let length_enabled = self.length_required.is_some() || self.length_limit.is_some();
        if phred == PhredEncoding::Phred33 && !length_enabled {
            return Ok(None);
        }
        Ok(Some(FilterConfig {
            phred,
            quality_enabled: false,
            length_enabled,
            length_required: self.length_required.unwrap_or(0),
            length_limit: self.length_limit.unwrap_or(0),
            complexity_threshold: None,
            ..FilterConfig::default()
        }))
    }
}

#[derive(Debug, Clone, Args)]
#[command(next_help_heading = "Length filtering")]
struct LengthOptions {
    /// Minimum accepted read length.
    #[arg(short = 'l', long, default_value_t = 15)]
    length_required: usize,

    /// Maximum accepted read length; zero means unlimited.
    #[arg(long, default_value_t = 0)]
    length_limit: usize,

    /// Disable both minimum and maximum length filtering.
    #[arg(short = 'L', long)]
    disable_length_filtering: bool,
}

impl LengthOptions {
    fn apply(&self, config: &mut FilterConfig) {
        config.length_enabled = !self.disable_length_filtering;
        config.length_required = self.length_required;
        config.length_limit = self.length_limit;
    }
}

#[derive(Debug, Clone, Args)]
#[command(next_help_heading = "Filtering")]
struct FilterOptions {
    /// Quality encoding; accepted values are 33 and 64.
    #[arg(long, default_value_t = 33, value_parser = clap::value_parser!(u8).range(33..=64))]
    phred_offset: u8,

    /// Minimum Phred score for a qualified base.
    #[arg(long, default_value_t = 15, value_parser = clap::value_parser!(u8).range(0..=93))]
    qualified_quality_phred: u8,

    /// Maximum percentage of unqualified bases.
    #[arg(long, default_value_t = 40, value_parser = clap::value_parser!(u8).range(0..=100))]
    unqualified_percent_limit: u8,

    /// Minimum integer mean Phred score; zero disables this filter.
    #[arg(long, default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93))]
    average_quality: u8,

    /// Maximum N bases per read.
    #[arg(long, default_value_t = 5)]
    n_base_limit: usize,

    /// Disable quality and N-base filtering.
    #[arg(short = 'Q', long)]
    disable_quality_filtering: bool,

    /// Enable fastp-style adjacent-base complexity filtering.
    #[arg(short = 'y', long)]
    low_complexity_filter: bool,

    /// Minimum adjacent-base change percentage.
    #[arg(short = 'Y', long, default_value_t = 30, value_parser = clap::value_parser!(u8).range(0..=100))]
    complexity_threshold: u8,

    #[command(flatten)]
    length: LengthOptions,
}

impl FilterOptions {
    fn build(&self) -> Result<FilterConfig> {
        let phred = parse_phred(self.phred_offset)?;
        let mut config = FilterConfig {
            phred,
            quality_enabled: !self.disable_quality_filtering,
            qualified_quality_phred: self.qualified_quality_phred,
            unqualified_percent_limit: self.unqualified_percent_limit,
            average_quality: self.average_quality,
            n_base_limit: self.n_base_limit,
            complexity_threshold: self
                .low_complexity_filter
                .then_some(self.complexity_threshold),
            ..FilterConfig::default()
        };
        self.length.apply(&mut config);
        Ok(config)
    }
}

fn parse_phred(offset: u8) -> Result<PhredEncoding> {
    match offset {
        33 => Ok(PhredEncoding::Phred33),
        64 => Ok(PhredEncoding::Phred64),
        other => Err(RsomicsError::ConfigError(format!(
            "--phred-offset must be 33 or 64, got {other}"
        ))),
    }
}

impl Cli {
    pub fn execute(self) -> Result<PreprocessReport> {
        let json = self.output.json;
        let (io, config) = match self.command {
            Command::Run(args) => (
                args.io.build(),
                PipelineConfig::run(args.trim.build()?, args.filter.build()?),
            ),
            Command::Trim(args) => (
                args.io.build(),
                PipelineConfig::trim(args.trim.build()?, args.output.build()?),
            ),
            Command::Filter(args) => (
                args.io.build(),
                PipelineConfig::filter(args.filter.build()?),
            ),
        };
        if json && io.output_r1 == "-" {
            return Err(RsomicsError::ConfigError(
                "--json requires FASTQ data to be written to file, not stdout".into(),
            ));
        }
        execute(&io, &config)
    }
}

#[cfg(test)]
mod tests {
    use clap::{CommandFactory, Parser as _};

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn thread_limit_uses_a_local_pool() {
        let cli =
            Cli::try_parse_from(["rsomics-fastq-preprocess", "--threads", "1", "filter"]).unwrap();
        let pool = cli.threads.build().unwrap();
        assert_eq!(pool.install(rayon::current_num_threads), 1);
    }

    #[test]
    fn paired_fixed_trim_inherits_only_when_r2_is_absent() {
        let inherited = Cli::try_parse_from([
            "rsomics-fastq-preprocess",
            "trim",
            "--trim-front1",
            "3",
            "--trim-tail1",
            "2",
        ])
        .unwrap();
        let explicit_zero = Cli::try_parse_from([
            "rsomics-fastq-preprocess",
            "trim",
            "--trim-front1",
            "3",
            "--trim-tail1",
            "2",
            "--trim-front2",
            "0",
            "--trim-tail2",
            "0",
        ])
        .unwrap();

        let Command::Trim(inherited) = inherited.command else {
            unreachable!()
        };
        let Command::Trim(explicit_zero) = explicit_zero.command else {
            unreachable!()
        };
        assert_eq!(
            inherited.trim.build().unwrap().fixed_r2,
            FixedTrim { front: 3, tail: 2 }
        );
        assert_eq!(
            explicit_zero.trim.build().unwrap().fixed_r2,
            FixedTrim { front: 0, tail: 0 }
        );
    }

    #[test]
    fn trim_length_filter_is_opt_in() {
        let default =
            Cli::try_parse_from(["rsomics-fastq-preprocess", "trim"]).expect("valid trim CLI");
        let explicit = Cli::try_parse_from([
            "rsomics-fastq-preprocess",
            "trim",
            "--length-required",
            "15",
        ])
        .expect("valid trim CLI");
        let Command::Trim(default) = default.command else {
            unreachable!()
        };
        let Command::Trim(explicit) = explicit.command else {
            unreachable!()
        };
        assert!(default.output.build().unwrap().is_none());
        let filter = explicit.output.build().unwrap().unwrap();
        assert!(filter.length_enabled);
        assert_eq!(filter.length_required, 15);
    }
}
