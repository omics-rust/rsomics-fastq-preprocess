use std::fs::File;
use std::io::Stdout;
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use rsomics_seqio::{Compression, Format, OwnedRecord, Writer};
use tempfile::NamedTempFile;

use crate::parallel_gzip::ParallelGzipWriter;

pub(crate) enum OutputSink {
    Stdout(Option<Writer<Stdout>>),
    File(TransactionOutput),
}

impl OutputSink {
    pub(crate) fn create(path: &str, gzip_level: u32) -> Result<Self> {
        if path == "-" {
            Ok(Self::Stdout(Some(Writer::new(
                std::io::stdout(),
                Format::Fastq,
            ))))
        } else {
            TransactionOutput::create(Path::new(path), gzip_level).map(Self::File)
        }
    }

    pub(crate) fn write(&mut self, record: &OwnedRecord) -> Result<()> {
        match self {
            Self::Stdout(writer) => writer
                .as_mut()
                .expect("stdout writer is present until finish")
                .write_owned(record),
            Self::File(output) => output.write(record),
        }
    }

    pub(crate) fn finish(mut self) -> Result<()> {
        match &mut self {
            Self::Stdout(writer) => writer
                .take()
                .expect("stdout writer is present until finish")
                .finish(),
            Self::File(output) => {
                output.prepare()?;
                output.commit()
            }
        }
    }
}

pub(crate) struct TransactionOutput {
    final_path: PathBuf,
    temporary: Option<NamedTempFile>,
    writer: Option<TransactionWriter>,
}

enum TransactionWriter {
    Plain(Writer<File>),
    Gzip(Writer<ParallelGzipWriter<File>>),
}

impl TransactionWriter {
    fn new(file: File, compression: Compression) -> Result<Self> {
        match compression {
            Compression::Plain => Ok(Self::Plain(Writer::new(file, Format::Fastq))),
            Compression::Gzip { level } => {
                let gzip = ParallelGzipWriter::new(file, level).map_err(RsomicsError::Io)?;
                Ok(Self::Gzip(Writer::new(gzip, Format::Fastq)))
            }
        }
    }

    fn write(&mut self, record: &OwnedRecord) -> Result<()> {
        match self {
            Self::Plain(writer) => writer.write_owned(record),
            Self::Gzip(writer) => writer.write_owned(record),
        }
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(writer) => writer.finish(),
            Self::Gzip(writer) => writer
                .finish_into_inner()?
                .finish()
                .map(drop)
                .map_err(RsomicsError::Io),
        }
    }
}

impl TransactionOutput {
    pub(crate) fn create(final_path: &Path, gzip_level: u32) -> Result<Self> {
        reject_occupied(final_path)?;
        let parent = final_path.parent().unwrap_or_else(|| Path::new("."));
        let temporary = NamedTempFile::new_in(parent).rs_with_context(|| {
            format!("creating temporary output beside {}", final_path.display())
        })?;
        let file = temporary
            .reopen()
            .rs_with_context(|| format!("opening temporary output for {}", final_path.display()))?;
        let compression = compression_for(final_path, gzip_level);
        let writer = TransactionWriter::new(file, compression)?;
        Ok(Self {
            final_path: final_path.to_owned(),
            temporary: Some(temporary),
            writer: Some(writer),
        })
    }

    pub(crate) fn write(&mut self, record: &OwnedRecord) -> Result<()> {
        self.writer
            .as_mut()
            .expect("transaction writer is present until prepare")
            .write(record)
    }

    pub(crate) fn prepare(&mut self) -> Result<()> {
        self.writer
            .take()
            .expect("transaction writer is prepared exactly once")
            .finish()
            .rs_with_context(|| {
                format!(
                    "finishing temporary output for {}",
                    self.final_path.display()
                )
            })
    }

    pub(crate) fn commit(&mut self) -> Result<()> {
        let temporary = self
            .temporary
            .take()
            .expect("transaction output is committed exactly once");
        temporary
            .persist_noclobber(&self.final_path)
            .map(drop)
            .map_err(|error| {
                RsomicsError::Io(std::io::Error::new(
                    error.error.kind(),
                    format!("committing {}: {}", self.final_path.display(), error.error),
                ))
            })
    }

    pub(crate) fn final_path(&self) -> &Path {
        &self.final_path
    }
}

pub(crate) fn finish_pair(
    mut first: TransactionOutput,
    mut second: TransactionOutput,
) -> Result<()> {
    first.prepare()?;
    second.prepare()?;
    first.commit()?;
    if let Err(error) = second.commit() {
        let rollback = std::fs::remove_file(first.final_path());
        return match rollback {
            Ok(()) => Err(error),
            Err(rollback_error) => Err(RsomicsError::Io(std::io::Error::new(
                rollback_error.kind(),
                format!(
                    "{error}; additionally failed to roll back {}: {rollback_error}",
                    first.final_path().display()
                ),
            ))),
        };
    }
    Ok(())
}

pub(crate) fn validate_output_path(path: &str) -> Result<()> {
    if path != "-" {
        reject_occupied(Path::new(path))?;
    }
    Ok(())
}

fn reject_occupied(path: &Path) -> Result<()> {
    match std::fs::symlink_metadata(path) {
        Ok(_) => Err(RsomicsError::ConfigError(format!(
            "output {} already exists; refusing to overwrite",
            path.display()
        ))),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(RsomicsError::Io(std::io::Error::new(
            error.kind(),
            format!("checking output {}: {error}", path.display()),
        ))),
    }
}

fn compression_for(path: &Path, gzip_level: u32) -> Compression {
    if path.extension().is_some_and(|extension| extension == "gz") {
        Compression::Gzip { level: gzip_level }
    } else {
        Compression::Plain
    }
}
