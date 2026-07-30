//! Ordered parallel gzip compression using concatenated members.

use std::io::{self, BufWriter, Write};

use libdeflater::{CompressionLvl, Compressor};
use rayon::prelude::*;

const CHUNK_BYTES: usize = 256 * 1024;
const MAX_PENDING_CHUNKS: usize = 16;

// Independent gzip members let the existing Rayon pool compress chunks
// concurrently without adding background threads or changing CLI thread limits.
pub(crate) struct ParallelGzipWriter<W: Write> {
    inner: BufWriter<W>,
    buffer: Vec<u8>,
    pending: Vec<Vec<u8>>,
    level: CompressionLvl,
    wrote_member: bool,
}

impl<W: Write> ParallelGzipWriter<W> {
    pub(crate) fn new(inner: W, level: u32) -> io::Result<Self> {
        let level = i32::try_from(level)
            .ok()
            .and_then(|level| CompressionLvl::new(level).ok())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid libdeflate gzip level {level}"),
                )
            })?;
        Ok(Self {
            inner: BufWriter::with_capacity(CHUNK_BYTES * 2, inner),
            buffer: Vec::with_capacity(CHUNK_BYTES),
            pending: Vec::with_capacity(MAX_PENDING_CHUNKS),
            level,
            wrote_member: false,
        })
    }

    pub(crate) fn finish(mut self) -> io::Result<W> {
        if self.buffer.is_empty() && self.pending.is_empty() && !self.wrote_member {
            self.pending.push(Vec::new());
        }
        self.flush()?;
        self.inner.into_inner().map_err(|error| error.into_error())
    }

    fn queue_buffer(&mut self) {
        if !self.buffer.is_empty() {
            self.pending.push(std::mem::take(&mut self.buffer));
            self.buffer = Vec::with_capacity(CHUNK_BYTES);
        }
    }

    fn drain_pending(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }

        let chunks = std::mem::take(&mut self.pending);
        let level = self.level;
        let compressed: Vec<io::Result<Vec<u8>>> = chunks
            .into_par_iter()
            .map_init(
                || Compressor::new(level),
                |compressor, chunk| {
                    let bound = compressor.gzip_compress_bound(chunk.len());
                    let mut output = vec![0; bound];
                    let written = compressor
                        .gzip_compress(&chunk, &mut output)
                        .map_err(io::Error::other)?;
                    output.truncate(written);
                    Ok(output)
                },
            )
            .collect();

        // Indexed parallel collection preserves source order, which is required
        // because concatenated members reconstruct one byte stream.
        for chunk in compressed {
            self.inner.write_all(&chunk?)?;
            self.wrote_member = true;
        }
        Ok(())
    }
}

impl<W: Write> Write for ParallelGzipWriter<W> {
    fn write(&mut self, mut input: &[u8]) -> io::Result<usize> {
        let input_len = input.len();
        while !input.is_empty() {
            let available = CHUNK_BYTES - self.buffer.len();
            let take = available.min(input.len());
            self.buffer.extend_from_slice(&input[..take]);
            input = &input[take..];

            if self.buffer.len() == CHUNK_BYTES {
                self.queue_buffer();
                if self.pending.len() == MAX_PENDING_CHUNKS {
                    self.drain_pending()?;
                }
            }
        }
        Ok(input_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.queue_buffer();
        self.drain_pending()?;
        self.inner.flush()
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read};

    use flate2::read::MultiGzDecoder;

    use super::*;

    fn decode(bytes: &[u8]) -> Vec<u8> {
        let mut decoded = Vec::new();
        MultiGzDecoder::new(Cursor::new(bytes))
            .read_to_end(&mut decoded)
            .unwrap();
        decoded
    }

    #[test]
    fn arbitrary_write_boundaries_round_trip() {
        let mut writer = ParallelGzipWriter::new(Vec::new(), 4).unwrap();
        writer.write_all(b"first").unwrap();
        writer.write_all(b"-second").unwrap();
        let compressed = writer.finish().unwrap();
        assert_eq!(decode(&compressed), b"first-second");
    }

    #[test]
    fn writes_after_flush_form_one_decoded_stream() {
        let mut writer = ParallelGzipWriter::new(Vec::new(), 4).unwrap();
        writer.write_all(b"first").unwrap();
        writer.flush().unwrap();
        writer.write_all(b"-second").unwrap();
        let compressed = writer.finish().unwrap();
        assert_eq!(decode(&compressed), b"first-second");
    }

    #[test]
    fn empty_output_is_a_valid_gzip_stream() {
        let writer = ParallelGzipWriter::new(Vec::new(), 4).unwrap();
        let compressed = writer.finish().unwrap();
        assert!(!compressed.is_empty());
        assert_eq!(decode(&compressed), b"");
    }

    #[test]
    fn multiple_parallel_batches_preserve_order() {
        let input_len = CHUNK_BYTES * (MAX_PENDING_CHUNKS + 2) + 17;
        let input: Vec<u8> = (0..input_len).map(|index| (index % 251) as u8).collect();
        let mut writer = ParallelGzipWriter::new(Vec::new(), 4).unwrap();
        writer.write_all(&input).unwrap();
        let compressed = writer.finish().unwrap();
        assert_eq!(decode(&compressed), input);
    }

    #[test]
    fn invalid_level_is_rejected() {
        assert_eq!(
            ParallelGzipWriter::new(Vec::new(), 13)
                .err()
                .unwrap()
                .kind(),
            io::ErrorKind::InvalidInput
        );
    }

    #[derive(Debug)]
    struct FailingWriter;

    impl Write for FailingWriter {
        fn write(&mut self, _input: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "injected output failure",
            ))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn downstream_write_failure_is_reported() {
        let mut writer = ParallelGzipWriter::new(FailingWriter, 4).unwrap();
        writer.write_all(b"not silently discarded").unwrap();
        assert_eq!(
            writer.finish().unwrap_err().kind(),
            io::ErrorKind::BrokenPipe
        );
    }
}
