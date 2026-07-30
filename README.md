# rsomics-fastq-preprocess

`rsomics-fastq-preprocess` is one coherent FASTQ preprocessing product. Its
first consumer slice exposes:

- `run`: fixed/poly-tail trimming and read filtering in one traversal;
- `trim`: fixed, poly-G, and poly-X trimming;
- `filter`: quality, average-quality, N, length, and low-complexity filtering.

All three commands share the same strict single-end/paired-end engine. Input
compression is detected from content; output paths ending in `.gz` use gzip.

## Usage

```text
# Single-end pipe, filtering disabled (strict validation and identity output)
cat reads.fastq | rsomics-fastq-preprocess filter -Q -L > checked.fastq

# Fixed and poly-tail trim
rsomics-fastq-preprocess trim \
  -i reads.fastq.gz -o trimmed.fastq.gz \
  --trim-front1 5 --trim-tail1 3 --trim-poly-g --trim-poly-x

# Paired single-pass trim + filter
rsomics-fastq-preprocess run \
  -i R1.fastq.gz -I R2.fastq.gz \
  -o clean.R1.fastq.gz -O clean.R2.fastq.gz \
  --trim-poly-g --qualified-quality-phred 20 --length-required 50

# Common JSON envelope on stdout, FASTQ data in a file
rsomics-fastq-preprocess filter --json \
  -i reads.fastq -o clean.fastq
```

## Stable semantics

The reader is strict FASTQ: malformed records, invalid bytes, sequence/quality
length mismatches, FASTA input, truncated gzip, and paired-record divergence
return non-zero.

Single-end mode accepts `-` for stdin/stdout. Paired mode requires two distinct
file inputs and two distinct file outputs; exact, normalized, hard-link, and
symbolic-link aliases are rejected. Interleaved FASTQ is not part of this
slice. Explicit mate roles must be R1 then R2, whether encoded as `/1` and `/2`
suffixes or CASAVA `1:...` and `2:...` comments. Same-base identifiers without
an explicit role remain valid, and accepted comments are preserved byte for
byte.

File outputs are transactional and no-clobber:

- records are written to temporary files beside their final destinations;
- both PE writers finish before either final name is committed;
- an existing output is never overwritten;
- input/output aliases therefore cannot truncate source data.

The two-file commit is coordinated best-effort, not filesystem-atomic across
both names. If the second commit loses a race, the first is removed; a
concurrent filesystem or rollback failure can still leave one final path and
is reported non-zero.

FASTQ and JSON stdout write and flush failures propagate non-zero, but stdout
cannot be rolled back after a downstream pipe accepts bytes. `--json` requires
FASTQ data to use a file output so reports never corrupt data stdout.

### Trimming

The implemented order is fixed front/tail trimming, poly-G, then poly-X,
matching the retained subset of fastp's processing order. Poly-G and poly-X
defaults use fastp v1.3.6's uppercase byte semantics, mismatch budgets, cut
position, and A/T/C/G tie order. User-supplied mismatch-budget values are an
rsomics extension and are not claimed byte-compatible with fastp. Poly-G is
explicit; instrument-name auto-detection is not performed. In paired mode,
omitted R2 fixed front/tail values inherit R1; explicitly supplied zero remains
zero.

`trim` performs no length filtering by default. Supplying
`--length-required` or `--length-limit` explicitly enables the requested
post-trim length gate.

Adapter auto-detection, static adapter matching, and paired-overlap trimming
are intentionally excluded. The historical rsomics implementations did not
establish exact fastp compatibility and are not advertised by this slice.

### Filtering

Filtering follows fastp v1.3.6 first-failure precedence:

1. unqualified-base percentage;
2. integer mean quality;
3. N-base count;
4. minimum then maximum length;
5. adjacent-base complexity.

Only uppercase `N` contributes to the exact fastp-compatible N limit;
lowercase `n` does not. Percentage boundaries use strict `>` for rejection.
Complexity is
`adjacent_changes / (length - 1)` and rejects length-zero and length-one reads
when enabled. Phred encoding is explicit (`--phred-offset 33` or `64`);
quality bytes below the selected offset are errors. Phred+64 input is always
serialized as Phred+33 by subtracting 31 from every emitted quality byte,
matching `fastp -6`.

For paired input, both mates are evaluated independently, but the pair is
emitted only when both pass. JSON failure counters count individual mates,
while `pairs_in` and `pairs_out` make pair retention explicit. This per-mate
failure accounting is an rsomics JSON contract, not fastp's report schema.

## Compatibility and performance

Frozen FASTQ goldens were generated with fastp 1.3.6. CI checks out exact
upstream commit `23d6211d4f05d61f561899f1b7702435a4b5d408`, builds it from
source, and asserts `fastp --version` before tests. Tests run live
byte-for-byte differential for:

- default quality/N/length filtering;
- low-complexity filtering including zero- and one-base reads;
- fixed, poly-G, and poly-X trimming, including allowed internal and trailing
  poly-G mismatches;
- Phred+64 boundary conversion, uppercase-only N filtering, maximum length,
  PE fixed-trim inheritance/explicit zero, and mixed PE failures.

Revision `de07879d1d5ddaab9c5534e50d161ca660ba44e9` replaces the serial
gzip sink with ordered, independently compressed gzip members. Compression
uses the existing Rayon pool, so it does not add background worker threads
beyond the requested processing pool. Exact-head CI is green on native Linux
and macOS for both `x86_64` and `aarch64`.

A real compressed-output gate used SRR341550 on Ubuntu 22.04 / Linux 6.8,
with two Intel Xeon Gold 6238R CPUs. Input SHA-256 values were
`d7a15c1762d64a5434ced0cc665d7f5d167ca81a71e239f8237b9cd490dd7683`
for R1 and
`18a8e61af21d276dfaf12035307d673e3f52c9f3ac57658ee2f593d1aabeb1a4`
for R2. Times are means and sample standard deviations; RSS is from a
separate `/usr/bin/time -v` run.

| Mode | Threads | rsomics | fastp 1.3.6 | Peak RSS, rsomics / fastp |
|---|---:|---:|---:|---:|
| paired | 1 | 22.308 ± 0.610 s | 39.091 ± 0.894 s | 31.5 / 88.7 MiB |
| paired | 4 | 10.863 ± 0.298 s | 13.891 ± 0.447 s | 31.5 / 101.9 MiB |
| single | 1 | 9.910 ± 0.186 s | 6.849 ± 0.090 s | not recorded |
| single | 4 | 5.360 ± 0.075 s | 4.937 ± 0.721 s | 19.6 / 52.9 MiB |

The paired hot path is faster and uses substantially less memory. The
single-end hot path is not a throughput win on this host; its demonstrated
advantage is lower peak memory. Decompressed single-end and paired outputs are
byte-identical to the aligned fastp slice. The paired gzip files are about
0.07% larger than fastp's and about 1.3% larger than the previous serial
zlib-rs output. `gzip -t`, SeqKit 2.13.0, and fastp 1.3.6 all read the
concatenated-member output.

The Criterion benchmark and `scripts/perf.sh` remain smoke scaffolds rather
than substitutes for this representative external measurement.

## Current exclusions

- adapter sequence/auto-detection and PE overlap;
- sliding-window/leading/trailing quality trimming;
- automatic sequencing-instrument detection;
- interleaved PE stdin/stdout;
- UMI, correction, merge, deduplication, and BBDuk-style filtering;
- fastp's complete JSON/HTML report schema.

The JSON contract here is the smaller rsomics operation report inside the
`rsomics-common` envelope.

## Origin and license

Historical rsomics sources are team-owned. Exact revisions, archive checksum,
dispositions, upstream source, and foundation feedback are in
[`PROVENANCE.md`](PROVENANCE.md).

fastp is MIT licensed. This crate is MIT OR Apache-2.0.
