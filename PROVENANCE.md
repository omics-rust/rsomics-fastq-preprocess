# Migration and compatibility record

## Source assets

| Target | Historical source | Revision / identity | Disposition |
|---|---|---|---|
| strict chunk engine | `rsomics-fastp` 0.2.0 registry archive | SHA-256 `f47ccbf47c416889772494777592347c93f0ef72af9114f4879b7f69f5bb4ab2` | architecture seed only; replace parser, writer, filter precedence, PE verdict, and reporting |
| fixed/poly trim | `rsomics-fastq-trim` | `30b07d34d02864518efc02f1e7a5f881769fcd3f` | refactor pure transforms; re-establish fastp behavior |
| whole-read filter | `rsomics-fastq-filter` | `bdc2824778ed1b8fc185b173d841f63d8b9e6759` | retain boundary fixtures; correct precedence and transactional I/O |
| quality trimming | `rsomics-fastq-quality` | `711e32cac0f45072d038171f42d6d4a831112cff` | evidence only; excluded until Trimmomatic versus fastp semantics have separate oracles |
| complexity | `rsomics-fastq-complexity` | `22948ed73a0f1b1b26cba7494b25c79fb8067811` | refactor adjacent-change metric; correct length-zero/one behavior |
| parallel gzip output | `rsomics-fqgz` | implementation introduced at `0ebfe46`; retained source head `c5e7de1` | refactor then merge as a private product module; retain strict seqio serialization and transactional output |
| old combined CLI/report | `rsomics-fastp` | `2d1d68d25d10c087aae4548b922b9b1b4e9b80ad` | fixtures and option inventory only |

All listed rsomics code is team-owned.

Foundation revisions exercised by this slice:

- `rsomics-seqio`
  `ce9c5514c23573a64406e1ff9ad02edfa4d02d31`;
- `rsomics-common`
  `1c51f7d0b356683697942d9c6a0f60585e0dc8a9`.

`rsomics-kmer` is deliberately not a dependency. fastp's selected complexity
contract is adjacent-base change fraction, not a k-mer algorithm; adding the
foundation would create an unused or semantically incorrect dependency.

## Upstream oracle

- Binary: fastp 1.3.6 (`fastp --version`).
- Exact source/tag commit:
  `OpenGene/fastp` `v1.3.6`,
  `23d6211d4f05d61f561899f1b7702435a4b5d408`.
- Trim implementation:
  `OpenGene/fastp` tag `v1.3.6`, `src/polyx.cpp`.
- Filter implementation:
  `OpenGene/fastp` tag `v1.3.6`, `src/filter.cpp`.
- Live command shapes are encoded in `tests/fastp_compat.rs`.
- Frozen byte-exact outputs are under `tests/golden/`.

Verified operations are limited to explicit fixed/poly-G/poly-X trimming and
quality/N/length/complexity filtering. No compatibility claim is inherited
from an old README.

## Corrected inherited behavior

- Empty and one-base reads fail enabled complexity filtering.
- Filter precedence is unqualified percentage, integer average quality, N
  count, minimum/maximum length, then complexity.
- PE inputs validate record count and normalized mate identity.
- Both mates are classified independently; R1 does not silently decide a
  pair-wide failure bucket.
- File output is staged and committed only after successful parse/write.
- Existing output is no-clobber, preventing input/output aliases.
- Metrics never saturate or substitute `u64::MAX`; supported targets are
  64-bit and conversions assert that static platform invariant.
- Phred encoding is explicit and bytes below the selected offset fail.
- Phred+64 output is normalized to Phred+33 by subtracting 31 from every
  emitted quality byte.
- Omitted R2 fixed-trim values inherit R1; explicit zero is retained.
- Exact fastp N filtering counts uppercase `N` only.
- `trim` has no implicit length gate; post-trim length filtering is enabled
  only by an explicit length option.
- PE input aliases are rejected across exact, normalized, hard-link, and
  symbolic-link paths. Explicit mate roles must be R1 then R2; valid CASAVA
  comments are preserved.

## Explicit exclusions

Adapter matching and overlap code from the old repositories is not migrated in
this slice. Its defaults and cut choices were not exact fastp behavior.
Sliding-window quality trimming is also excluded because the historical crate
mixed Trimmomatic and fastp contracts.

UMI, correction, merge, deduplication, BBDuk-style filtering, interleaved PE,
instrument auto-detection, and the full fastp report remain outside this
checkpoint.

## Foundation feedback

- `rsomics-seqio` correctly owns strict parsing, content-based gzip/BGZF
  detection, records, and FASTQ serialization.
- Coordinated PE synchronization, output alias policy, and two-output
  transactions remain product policy and are implemented here.
- A future shared quality-encoding type requires the concrete
  `rsomics-fastq-qc` consumer before promotion to `rsomics-seqio`.
- `rsomics-common` provides execution flags and a fallible JSON envelope path;
  serialization, write, newline, and flush failures propagate non-zero.
- FASTQ stdout errors are propagated by this product's `rsomics-seqio` writer.
- Two PE final paths are coordinated with prepare-before-commit, no-clobber,
  and best-effort rollback. They are not an atomic multi-file transaction.
- Per-mate PE failure counters are the rsomics report contract and are not
  presented as fastp JSON compatibility.
- The current `rsomics-help` model duplicates command metadata instead of
  deriving the nested Clap tree; this product retains one Clap source of
  truth.

No performance result is inherited as a pass.

## Measured compressed-output gate

The parallel gzip implementation was re-established rather than accepted from
the historical README. Revision
`de07879d1d5ddaab9c5534e50d161ca660ba44e9` uses `libdeflater` 1.25.2
and the existing product Rayon pool to compress ordered 256 KiB gzip members.
It preserves `rsomics-seqio::Writer` validation and serialization, empty gzip
validity, flush-then-write ordering, downstream error propagation, and the
existing no-clobber transaction.

The Linux `x86_64` gate used Rust 1.91.0, fastp source commit
`23d6211d4f05d61f561899f1b7702435a4b5d408`, and SRR341550 inputs:

- R1 SHA-256:
  `d7a15c1762d64a5434ced0cc665d7f5d167ca81a71e239f8237b9cd490dd7683`;
- R2 SHA-256:
  `18a8e61af21d276dfaf12035307d673e3f52c9f3ac57658ee2f593d1aabeb1a4`;
- filtered R1 byte-stream SHA-256:
  `f13cb655feedf78cf1f3c512675ad73323409f5862b0b3a6e5e3d48e21e6e365`;
- filtered R2 byte-stream SHA-256:
  `452c78a98878e56bf1e5e7728b749e0277e0e14607fa465f7da3e83e551c078c`.

At four threads, ten measured paired runs were
`10.863 ± 0.298 s` versus fastp's `13.891 ± 0.447 s`; peak RSS was
31.5 MiB versus 101.9 MiB. At one thread, five measured paired runs were
`22.308 ± 0.610 s` versus `39.091 ± 0.894 s`; peak RSS was 31.5 MiB
versus 88.7 MiB.

Single-end output remained byte-identical but was slower: five-run means were
`9.910 ± 0.186 s` versus `6.849 ± 0.090 s` at one thread and
`5.360 ± 0.075 s` versus `4.937 ± 0.721 s` at four threads. The measured
four-thread RSS was 19.6 MiB versus fastp's 52.9 MiB. This is a memory
advantage, not a single-end throughput claim.

The final gzip members passed `gzip -t` and were consumed by SeqKit 2.13.0 and
fastp 1.3.6. Their paired file sizes were approximately 0.07% above fastp and
1.3% above the previous serial zlib-rs output. Exact-head CI run
`30551968781` passed native Linux and macOS on both `x86_64` and `aarch64`.

The backend remains private to this product. It is not a new foundation API;
promotion requires a second concrete product consumer with the same contract.
`libdeflater` and its bundled libdeflate C library are MIT licensed.
