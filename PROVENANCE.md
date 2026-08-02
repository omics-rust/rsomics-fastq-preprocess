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

- `rsomics-common` 0.11.0,
  `5bac25e251cc74c6a43e8302a3a6cc150886a340`;
- `rsomics-help` 0.4.0,
  `61dd6f2ce0cef6d9b4e349af5f96f96a7c95a013`;
- `rsomics-seqio` 0.4.0,
  `0c6ce988d8c90c5bfdaea00c1bcf53ae4aa443dd`.

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
- `rsomics-common` provides the JSON output contract and fallible envelope path;
  serialization, write, newline, and flush failures propagate non-zero.
- FASTQ stdout errors are propagated by this product's `rsomics-seqio` writer.
- Two PE final paths are coordinated with prepare-before-commit, no-clobber,
  and best-effort rollback. They are not an atomic multi-file transaction.
- Per-mate PE failure counters are the rsomics report contract and are not
  presented as fastp JSON compatibility.
- `rsomics-help` renders the product's real nested Clap tree without a second
  help model.

The release API review also made arbitrary caller input safe: public trim and
filter functions reject missing or mismatched qualities, invalid quality
bytes, and invalid configuration rather than panicking or truncating through
`zip`. `PipelineConfig` constructors bind the reported operation to the active
stages. The already parsed internal record path avoids repeating those checks
inside the production hot loop. An integration test proves that `trim` piped
into `filter` produces the same FASTQ stream as `run` for the same stages.

No performance result is inherited as a pass.

## Measured compressed-output gate

The parallel gzip implementation was re-established rather than accepted from
the historical README. Revision
`de07879d1d5ddaab9c5534e50d161ca660ba44e9` uses `libdeflater` 1.25.2
and the existing product Rayon pool to compress ordered 256 KiB gzip members.
It preserves `rsomics-seqio::Writer` validation and serialization, empty gzip
validity, flush-then-write ordering, downstream error propagation, and the
existing no-clobber transaction.

The final production-code gate used revision
`fd04e662426d98f414c51d16a84a2e0eb643e010`, Rust 1.91.0, fastp source
commit `23d6211d4f05d61f561899f1b7702435a4b5d408`, and SRR341550 inputs on
Ubuntu 22.04, Linux 6.8.0, and a two-socket Intel Xeon Gold 6238R host:

- R1 SHA-256:
  `d7a15c1762d64a5434ced0cc665d7f5d167ca81a71e239f8237b9cd490dd7683`;
- R2 SHA-256:
  `18a8e61af21d276dfaf12035307d673e3f52c9f3ac57658ee2f593d1aabeb1a4`;
- filtered R1 byte-stream SHA-256:
  `f13cb655feedf78cf1f3c512675ad73323409f5862b0b3a6e5e3d48e21e6e365`;
- filtered R2 byte-stream SHA-256:
  `452c78a98878e56bf1e5e7728b749e0277e0e14607fa465f7da3e83e551c078c`.
- filtered single-end byte-stream SHA-256:
  `9cc5172922740e7291bdf9fdfadc3d03370665fb0a8d4d4c4c5d4b930c800b58`.

At four threads, five measured paired runs were
`10.914 ± 0.493 s` versus fastp's `14.690 ± 0.715 s`; peak RSS was
31.5 MiB versus 99.2 MiB. The paired slice was 1.35 times faster and used
68% less peak memory on this host.

Single-end output remained byte-identical but was slower: five-run means were
`5.969 ± 0.431 s` versus `5.503 ± 0.862 s` at four threads. Peak RSS was
18.0 MiB versus fastp's 51.1 MiB. This is a memory advantage, not a
single-end throughput claim.

The measured rsomics binary SHA-256 was
`80aae1d1395627ad845f232eeda0652ab7edbcff012cd151ddf7dbf3c772422b`;
the fastp binary SHA-256 was
`8b0521f3d246e13178c49235c0a76230e5ee930fafcaf0db647a4210a4a65966`.
Raw Hyperfine JSON and `/usr/bin/time -v` records are retained under
`benchmarks/linux-x86_64-fastp-1.3.6`. Revision `fd5e1ec` changes only the
broken-pipe integration-test harness, so the production source measured above
is unchanged. Exact-head CI run `30726244422` passed native Linux and macOS on
both `x86_64` and `aarch64` for that production and test tree.

The backend remains private to this product. It is not a new foundation API;
promotion requires a second concrete product consumer with the same contract.
`libdeflater` and its bundled libdeflate C library are MIT licensed.
