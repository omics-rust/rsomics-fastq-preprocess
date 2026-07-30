# Migration and compatibility record

## Source assets

| Target | Historical source | Revision / identity | Disposition |
|---|---|---|---|
| strict chunk engine | `rsomics-fastp` 0.2.0 registry archive | SHA-256 `f47ccbf47c416889772494777592347c93f0ef72af9114f4879b7f69f5bb4ab2` | architecture seed only; replace parser, writer, filter precedence, PE verdict, and reporting |
| fixed/poly trim | `rsomics-fastq-trim` | `30b07d34d02864518efc02f1e7a5f881769fcd3f` | refactor pure transforms; re-establish fastp behavior |
| whole-read filter | `rsomics-fastq-filter` | `bdc2824778ed1b8fc185b173d841f63d8b9e6759` | retain boundary fixtures; correct precedence and transactional I/O |
| quality trimming | `rsomics-fastq-quality` | `711e32cac0f45072d038171f42d6d4a831112cff` | evidence only; excluded until Trimmomatic versus fastp semantics have separate oracles |
| complexity | `rsomics-fastq-complexity` | `22948ed73a0f1b1b26cba7494b25c79fb8067811` | refactor adjacent-change metric; correct length-zero/one behavior |
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
