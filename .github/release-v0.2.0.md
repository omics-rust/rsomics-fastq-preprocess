`rsomics-fastq-preprocess` trims and filters single-end and paired-end FASTQ
with `run`, `trim`, and `filter`.

- `run` and `trim` support opt-in `--cut-front`, `--cut-tail`, and `--cut-right`
  quality windows compatible with the tested fastp 1.3.6 workflows.
- Phred+64 input is validated before trimming. Invalid quality bytes now
  fail even if the selected trim would remove them.
- Gzip output, mate validation, and no-clobber file output remain supported.
  Adapter detection and paired-end overlap correction are out of scope.

Download the archive for Linux or macOS on x86_64 or arm64. Check it against
`SHA256SUMS`, extract it, and run `rsomics-fastq-preprocess --version`.
The source `.crate` is also attached. This is a GitHub release, separate
from crates.io publication.

On a 1M-pair synthetic workload on Apple M2, the measured 0.2.0 production
source matched fastp 1.3.6 output. Single-end stdout throughput was about
2×; paired gzip took about 23% longer with 68% lower peak memory. These
figures describe the [measured workload and build](https://github.com/omics-rust/rsomics-world/blob/3ab7a8b61c6aa34cc137e784ad7dabc21d4140b2/.autopilot/state/fastq-phred64-validation-2026-11-05.md),
not a speedup claim for every workload or shipped binary.
