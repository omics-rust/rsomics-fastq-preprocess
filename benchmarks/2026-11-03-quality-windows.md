# Quality-window workflow measurements

Apple M2,8 cores,8GiB RAM, macOS arm64, Rust1.91 release; actual fastp1.3.6
built from `23d6211d4f05d61f561899f1b7702435a4b5d408` using existing highway,
isa-l and libdeflate. Neither side uses a substituted oracle version. All
inputs, outputs, scratch and caches are external. No RUSTFLAGS override.

Deterministic fixture:1M read pairs,150 bases per mate,317,888,890 bytes per
input. A fixed-seed xorshift emits independent A/C/G/T bases. Five quality
patterns cover good, low leading/trailing, internal low-quality windows and
combined leading/trailing regions; every25th pair has all-low quality. Both
implementations emit960,000 reads per tested mate. This is a synthetic mixed-
quality workload, not a claim about a biological cohort.

Both commands enable front+right windows4/mean20, disable adapter/poly/ordinary
quality/length filters, and use the same explicit threads. SE writes plain
stdout with1 worker. PE writes two gzip streams at level4 with4 workers. Source
HTML/JSON reports remain outside FASTQ output. Full argv and binary/input/source
hashes are in `2026-11-03-quality-samples.json`.

One paired warmup and five alternating-order pairs per workflow, no concurrent
build or benchmark. SE warmup retains full outputs for byte-equivalent hashing;
timed SE stdout is discarded identically for both tools. Every PE pair is
compared by complete decompressed hashes. All checked sequence, quality and
record-order bytes agree. Original outputs and time/stderr logs are preserved.

| Workflow | Native wall median s | fastp wall median s | Native RSS bytes | fastp RSS bytes | Native/source median user CPU s |
|---|---:|---:|---:|---:|---:|
| SE stdout,1 worker |0.455374|0.928148|6,684,672|40,861,696|0.42/1.11|
| PE gzip,4 workers |2.574512|2.022413|32,210,944|100,564,992|3.76/6.56|

Resource gate passes: every measured native RSS is below every measured source
RSS for its workflow. Median reductions are approximately84% and68%. SE wall
throughput improves approximately2.04x. PE gzip wall is approximately27% slower,
despite lower CPU consumption and memory; do not advertise a general speedup.
No fixture/entropy variation was used to hide that negative result.

Native executable SHA-256:
`40f0b95859bc12d41e598531b76c1f26028fb0273a23e1df441dbb5ac90906a4`.
Oracle SHA-256:
`7fb979ed4d571be495a0f5cc210122cc98a4e8f1bb1224d2aa5120947f34401d`.
The measured candidate was built on89d4f53 with working quality changes, before
the0.2 root-version update. No future release binary is independently measured
by this record. All production bodies at measurement are individually hashed.

Original packet:
`/Volumes/KIOXIA/Developments/tmp/fastq-quality-perf-20261103.q8zY938x`.
RESULT SHA-256:
`433749a0327fa291cae0788cb92d67f82b7cf8e7120cdba67675415e56192002`.
