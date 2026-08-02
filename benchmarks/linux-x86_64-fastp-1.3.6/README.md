# Linux x86_64 release gate

This directory retains the raw four-thread Hyperfine and `/usr/bin/time -v`
records for the release-candidate production tree at
`fd04e662426d98f414c51d16a84a2e0eb643e010`.

- Host: Ubuntu 22.04, Linux 6.8.0, two Intel Xeon Gold 6238R CPUs
- Rust: 1.91.0
- Oracle: fastp 1.3.6 at `23d6211d4f05d61f561899f1b7702435a4b5d408`
- Warmup: one run per command
- Samples: five runs per command
- Input: SRR341550 compressed paired FASTQ

Input SHA-256:

- R1: `d7a15c1762d64a5434ced0cc665d7f5d167ca81a71e239f8237b9cd490dd7683`
- R2: `18a8e61af21d276dfaf12035307d673e3f52c9f3ac57658ee2f593d1aabeb1a4`

Decompressed outputs from rsomics and fastp matched at each corresponding
position:

- paired R1: `f13cb655feedf78cf1f3c512675ad73323409f5862b0b3a6e5e3d48e21e6e365`
- paired R2: `452c78a98878e56bf1e5e7728b749e0277e0e14607fa465f7da3e83e551c078c`
- single R1: `9cc5172922740e7291bdf9fdfadc3d03370665fb0a8d4d4c4c5d4b930c800b58`

Binary SHA-256:

- rsomics: `80aae1d1395627ad845f232eeda0652ab7edbcff012cd151ddf7dbf3c772422b`
- fastp: `8b0521f3d246e13178c49235c0a76230e5ee930fafcaf0db647a4210a4a65966`

The absolute paths in the raw records identify the isolated external-disk
benchmark workspace. Output FASTQ files are omitted because the source input
is public, the commands are preserved in the raw records, and the decompressed
content checksums above bind the compatibility result.
