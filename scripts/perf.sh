#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <representative.fastq-or-fastq.gz>" >&2
  exit 2
fi

fixture=$1
for tool in hyperfine fastp shasum; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "required tool missing: $tool" >&2
    exit 2
  fi
done

: "${CARGO_TARGET_DIR:?set CARGO_TARGET_DIR to the external rsomics target directory}"
binary="$CARGO_TARGET_DIR/release/rsomics-fastq-preprocess"
if [[ ! -x "$binary" ]]; then
  echo "release binary missing at $binary" >&2
  exit 2
fi

checksum=$(shasum -a 256 "$fixture" | awk '{print $1}')
echo "fixture=$fixture"
echo "sha256=$checksum"
echo "fastp=$(fastp --version 2>&1)"
echo "threads=1"
echo "compression=none (stdout discarded)"
uname -a

hyperfine --warmup 3 --runs 10 \
  "$binary --threads 1 run -i '$fixture' -o - >/dev/null" \
  "fastp -w 1 -i '$fixture' --stdout -A -G --dont_eval_duplication -j /dev/null -h /dev/null >/dev/null 2>/dev/null"

echo "Wall time only; collect peak RSS before any release decision." >&2
