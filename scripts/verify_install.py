import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess


def identity(path):
    with path.open("rb") as stream:
        return {"bytes": path.stat().st_size, "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}


def verify(binary, oracle, delivery):
    source = json.loads((delivery / "SOURCE.json").read_bytes())
    work = delivery / "lifecycle"
    work.mkdir()
    calls = []

    def save(name, value):
        (work / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")

    def run(program, *args, stdin=None, ok=True):
        argv = [str(program), *args]
        result = subprocess.run(argv, cwd=work, input=stdin, capture_output=True, timeout=120,
                                env=dict(os.environ, NO_COLOR="1"))
        calls.append({"argv": argv, "returncode": result.returncode})
        (work / f"{len(calls)}.stdout").write_bytes(result.stdout)
        (work / f"{len(calls)}.stderr").write_bytes(result.stderr)
        save("commands.json", calls)
        if not (result.returncode == 0 if ok else result.returncode > 0):
            raise ValueError(f"unexpected exit {result.returncode}: {argv}; {result.stderr.decode(errors='replace')}")
        return result

    assert source["product"] == "rsomics-fastq-preprocess"
    version = f"rsomics-fastq-preprocess {source['version']}".encode()
    assert run(binary, "--version").stdout.strip() == version
    for command in ((), ("run",), ("trim",), ("filter",)):
        help_text = run(binary, *command, "--help").stdout
        if command in (("run",), ("trim",)):
            assert all(flag in help_text for flag in (b"--cut-front", b"--cut-tail", b"--cut-right"))
    upstream = run(oracle, "--version")
    assert (upstream.stdout or upstream.stderr).strip() == b"fastp 1.3.6"

    records = [
        ("good", "IIIIIIIIIIIIIIII"),
        ("front", "!!!!IIIIIIIIIIII"),
        ("middle", "IIII!!!!IIIIIIII"),
        ("tail", "IIIIIIII!!!!!!!!"),
    ]

    def fastq(mate="", offset=0):
        return "".join(f"@{name}{mate}\nACGTACGTACGTACGT\n+\n"
                       + "".join(chr(ord(q) + offset) for q in quality) + "\n"
                       for name, quality in records).encode()

    raw = fastq()
    (work / "input.fastq").write_bytes(raw)
    oracle_flags = ("-A", "-G", "-Q", "-L", "--dont_eval_duplication", "--thread", "1",
                    "--cut_front", "--cut_right")
    quality_flags = ("-Q", "-L", "--cut-front", "--cut-right")
    run(oracle, "-i", "input.fastq", "-o", "expected.fastq", *oracle_flags,
        "--json", "source.json", "--html", "source.html")
    expected = (work / "expected.fastq").read_bytes()
    assert b"@good\n" in expected and b"@middle\n" in expected
    for threads in ("1", "4"):
        output = f"quality-{threads}.fastq"
        run(binary, "--threads", threads, "run", *quality_flags, "-i", "input.fastq", "-o", output)
        assert (work / output).read_bytes() == expected
    assert run(binary, "trim", "--cut-front", "--cut-right", "--phred-offset", "64",
               stdin=fastq(offset=31)).stdout == expected

    for mate in (1, 2):
        (work / f"input{mate}.fastq").write_bytes(fastq(mate=f"/{mate}"))
    run(oracle, "-i", "input1.fastq", "-I", "input2.fastq", "-o", "expected1.fastq",
        "-O", "expected2.fastq", *oracle_flags, "--json", "pair.json", "--html", "pair.html")
    run(binary, "--threads", "4", "run", *quality_flags, "-i", "input1.fastq", "-I", "input2.fastq",
        "-o", "output1.fastq.gz", "-O", "output2.fastq.gz")
    for mate in (1, 2):
        with gzip.open(work / f"output{mate}.fastq.gz", "rb") as stream:
            assert stream.read() == (work / f"expected{mate}.fastq").read_bytes()

    trimmed = run(binary, "trim", "--cut-front", "--cut-tail", stdin=raw).stdout
    filtered = run(binary, "filter", "-Q", "-L", stdin=trimmed).stdout
    assert filtered == run(binary, "run", "-Q", "-L", "--cut-front", "--cut-tail", stdin=raw).stdout
    (work / "poly.fastq").write_bytes(b"@poly\nGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIII\n")
    result = run(binary, "--json", "run", "-Q", "-L", "--cut-right", "--trim-poly-g",
                 "--poly-g-min-len", "4", "-i", "poly.fastq", "-o", "poly-output.fastq")
    report = json.loads(result.stdout)
    assert report["status"] == "ok" and report["result"]["reads_out"] == 0
    assert report["result"]["filtering"]["reads_failed_too_short"] == 1
    assert report["result"]["trimming"]["poly_g_trimmed_bases"] == 16
    assert report["result"]["trimming"]["quality_trimmed_bases"] == 0
    assert (work / "poly-output.fastq").read_bytes() == b""

    (work / "malformed.fastq").write_bytes(raw + b"@bad\nACGT\n+\nIII\n")
    failed = run(binary, "--threads", "4", "trim", "--cut-right", "-i", "malformed.fastq",
                 "-o", "uncommitted.fastq", ok=False)
    assert b"truncated FASTQ quality" in failed.stderr
    assert not (work / "uncommitted.fastq").exists()
    save("VERIFIED.json", {
        "binary": identity(binary), "oracle": identity(oracle), "oracle_version": "fastp 1.3.6",
        "source_head": source["head_sha"], "version": source["version"],
        "native_calls": sum(row["argv"][0] == str(binary) for row in calls), "expected_failures": 1,
        "whole_output_equivalence": True, "quality_report_and_transaction_verified": True,
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for option in ("binary", "oracle", "delivery"):
        parser.add_argument("--" + option, type=Path, required=True)
    options = parser.parse_args()
    verify(*(p.resolve(strict=True) for p in (options.binary, options.oracle, options.delivery)))
