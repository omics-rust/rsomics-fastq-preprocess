use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
}

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn fixture(name: &str) -> PathBuf {
    root().join("tests/golden").join(name)
}

fn run(args: &[&str]) -> std::process::Output {
    Command::new(binary())
        .args(args)
        .current_dir(root())
        .output()
        .expect("run rsomics-fastq-preprocess")
}

#[test]
fn help_is_available_for_the_product_and_each_operation() {
    for args in [
        &["--help"][..],
        &["run", "--help"],
        &["trim", "--help"],
        &["filter", "--help"],
    ] {
        let output = run(args);
        assert!(
            output.status.success(),
            "args={args:?}\nstderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let top = String::from_utf8(run(&["--help"]).stdout).unwrap();
    assert!(top.contains("Global options:"));
    assert!(top.contains("--threads"));
    assert!(top.contains("--json"));
    for absent in ["--seed", "--quiet", "--verbose"] {
        assert!(!top.contains(absent));
    }

    let run_help = String::from_utf8(run(&["run", "--help"]).stdout).unwrap();
    for heading in [
        "Input/output:",
        "Trimming:",
        "Filtering:",
        "Length filtering:",
        "Global options:",
    ] {
        assert!(run_help.contains(heading), "missing {heading:?}");
    }

    let trim_help = String::from_utf8(run(&["trim", "--help"]).stdout).unwrap();
    assert!(trim_help.contains("Output filtering:"));
}

#[test]
fn filter_matches_frozen_fastp_default_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("filtered.fastq");
    let output = run(&[
        "filter",
        "-i",
        fixture("filter.fastq").to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read(output_path).unwrap(),
        std::fs::read(fixture("filter.default.fastq")).unwrap()
    );
}

#[test]
fn trim_and_single_pass_run_match_frozen_fastp_output() {
    for operation in ["trim", "run"] {
        let directory = tempfile::tempdir().unwrap();
        let output_path = directory.path().join(format!("{operation}.fastq"));
        let input_path = fixture("trim.fastq");
        let mut args = vec![
            operation,
            "-i",
            input_path.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
            "-f",
            "2",
            "--trim-tail1",
            "2",
            "-g",
            "-x",
        ];
        if operation == "run" {
            args.push("-Q");
            args.push("-L");
        }
        let output = run(&args);
        assert!(
            output.status.success(),
            "operation={operation}\nstderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            std::fs::read(output_path).unwrap(),
            std::fs::read(fixture("trim.fastp-1.3.6.fastq")).unwrap()
        );
    }
}

#[test]
fn complexity_short_read_boundary_matches_frozen_fastp_output() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("complexity.fastq");
    let output = run(&[
        "filter",
        "-i",
        fixture("complexity.fastq").to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
        "-Q",
        "-L",
        "-y",
        "-Y",
        "30",
    ]);
    assert!(output.status.success());
    assert_eq!(
        std::fs::read(output_path).unwrap(),
        std::fs::read(fixture("complexity.fastp-1.3.6.fastq")).unwrap()
    );
}

#[test]
fn paired_filter_is_synchronized_and_json_is_populated() {
    let directory = tempfile::tempdir().unwrap();
    let output_r1 = directory.path().join("filtered.r1.fastq");
    let output_r2 = directory.path().join("filtered.r2.fastq");
    let output = run(&[
        "--json",
        "filter",
        "-i",
        fixture("paired.r1.fastq").to_str().unwrap(),
        "-I",
        fixture("paired.r2.fastq").to_str().unwrap(),
        "-o",
        output_r1.to_str().unwrap(),
        "-O",
        output_r2.to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        std::fs::read(output_r1).unwrap(),
        std::fs::read(fixture("paired.default.r1.fastq")).unwrap()
    );
    assert_eq!(
        std::fs::read(output_r2).unwrap(),
        std::fs::read(fixture("paired.default.r2.fastq")).unwrap()
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema_version"], "1.0");
    assert_eq!(value["result"]["mode"], "PE");
    assert_eq!(value["result"]["reads_in"], 8);
    assert_eq!(value["result"]["reads_out"], 2);
    assert_eq!(value["result"]["pairs_in"], 4);
    assert_eq!(value["result"]["pairs_out"], 1);
    assert_eq!(value["result"]["filtering"]["reads_failed_quality"], 4);
}

#[test]
fn thread_counts_preserve_single_and_paired_results() {
    for paired in [false, true] {
        let directory = tempfile::tempdir().unwrap();
        let input_r1 = if paired {
            fixture("paired.r1.fastq")
        } else {
            fixture("filter.fastq")
        };
        let input_r2 = paired.then(|| fixture("paired.r2.fastq"));
        let mut results = Vec::new();
        for threads in ["1", "4"] {
            let output_r1 = directory.path().join(format!("{threads}.r1.fastq"));
            let output_r2 = directory.path().join(format!("{threads}.r2.fastq"));
            let mut args = vec![
                "--threads",
                threads,
                "--json",
                "filter",
                "-i",
                input_r1.to_str().unwrap(),
                "-o",
                output_r1.to_str().unwrap(),
            ];
            if let Some(input_r2) = input_r2.as_ref() {
                args.extend([
                    "-I",
                    input_r2.to_str().unwrap(),
                    "-O",
                    output_r2.to_str().unwrap(),
                ]);
            }
            let output = run(&args);
            assert!(
                output.status.success(),
                "paired={paired}, threads={threads}, stderr={}",
                String::from_utf8_lossy(&output.stderr)
            );
            let mut report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
            let result = report["result"].as_object_mut().unwrap();
            result.remove("output_r1");
            result.remove("output_r2");
            results.push((
                std::fs::read(output_r1).unwrap(),
                paired.then(|| std::fs::read(output_r2).unwrap()),
                report,
            ));
        }
        assert_eq!(results[0], results[1], "paired={paired}");
    }
}

#[test]
fn single_end_stdin_and_stdout_form_an_identity_pipeline() {
    let bytes = std::fs::read(fixture("filter.fastq")).unwrap();
    let mut child = Command::new(binary())
        .args(["filter", "-Q", "-L"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "stderr={}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, bytes);
}

#[test]
fn trim_and_filter_stream_composition_matches_run() {
    let trim = run(&[
        "trim",
        "-i",
        fixture("trim.fastq").to_str().unwrap(),
        "-f",
        "2",
        "--trim-tail1",
        "2",
        "-g",
        "-x",
    ]);
    assert!(trim.status.success());

    let mut filter = Command::new(binary())
        .args(["filter", "-Q", "-L"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    filter
        .stdin
        .take()
        .unwrap()
        .write_all(&trim.stdout)
        .unwrap();
    let filtered = filter.wait_with_output().unwrap();
    assert!(filtered.status.success());

    let combined = run(&[
        "run",
        "-i",
        fixture("trim.fastq").to_str().unwrap(),
        "-f",
        "2",
        "--trim-tail1",
        "2",
        "-g",
        "-x",
        "-Q",
        "-L",
    ]);
    assert!(combined.status.success());
    assert_eq!(filtered.stdout, combined.stdout);
}

#[test]
fn malformed_input_leaves_no_final_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("malformed.fastq");
    let output_path = directory.path().join("must-not-exist.fastq");
    std::fs::write(&input, b"@one\nACGT\n+\n!!!\n").unwrap();
    let output = run(&[
        "filter",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(!output_path.exists());
    assert!(String::from_utf8_lossy(&output.stderr).contains("truncated FASTQ quality"));
}

#[test]
fn paired_identifier_or_count_mismatch_rolls_back_both_outputs() {
    for second_bytes in [
        b"@other/2\nACGT\n+\nIIII\n".as_slice(),
        b"@pair/2\nACGT\n+\nIIII\n@extra/2\nACGT\n+\nIIII\n",
    ] {
        let directory = tempfile::tempdir().unwrap();
        let first = directory.path().join("r1.fastq");
        let second = directory.path().join("r2.fastq");
        let output_r1 = directory.path().join("out.r1.fastq");
        let output_r2 = directory.path().join("out.r2.fastq");
        std::fs::write(&first, b"@pair/1\nACGT\n+\nIIII\n").unwrap();
        std::fs::write(&second, second_bytes).unwrap();
        let output = run(&[
            "filter",
            "-i",
            first.to_str().unwrap(),
            "-I",
            second.to_str().unwrap(),
            "-o",
            output_r1.to_str().unwrap(),
            "-O",
            output_r2.to_str().unwrap(),
            "-Q",
            "-L",
        ]);
        assert!(!output.status.success());
        assert!(!output_r1.exists());
        assert!(!output_r2.exists());
    }
}

#[test]
fn existing_or_aliased_output_is_never_truncated() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("existing.fastq");
    std::fs::write(&output_path, b"keep this").unwrap();
    let output = run(&[
        "filter",
        "-i",
        fixture("filter.fastq").to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert_eq!(std::fs::read(&output_path).unwrap(), b"keep this");

    let aliased = fixture("filter.fastq");
    let before = std::fs::read(&aliased).unwrap();
    let output = run(&[
        "filter",
        "-i",
        aliased.to_str().unwrap(),
        "-o",
        aliased.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert_eq!(std::fs::read(aliased).unwrap(), before);
}

#[test]
fn gzip_output_round_trips_through_seqio_content_detection() {
    let directory = tempfile::tempdir().unwrap();
    let output_path = directory.path().join("filtered.data.gz");
    let output = run(&[
        "filter",
        "-i",
        fixture("filter.fastq").to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
    ]);
    assert!(output.status.success());
    assert_eq!(&std::fs::read(&output_path).unwrap()[..2], &[0x1f, 0x8b]);
    let mut reader = rsomics_seqio::open_path(&output_path).unwrap();
    let mut count = 0;
    while reader.read_record().unwrap().is_some() {
        count += 1;
    }
    assert_eq!(count, 3);
}

#[test]
fn fully_filtered_gzip_output_is_valid_and_empty() {
    let directory = tempfile::tempdir().unwrap();
    let input_path = directory.path().join("short.fastq");
    let output_path = directory.path().join("filtered.fastq.gz");
    std::fs::write(&input_path, b"@short\nACGT\n+\nIIII\n").unwrap();

    let output = run(&[
        "filter",
        "-i",
        input_path.to_str().unwrap(),
        "-o",
        output_path.to_str().unwrap(),
        "--length-required",
        "5",
    ]);
    assert!(output.status.success());
    let compressed = std::fs::read(&output_path).unwrap();
    assert_eq!(&compressed[..2], &[0x1f, 0x8b]);
    let mut decoded = Vec::new();
    flate2::read::MultiGzDecoder::new(compressed.as_slice())
        .read_to_end(&mut decoded)
        .unwrap();
    assert!(decoded.is_empty());
}

#[test]
fn json_and_fastq_cannot_share_stdout() {
    let output = run(&["--json", "filter", "-i", "tests/golden/filter.fastq"]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("stdout"));
}

#[test]
fn trim_is_pure_by_default_and_length_filtering_is_explicit() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("short.fastq");
    let pure_output = directory.path().join("pure.fastq");
    let filtered_output = directory.path().join("filtered.fastq");
    let bytes = b"@short\nACGT\n+\nIIII\n";
    std::fs::write(&input, bytes).unwrap();

    let pure = run(&[
        "trim",
        "-i",
        input.to_str().unwrap(),
        "-o",
        pure_output.to_str().unwrap(),
    ]);
    assert!(pure.status.success());
    assert_eq!(std::fs::read(pure_output).unwrap(), bytes);

    let filtered = run(&[
        "trim",
        "-i",
        input.to_str().unwrap(),
        "-o",
        filtered_output.to_str().unwrap(),
        "--length-required",
        "5",
    ]);
    assert!(filtered.status.success());
    assert!(std::fs::read(filtered_output).unwrap().is_empty());
}

#[test]
fn phred64_output_is_normalized_to_phred33() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("phred64.fastq");
    std::fs::write(&input, b"@boundary\nACGT\n+\n@A_~\n").unwrap();

    for (index, args) in [
        &["filter", "--phred-offset", "64", "-Q", "-L"][..],
        &["trim", "--phred-offset", "64"][..],
        &["run", "--phred-offset", "64", "-Q", "-L"][..],
    ]
    .into_iter()
    .enumerate()
    {
        let output_path = directory.path().join(format!("{index}.phred33.fastq"));
        let mut command = args.to_vec();
        command.extend([
            "-i",
            input.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ]);
        let output = run(&command);
        assert!(
            output.status.success(),
            "args={args:?}, stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            std::fs::read(output_path).unwrap(),
            b"@boundary\nACGT\n+\n!\"@_\n"
        );
    }
}

#[test]
fn paired_input_aliases_are_rejected_before_output_creation() {
    for alias in ["exact", "normalized", "hardlink", "symlink"] {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("r1.fastq");
        std::fs::write(&input, b"@pair/1\nACGT\n+\nIIII\n").unwrap();
        let second = match alias {
            "exact" => input.clone(),
            "normalized" => directory.path().join(".").join("r1.fastq"),
            "hardlink" => {
                let path = directory.path().join("r2.fastq");
                std::fs::hard_link(&input, &path).unwrap();
                path
            }
            "symlink" => {
                #[cfg(unix)]
                {
                    let path = directory.path().join("r2.fastq");
                    std::os::unix::fs::symlink(&input, &path).unwrap();
                    path
                }
                #[cfg(not(unix))]
                {
                    continue;
                }
            }
            _ => unreachable!(),
        };
        let output_r1 = directory.path().join("out.r1.fastq");
        let output_r2 = directory.path().join("out.r2.fastq");
        let result = run(&[
            "filter",
            "-i",
            input.to_str().unwrap(),
            "-I",
            second.to_str().unwrap(),
            "-o",
            output_r1.to_str().unwrap(),
            "-O",
            output_r2.to_str().unwrap(),
        ]);
        assert_eq!(
            result.status.code(),
            Some(2),
            "alias={alias}, stderr={}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(!output_r1.exists());
        assert!(!output_r2.exists());
    }
}

#[test]
fn mate_roles_are_enforced_and_casava_comments_are_preserved() {
    for (left_id, right_id) in [
        ("pair/2", "pair/1"),
        ("pair/1", "pair/1"),
        ("pair 2:N:0:1", "pair 1:N:0:1"),
    ] {
        let directory = tempfile::tempdir().unwrap();
        let input_r1 = directory.path().join("r1.fastq");
        let input_r2 = directory.path().join("r2.fastq");
        let output_r1 = directory.path().join("out.r1.fastq");
        let output_r2 = directory.path().join("out.r2.fastq");
        std::fs::write(&input_r1, format!("@{left_id}\nACGT\n+\nIIII\n")).unwrap();
        std::fs::write(&input_r2, format!("@{right_id}\nTGCA\n+\nIIII\n")).unwrap();
        let result = run(&[
            "filter",
            "-i",
            input_r1.to_str().unwrap(),
            "-I",
            input_r2.to_str().unwrap(),
            "-o",
            output_r1.to_str().unwrap(),
            "-O",
            output_r2.to_str().unwrap(),
            "-Q",
            "-L",
        ]);
        assert!(!result.status.success());
        assert!(!output_r1.exists());
        assert!(!output_r2.exists());
    }

    let directory = tempfile::tempdir().unwrap();
    let input_r1 = directory.path().join("casava.r1.fastq");
    let input_r2 = directory.path().join("casava.r2.fastq");
    let output_r1 = directory.path().join("out.r1.fastq");
    let output_r2 = directory.path().join("out.r2.fastq");
    let r1 = b"@instrument:1:flowcell:1 1:N:0:ACGT\nACGT\n+\nIIII\n";
    let r2 = b"@instrument:1:flowcell:1 2:N:0:ACGT\nTGCA\n+\nIIII\n";
    std::fs::write(&input_r1, r1).unwrap();
    std::fs::write(&input_r2, r2).unwrap();
    let result = run(&[
        "filter",
        "-i",
        input_r1.to_str().unwrap(),
        "-I",
        input_r2.to_str().unwrap(),
        "-o",
        output_r1.to_str().unwrap(),
        "-O",
        output_r2.to_str().unwrap(),
        "-Q",
        "-L",
    ]);
    assert!(result.status.success());
    assert_eq!(std::fs::read(output_r1).unwrap(), r1);
    assert_eq!(std::fs::read(output_r2).unwrap(), r2);
}

#[test]
fn paired_json_counts_each_mates_distinct_failure() {
    let directory = tempfile::tempdir().unwrap();
    let input_r1 = directory.path().join("mixed.r1.fastq");
    let input_r2 = directory.path().join("mixed.r2.fastq");
    let output_r1 = directory.path().join("out.r1.fastq");
    let output_r2 = directory.path().join("out.r2.fastq");
    std::fs::write(
        &input_r1,
        b"@pair/1\nACGTACGTACGTACGT\n+\n!!!!!!!!!!!!!!!!\n",
    )
    .unwrap();
    std::fs::write(
        &input_r2,
        b"@pair/2\nNNNNNNAAAAAAAAAA\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let result = run(&[
        "--json",
        "filter",
        "-i",
        input_r1.to_str().unwrap(),
        "-I",
        input_r2.to_str().unwrap(),
        "-o",
        output_r1.to_str().unwrap(),
        "-O",
        output_r2.to_str().unwrap(),
    ]);
    assert!(result.status.success());
    assert!(std::fs::read(output_r1).unwrap().is_empty());
    assert!(std::fs::read(output_r2).unwrap().is_empty());
    let envelope: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(envelope["result"]["filtering"]["reads_failed_quality"], 1);
    assert_eq!(envelope["result"]["filtering"]["reads_failed_n_bases"], 1);
    assert_eq!(envelope["result"]["pairs_in"], 1);
    assert_eq!(envelope["result"]["pairs_out"], 0);
}

#[test]
fn gzip_stdin_is_decoded_and_truncation_leaves_no_output() {
    let raw = b"@one\nACGT\n+\nIIII\n";
    let encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    let mut writer = rsomics_seqio::Writer::new(encoder, rsomics_seqio::Format::Fastq);
    writer
        .write_record(rsomics_seqio::Record {
            id: b"one",
            seq: b"ACGT",
            qual: Some(b"IIII"),
        })
        .unwrap();
    let encoded = writer.finish_into_inner().unwrap().finish().unwrap();

    let directory = tempfile::tempdir().unwrap();
    let decoded = directory.path().join("decoded.fastq");
    let mut child = Command::new(binary())
        .args(["filter", "-Q", "-L", "-o"])
        .arg(&decoded)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(&encoded).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    assert_eq!(std::fs::read(decoded).unwrap(), raw);

    let truncated_output = directory.path().join("truncated.fastq");
    let mut child = Command::new(binary())
        .args(["filter", "-Q", "-L", "-o"])
        .arg(&truncated_output)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&encoded[..encoded.len() - 6])
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(!truncated_output.exists());
}

#[cfg(unix)]
#[test]
fn hardlink_and_symlink_output_aliases_preserve_input() {
    for alias in ["hardlink", "symlink"] {
        let directory = tempfile::tempdir().unwrap();
        let input = directory.path().join("input.fastq");
        let output_path = directory.path().join("output.fastq");
        let bytes = b"@one\nACGT\n+\nIIII\n";
        std::fs::write(&input, bytes).unwrap();
        match alias {
            "hardlink" => std::fs::hard_link(&input, &output_path).unwrap(),
            "symlink" => std::os::unix::fs::symlink(&input, &output_path).unwrap(),
            _ => unreachable!(),
        }
        let result = run(&[
            "filter",
            "-i",
            input.to_str().unwrap(),
            "-o",
            output_path.to_str().unwrap(),
        ]);
        assert_eq!(result.status.code(), Some(2));
        assert_eq!(std::fs::read(input).unwrap(), bytes);
    }
}

#[cfg(unix)]
#[test]
fn fastq_stdout_write_failure_is_nonzero() {
    let mut child = Command::new(binary())
        .args([
            "filter",
            "-i",
            fixture("filter.fastq").to_str().unwrap(),
            "-Q",
            "-L",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Broken pipe")
            || String::from_utf8_lossy(&output.stderr).contains("broken pipe")
    );
}
