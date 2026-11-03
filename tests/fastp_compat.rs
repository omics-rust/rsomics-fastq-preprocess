use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn installed_version() -> Option<String> {
    let output = Command::new("fastp").arg("--version").output().ok()?;
    output.status.success().then(|| {
        String::from_utf8_lossy(if output.stdout.is_empty() {
            &output.stderr
        } else {
            &output.stdout
        })
        .trim()
        .to_owned()
    })
}

fn require_or_skip() -> bool {
    if let Some(version) = installed_version() {
        assert_eq!(version, "fastp 1.3.6", "live oracle version drifted");
        return true;
    }
    assert_ne!(
        std::env::var("RSOMICS_REQUIRE_FASTP").as_deref(),
        Ok("1"),
        "fastp is required but not on PATH"
    );
    eprintln!("skipping live fastp differential; frozen goldens still run");
    false
}

#[test]
fn quality_windows_preserve_actual_se_and_pe_reads() {
    if !require_or_skip() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let bytes = b"@good\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@front\nACGTACGTACGTACGT\n+\n!!!!IIIIIIIIIIII\n@middle\nACGTACGTACGTACGT\n+\nIIII!!!!IIIIIIII\n@tail\nACGTACGTACGTACGT\n+\nIIIIIIII!!!!!!!!\n@short\nACGT\n+\nIIII\n@poly\nACGTGGGGGGGGGGGG\n+\nIIIIIIII!!!!!!!!\n@poly_empty\nGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIII\n";
    let modes: &[&[&str]] = &[
        &["--cut_front"],
        &["--cut_tail"],
        &["--cut_right"],
        &["--cut_front", "--cut_right"],
        &["--cut_front", "--cut_tail"],
        &[
            "--cut_front",
            "--cut_tail",
            "--trim_front1",
            "2",
            "--trim_tail1",
            "1",
        ],
        &["--cut_right", "--trim_poly_g", "--poly_g_min_len", "4"],
        &["--cut_right", "--trim_poly_x", "--poly_x_min_len", "4"],
        &[
            "--cut_front",
            "--cut_tail",
            "--cut_right",
            "--cut_window_size",
            "3",
            "--cut_mean_quality",
            "15",
            "--cut_front_window_size",
            "2",
            "--cut_right_mean_quality",
            "20",
        ],
    ];
    for encoding in [33, 64] {
        let mut input_bytes = bytes.to_vec();
        if encoding == 64 {
            let mut line = 0;
            for byte in &mut input_bytes {
                if *byte == b'\n' {
                    line += 1;
                } else if line % 4 == 3 {
                    *byte += 31;
                }
            }
        }
        let input = directory.path().join(format!("input-{encoding}.fastq"));
        std::fs::write(&input, &input_bytes).unwrap();
        for (index, mode) in modes.iter().enumerate() {
            let theirs = directory
                .path()
                .join(format!("source-{encoding}-{index}.fastq"));
            let mut flags = vec!["-A", "-G", "-Q", "-L"];
            if mode.contains(&"--trim_poly_g") {
                flags.retain(|flag| *flag != "-G");
            }
            flags.extend_from_slice(mode);
            if encoding == 64 {
                flags.push("-6");
            }
            run_fastp(&input, &theirs, &flags);
            for workers in ["1", "4"] {
                let ours = directory
                    .path()
                    .join(format!("ours-{encoding}-{index}-{workers}.fastq"));
                let mut options = vec![
                    "--threads".to_owned(),
                    workers.to_owned(),
                    if workers == "1" { "run" } else { "trim" }.to_owned(),
                    "--phred-offset".to_owned(),
                    encoding.to_string(),
                ];
                if workers == "1" {
                    options.extend(["-Q".to_owned(), "-L".to_owned()]);
                }
                options.extend(mode.iter().map(|flag| flag.replace('_', "-")));
                run_ours(
                    &input,
                    &ours,
                    &options.iter().map(String::as_str).collect::<Vec<_>>(),
                );
                assert_eq!(
                    std::fs::read(&ours).unwrap(),
                    std::fs::read(&theirs).unwrap(),
                    "encoding={encoding}, mode={mode:?}, workers={workers}"
                );
            }
        }
    }
    let left = directory.path().join("left.fastq");
    let right = directory.path().join("right.fastq");
    let mut first = Vec::new();
    let mut second = Vec::new();
    for record in bytes.split(|byte| *byte == b'@').skip(1) {
        let end = record.iter().position(|byte| *byte == b'\n').unwrap();
        for (output, mate) in [(&mut first, b'1'), (&mut second, b'2')] {
            output.push(b'@');
            output.extend_from_slice(&record[..end]);
            output.extend_from_slice(&[b'/', mate]);
            output.extend_from_slice(&record[end..]);
        }
    }
    let second = String::from_utf8(second)
        .unwrap()
        .replace("IIIIIIIIIIIIIIII", "!!!!!!!!!!!!!!!!");
    std::fs::write(&left, first).unwrap();
    std::fs::write(&right, second).unwrap();
    let source1 = directory.path().join("pair-source1.fastq");
    let source2 = directory.path().join("pair-source2.fastq");
    run_fastp_pair(
        &left,
        &right,
        &source1,
        &source2,
        &["-A", "-G", "-Q", "-L", "--cut_front", "--cut_right"],
    );
    for workers in ["1", "4"] {
        let ours1 = directory.path().join(format!("pair-ours1-{workers}.fastq"));
        let ours2 = directory.path().join(format!("pair-ours2-{workers}.fastq"));
        let status = Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
            .args([
                "--threads",
                workers,
                "run",
                "-Q",
                "-L",
                "--cut-front",
                "--cut-right",
                "-i",
            ])
            .arg(&left)
            .arg("-I")
            .arg(&right)
            .arg("-o")
            .arg(&ours1)
            .arg("-O")
            .arg(&ours2)
            .status()
            .unwrap();
        assert!(status.success());
        assert_eq!(
            std::fs::read(ours1).unwrap(),
            std::fs::read(&source1).unwrap()
        );
        assert_eq!(
            std::fs::read(ours2).unwrap(),
            std::fs::read(&source2).unwrap()
        );
    }
    let poly_input = directory.path().join("poly-only.fastq");
    let poly_output = directory.path().join("poly-empty.fastq");
    std::fs::write(
        &poly_input,
        b"@poly\nGGGGGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
        .args([
            "--json",
            "run",
            "-Q",
            "-L",
            "--cut-right",
            "--trim-poly-g",
            "--poly-g-min-len",
            "4",
            "-i",
        ])
        .arg(poly_input)
        .arg("-o")
        .arg(&poly_output)
        .output()
        .unwrap();
    assert!(result.status.success());
    let report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["result"]["reads_out"], 0);
    assert_eq!(report["result"]["filtering"]["reads_failed_too_short"], 1);
    assert_eq!(report["result"]["trimming"]["poly_g_trimmed_bases"], 16);
    assert_eq!(report["result"]["trimming"]["quality_trimmed_bases"], 0);
    assert_eq!(report["result"]["trimming"]["quality_trimmed_reads"], 0);
    assert!(std::fs::read(poly_output).unwrap().is_empty());
}

#[test]
fn quality_cut_parse_failure_does_not_commit_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("malformed.fastq");
    let output = directory.path().join("output.fastq");
    std::fs::write(
        &input,
        b"@good\nACGTACGT\n+\nIIIIIIII\n@bad\nACGTACGT\n+\nIII\n",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
        .args(["--threads", "4", "trim", "--cut-right", "-i"])
        .arg(input)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
    assert!(!String::from_utf8_lossy(&result.stderr).contains("unexpected argument"));
    assert!(!result.stderr.is_empty());
}

fn run_fastp_pair(
    input_r1: &Path,
    input_r2: &Path,
    output_r1: &Path,
    output_r2: &Path,
    extra: &[&str],
) {
    let status = Command::new("fastp")
        .args(["-i"])
        .arg(input_r1)
        .args(["-I"])
        .arg(input_r2)
        .args(["-o"])
        .arg(output_r1)
        .args(["-O"])
        .arg(output_r2)
        .args(extra)
        .args(["--dont_eval_duplication", "-w", "1", "-j"])
        .arg(output_r1.with_extension("json"))
        .args(["-h"])
        .arg(output_r1.with_extension("html"))
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "fastp paired differential failed");
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn phred64_boundaries_and_uppercase_n_match_live_fastp() {
    if !require_or_skip() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    for (index, (bytes, fastp_args, our_args)) in [
        (
            b"@boundary\nACGT\n+\n@A_~\n".as_slice(),
            &["-6", "-A", "-G", "-Q", "-L"][..],
            &["filter", "--phred-offset", "64", "-Q", "-L"][..],
        ),
        (
            b"@lower\nnnnnAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIII\n@upper\nNNNNAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIII\n"
                .as_slice(),
            &["-A", "-G", "-L", "-n", "0"][..],
            &["filter", "-L", "--n-base-limit", "0"][..],
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let input = directory.path().join(format!("{index}.input.fastq"));
        let theirs = directory.path().join(format!("{index}.fastp.fastq"));
        let ours = directory.path().join(format!("{index}.rsomics.fastq"));
        std::fs::write(&input, bytes).unwrap();
        run_fastp(&input, &theirs, fastp_args);
        run_ours(&input, &ours, our_args);
        assert_eq!(
            std::fs::read(ours).unwrap(),
            std::fs::read(theirs).unwrap(),
            "live fastp boundary differential {index} failed"
        );
    }
}

#[test]
fn paired_fixed_trim_inheritance_and_explicit_zero_match_live_fastp() {
    if !require_or_skip() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let input_r1 = directory.path().join("input.r1.fastq");
    let input_r2 = directory.path().join("input.r2.fastq");
    std::fs::write(&input_r1, b"@pair/1\nAACCGGTT\n+\nIIIIIIII\n").unwrap();
    std::fs::write(&input_r2, b"@pair/2\nTTGGCCAA\n+\nIIIIIIII\n").unwrap();

    for (index, (fastp_trim, our_trim)) in [
        (
            &["-f", "2", "-t", "1"][..],
            &["--trim-front1", "2", "--trim-tail1", "1"][..],
        ),
        (
            &["-f", "2", "-t", "1", "-F", "0", "-T", "0"][..],
            &[
                "--trim-front1",
                "2",
                "--trim-tail1",
                "1",
                "--trim-front2",
                "0",
                "--trim-tail2",
                "0",
            ][..],
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let fastp_r1 = directory.path().join(format!("{index}.fastp.r1.fastq"));
        let fastp_r2 = directory.path().join(format!("{index}.fastp.r2.fastq"));
        let ours_r1 = directory.path().join(format!("{index}.ours.r1.fastq"));
        let ours_r2 = directory.path().join(format!("{index}.ours.r2.fastq"));
        let mut fastp_args = vec!["-A", "-G", "-Q", "-L"];
        fastp_args.extend_from_slice(fastp_trim);
        run_fastp_pair(&input_r1, &input_r2, &fastp_r1, &fastp_r2, &fastp_args);
        let status = Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
            .arg("trim")
            .args(our_trim)
            .args(["-i"])
            .arg(&input_r1)
            .args(["-I"])
            .arg(&input_r2)
            .args(["-o"])
            .arg(&ours_r1)
            .args(["-O"])
            .arg(&ours_r2)
            .status()
            .unwrap();
        assert!(status.success());
        assert_eq!(
            std::fs::read(ours_r1).unwrap(),
            std::fs::read(fastp_r1).unwrap()
        );
        assert_eq!(
            std::fs::read(ours_r2).unwrap(),
            std::fs::read(fastp_r2).unwrap()
        );
    }
}

#[test]
fn maximum_length_and_paired_mixed_failure_match_live_fastp() {
    if !require_or_skip() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();

    let length_input = directory.path().join("length.input.fastq");
    let fastp_length = directory.path().join("length.fastp.fastq");
    let ours_length = directory.path().join("length.ours.fastq");
    std::fs::write(
        &length_input,
        b"@four\nACGT\n+\nIIII\n@five\nACGTA\n+\nIIIII\n",
    )
    .unwrap();
    run_fastp(
        &length_input,
        &fastp_length,
        &["-A", "-G", "-Q", "-l", "0", "--length_limit", "4"],
    );
    run_ours(
        &length_input,
        &ours_length,
        &[
            "filter",
            "-Q",
            "--length-required",
            "0",
            "--length-limit",
            "4",
        ],
    );
    assert_eq!(
        std::fs::read(ours_length).unwrap(),
        std::fs::read(fastp_length).unwrap()
    );

    let input_r1 = directory.path().join("mixed.r1.fastq");
    let input_r2 = directory.path().join("mixed.r2.fastq");
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
    let fastp_r1 = directory.path().join("mixed.fastp.r1.fastq");
    let fastp_r2 = directory.path().join("mixed.fastp.r2.fastq");
    let ours_r1 = directory.path().join("mixed.ours.r1.fastq");
    let ours_r2 = directory.path().join("mixed.ours.r2.fastq");
    run_fastp_pair(&input_r1, &input_r2, &fastp_r1, &fastp_r2, &["-A", "-G"]);
    let status = Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
        .args(["filter", "-i"])
        .arg(&input_r1)
        .args(["-I"])
        .arg(&input_r2)
        .args(["-o"])
        .arg(&ours_r1)
        .args(["-O"])
        .arg(&ours_r2)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(
        std::fs::read(ours_r1).unwrap(),
        std::fs::read(fastp_r1).unwrap()
    );
    assert_eq!(
        std::fs::read(ours_r2).unwrap(),
        std::fs::read(fastp_r2).unwrap()
    );
}

fn run_fastp(input: &Path, output: &Path, extra: &[&str]) {
    let status = Command::new("fastp")
        .args(["-i"])
        .arg(input)
        .args(["-o"])
        .arg(output)
        .args(extra)
        .args(["--dont_eval_duplication", "-w", "1", "-j"])
        .arg(output.with_extension("json"))
        .args(["-h"])
        .arg(output.with_extension("html"))
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "fastp failed for {}", input.display());
}

fn run_ours(input: &Path, output: &Path, args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-preprocess"))
        .args(args)
        .args(["-i"])
        .arg(input)
        .args(["-o"])
        .arg(output)
        .status()
        .unwrap();
    assert!(status.success(), "rsomics failed for {}", input.display());
}

#[test]
fn selected_trim_and_filter_semantics_match_live_fastp() {
    if !require_or_skip() {
        return;
    }
    let directory = tempfile::tempdir().unwrap();
    let cases: &[(&str, &[&str], &[&str])] = &[
        ("filter.fastq", &["-A", "-G"], &["filter"]),
        (
            "complexity.fastq",
            &["-A", "-G", "-Q", "-L", "-y", "-Y", "30"],
            &["filter", "-Q", "-L", "-y", "-Y", "30"],
        ),
        (
            "trim.fastq",
            &["-A", "-g", "-x", "-f", "2", "-t", "2", "-Q", "-L"],
            &["trim", "-f", "2", "--trim-tail1", "2", "-g", "-x"],
        ),
    ];
    for (index, (fixture_name, fastp_args, our_args)) in cases.iter().enumerate() {
        let input = fixture(fixture_name);
        let theirs = directory.path().join(format!("{index}.fastp.fastq"));
        let ours = directory.path().join(format!("{index}.rsomics.fastq"));
        run_fastp(&input, &theirs, fastp_args);
        run_ours(&input, &ours, our_args);
        assert_eq!(
            std::fs::read(&ours).unwrap(),
            std::fs::read(&theirs).unwrap(),
            "live fastp differential failed for {fixture_name}"
        );
    }
}
