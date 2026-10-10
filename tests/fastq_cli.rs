use std::{
    fs::{self, File},
    io::Write,
    path::PathBuf,
    process::{Command, Output},
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

use flate2::{Compression, write::GzEncoder};

static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

struct TempFastq {
    path: PathBuf,
}

impl TempFastq {
    fn path_with_suffix(suffix: &str) -> PathBuf {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "seqsource-cli-{}-{timestamp}-{}{suffix}",
            std::process::id(),
            NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn plain(contents: &[u8]) -> Self {
        Self::plain_with_suffix(contents, ".fastq")
    }

    fn plain_with_suffix(contents: &[u8], suffix: &str) -> Self {
        let path = Self::path_with_suffix(suffix);
        fs::write(&path, contents).unwrap();
        Self { path }
    }

    fn gzipped(contents: &[u8]) -> Self {
        Self::gzipped_with_suffix(contents, ".fastq.gz")
    }

    fn gzipped_with_suffix(contents: &[u8], suffix: &str) -> Self {
        let path = Self::path_with_suffix(suffix);
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        encoder.write_all(contents).unwrap();
        encoder.finish().unwrap();
        Self { path }
    }
}

fn run_fastq(fastq: &TempFastq, options: &[&str]) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_seqsource"));
    command.arg("fastq");
    command.args(options).arg(&fastq.path).output().unwrap()
}

fn records(bytes: &[u8]) -> (csv::StringRecord, Vec<csv::StringRecord>) {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(bytes);
    let header = reader.headers().unwrap().clone();
    let rows = reader.records().map(Result::unwrap).collect();
    (header, rows)
}

impl Drop for TempFastq {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[test]
fn default_output_is_headered_tsv_for_plain_and_gzipped_fastq() {
    let contents = b"@A00119:1:FLOWCELL:1:1101:1000:1000\nACGT\n+\n!!!!\n";
    for fastq in [TempFastq::plain(contents), TempFastq::gzipped(contents)] {
        let output = run_fastq(&fastq, &[]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let (header, rows) = records(&output.stdout);
        assert_eq!(
            header,
            csv::StringRecord::from(vec![
                "id",
                "fastq",
                "instrument_guess",
                "code",
                "fastq_header",
                "fastq_check_success",
                "failure_reason"
            ])
        );
        assert_eq!(rows.len(), 1);
        let basename = fastq.path.file_name().unwrap().to_str().unwrap();
        let suffix = if basename.ends_with(".gz") {
            ".fastq.gz"
        } else {
            ".fastq"
        };
        assert_eq!(rows[0].get(0), basename.strip_suffix(suffix));
        assert_eq!(rows[0].get(1), fastq.path.to_str());
        assert_eq!(rows[0].get(2), Some("Illumina NovaSeq 6000"));
        assert_eq!(rows[0].get(3), Some("A00119"));
        assert_eq!(rows[0].get(4), Some("@A00119:1:FLOWCELL:1:1101:1000:1000"));
        assert_eq!(rows[0].get(5), Some("true"));
        assert_eq!(rows[0].get(6), Some(""));
        assert!(output.stderr.is_empty());
    }
}

#[test]
fn simple_output_is_prediction_only() {
    let contents = b"@A00119:1:FLOWCELL:1:1101:1000:1000\nACGT\n+\n!!!!\n";
    for fastq in [TempFastq::plain(contents), TempFastq::gzipped(contents)] {
        let output = run_fastq(&fastq, &["--simple"]);
        assert!(output.status.success());
        assert_eq!(output.stdout, b"Illumina NovaSeq 6000\n");
    }
}

#[test]
fn unknown_prediction_is_still_a_success() {
    let contents = b"@ZZ12345:1:FLOWCELL\nACGT\n+\n!!!!\n";
    for fastq in [TempFastq::plain(contents), TempFastq::gzipped(contents)] {
        let output = run_fastq(&fastq, &[]);
        assert!(output.status.success());
        let (_, rows) = records(&output.stdout);
        assert_eq!(rows[0].get(2), Some("unknown"));
        assert_eq!(rows[0].get(3), Some("unknown"));
        assert_eq!(rows[0].get(4), Some("@ZZ12345:1:FLOWCELL"));
        assert_eq!(rows[0].get(5), Some("true"));
        assert_eq!(rows[0].get(6), Some(""));
        assert_eq!(run_fastq(&fastq, &["--simple"]).stdout, b"unknown\n");
    }
}

#[test]
fn id_override_and_header_are_quoted_when_needed() {
    let contents = b"@A00119:1:FLOWCELL\tsample \"one\"\nACGT\n+\n!!!!\n";
    let fastq = TempFastq::plain(contents);
    let output = run_fastq(&fastq, &["--id", "sample\t\"one\""]);
    assert!(output.status.success());
    let (_, rows) = records(&output.stdout);
    assert_eq!(rows[0].get(0), Some("sample\t\"one\""));
    assert_eq!(rows[0].get(4), Some("@A00119:1:FLOWCELL\tsample \"one\""));
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .contains("\"sample\t\"\"one\"\"\"")
    );
}

#[test]
fn derives_ids_from_recognized_suffixes_only() {
    let contents = b"@A00119:1\n";
    for (fastq, suffix) in [
        (TempFastq::plain_with_suffix(contents, ".fastq"), ".fastq"),
        (TempFastq::plain_with_suffix(contents, ".fq"), ".fq"),
        (
            TempFastq::gzipped_with_suffix(contents, ".fastq.gz"),
            ".fastq.gz",
        ),
        (TempFastq::gzipped_with_suffix(contents, ".fq.gz"), ".fq.gz"),
    ] {
        let output = run_fastq(&fastq, &[]);
        assert!(output.status.success());
        let (_, rows) = records(&output.stdout);
        let basename = fastq.path.file_name().unwrap().to_str().unwrap();
        assert_eq!(rows[0].get(0), basename.strip_suffix(suffix));
    }
    for suffix in [".reads", ".FASTQ"] {
        let fastq = TempFastq::plain_with_suffix(contents, suffix);
        let output = run_fastq(&fastq, &[]);
        assert!(output.status.success());
        let (_, rows) = records(&output.stdout);
        assert_eq!(rows[0].get(0), fastq.path.file_name().unwrap().to_str());
    }
}

#[test]
fn invalid_fastq_returns_error_without_tsv() {
    let fastq = TempFastq::plain(b"");
    let output = run_fastq(&fastq, &[]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("FASTQ file is empty"));

    let missing = TempFastq {
        path: TempFastq::path_with_suffix(".fastq"),
    };
    let output = run_fastq(&missing, &[]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("could not access FASTQ file"));
}

#[test]
fn id_and_simple_conflict_and_empty_id_is_rejected() {
    let fastq = TempFastq::plain(b"@A00119:1\n");
    for options in [&["--id", "sample", "--simple"][..], &["--id", ""][..]] {
        let output = run_fastq(&fastq, options);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
    }
    let old_option = run_fastq(&fastq, &["--detailed"]);
    assert!(!old_option.status.success());
}
