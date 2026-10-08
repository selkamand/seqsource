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
        let path = Self::path_with_suffix(".fastq");
        fs::write(&path, contents).unwrap();
        Self { path }
    }

    fn gzipped(contents: &[u8]) -> Self {
        let path = Self::path_with_suffix(".fastq.gz");
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        encoder.write_all(contents).unwrap();
        encoder.finish().unwrap();
        Self { path }
    }
}

fn run_fastq(fastq: &TempFastq, detailed: bool) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_seqsource"));
    command.arg("fastq");
    if detailed {
        command.arg("--detailed");
    }
    let output = command.arg(&fastq.path).output().unwrap();
    assert!(
        output.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

impl Drop for TempFastq {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[test]
fn default_output_is_prediction_only_for_plain_and_gzipped_fastq() {
    let contents = b"@A00119:1:FLOWCELL:1:1101:1000:1000\nACGT\n+\n!!!!\n";
    for fastq in [TempFastq::plain(contents), TempFastq::gzipped(contents)] {
        assert_eq!(run_fastq(&fastq, false).stdout, b"Illumina NovaSeq 6000\n");
    }
}

#[test]
fn detailed_output_includes_code_and_header_for_plain_and_gzipped_fastq() {
    let contents = b"@A00119:1:FLOWCELL:1:1101:1000:1000\nACGT\n+\n!!!!\n";
    for fastq in [TempFastq::plain(contents), TempFastq::gzipped(contents)] {
        assert_eq!(
            run_fastq(&fastq, true).stdout,
            b"Illumina NovaSeq 6000\tA00119\t@A00119:1:FLOWCELL:1:1101:1000:1000\n"
        );
    }
}

#[test]
fn detailed_output_uses_unknown_placeholders_for_plain_and_gzipped_fastq() {
    let contents = b"@ZZ12345:1:FLOWCELL\nACGT\n+\n!!!!\n";
    for fastq in [TempFastq::plain(contents), TempFastq::gzipped(contents)] {
        assert_eq!(run_fastq(&fastq, false).stdout, b"unknown\n");
        assert_eq!(
            run_fastq(&fastq, true).stdout,
            b"unknown\tunknown\t@ZZ12345:1:FLOWCELL\n"
        );
    }
}

#[test]
fn detailed_output_quotes_header_with_tab_and_quote() {
    let contents = b"@A00119:1:FLOWCELL\tsample \"one\"\nACGT\n+\n!!!!\n";
    let fastq = TempFastq::plain(contents);
    assert_eq!(
        run_fastq(&fastq, true).stdout,
        b"Illumina NovaSeq 6000\tA00119\t\"@A00119:1:FLOWCELL\tsample \"\"one\"\"\"\n"
    );
}
