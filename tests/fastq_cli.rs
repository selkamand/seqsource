use std::{
    fs::{self, File},
    io::Write,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicU64, Ordering},
    time::{SystemTime, UNIX_EPOCH},
};

use flate2::{Compression, write::GzEncoder};

static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

struct TempFastq {
    path: PathBuf,
}

impl TempFastq {
    fn gzipped(contents: &[u8]) -> Self {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "seqsource-cli-{}-{timestamp}-{}.fastq.gz",
            std::process::id(),
            NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed)
        ));
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        encoder.write_all(contents).unwrap();
        encoder.finish().unwrap();
        Self { path }
    }
}

impl Drop for TempFastq {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[test]
fn identifies_instrument_from_gzipped_fastq() {
    let fastq = TempFastq::gzipped(b"@A00119:1:FLOWCELL:1:1101:1000:1000\nACGT\n+\n!!!!\n");
    let output = Command::new(env!("CARGO_BIN_EXE_seqsource"))
        .arg("fastq")
        .arg(&fastq.path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"Illumina NovaSeq 6000\n");
}
