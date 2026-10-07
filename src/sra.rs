use std::{
    fs,
    io::{self, BufRead, Cursor},
    path::PathBuf,
};

/// Checks if fastq-dump is available on path
pub fn is_fastqdump_available() -> bool {
    // Returns Result<PathBuf, Error>
    which::which("fastq-dump").is_ok()
}

/// Download the first FASTQ read for an SRA accession to `<accession>.fastq`.
///
/// The output file is created in the current directory. `fastq-dump` must be
/// installed and available on `PATH`.
pub fn download_first_read(accession: &str) -> io::Result<PathBuf> {
    if accession.is_empty() || !accession.bytes().all(|byte| byte.is_ascii_alphanumeric()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "accession must contain only ASCII letters and digits",
        ));
    }

    // Check we have fastq-dump on path somewhere
    if !is_fastqdump_available() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "fastq-dump is not available on PATH; install the NCBI SRA Toolkit and add its bin directory to PATH",
        ));
    }

    let output = std::process::Command::new("fastq-dump")
        .arg(accession)
        .args([
            "-N",
            "1",
            "-X",
            "1",
            "--split-spot",
            "--skip-technical",
            "--origfmt",
            "--stdout",
        ])
        .output()?;

    if !output.status.success() {
        return Err(io::Error::other(format!(
            "fastq-dump failed for {accession} ({}): {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        )));
    }

    let first_read = first_fastq_record(&output.stdout)?;
    let output_path = PathBuf::from(format!("{accession}.fastq"));
    fs::write(&output_path, first_read)?;
    Ok(output_path)
}

fn first_fastq_record(stdout: &[u8]) -> io::Result<Vec<u8>> {
    let mut reader = Cursor::new(stdout);
    let mut record = Vec::new();

    for line_number in 0..4 {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 || !line.ends_with(b"\n") {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "fastq-dump did not return a complete FASTQ read",
            ));
        }
        if (line_number == 0 && !line.starts_with(b"@"))
            || (line_number == 2 && !line.starts_with(b"+"))
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "fastq-dump returned an invalid FASTQ read",
            ));
        }
        record.extend_from_slice(&line);
    }

    Ok(record)
}

#[cfg(test)]
mod tests {
    use std::io::ErrorKind;

    use super::*;

    #[test]
    fn keeps_only_the_first_fastq_read() {
        let stdout = b"@FIRST\nACGT\n+\n!!!!\n@SECOND\nTGCA\n+\n####\n";
        assert_eq!(
            first_fastq_record(stdout).unwrap(),
            b"@FIRST\nACGT\n+\n!!!!\n"
        );
    }

    #[test]
    fn rejects_incomplete_or_invalid_fastq_output() {
        for stdout in [
            b"@FIRST\nACGT\n+\n".as_slice(),
            b"FIRST\nACGT\n+\n!!!!\n",
            b"@FIRST\nACGT\nNOT_PLUS\n!!!!\n",
        ] {
            assert_eq!(
                first_fastq_record(stdout).unwrap_err().kind(),
                ErrorKind::InvalidData
            );
        }
    }

    #[test]
    fn rejects_invalid_accession_before_running_fastq_dump() {
        for accession in ["", "../SRR41014568", "-X"] {
            assert_eq!(
                download_first_read(accession).unwrap_err().kind(),
                ErrorKind::InvalidInput
            );
        }
    }
}
