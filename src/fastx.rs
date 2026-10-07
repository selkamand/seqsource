use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::PathBuf,
};

/// Read the first FASTQ header without its line ending.
pub fn read_id(fastq: &PathBuf) -> io::Result<String> {
    let mut header = String::new();
    if BufReader::new(File::open(fastq)?).read_line(&mut header)? == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTQ file is empty",
        ));
    }

    if header.ends_with('\n') {
        header.pop();
        if header.ends_with('\r') {
            header.pop();
        }
    }

    if !header.starts_with('@') {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "first FASTQ line must start with '@'",
        ));
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        io::ErrorKind,
        path::PathBuf,
        sync::atomic::{AtomicU64, Ordering},
        time::{SystemTime, UNIX_EPOCH},
    };

    use super::read_id;
    use crate::core::identify_instrument;

    static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

    struct TempFastq {
        path: PathBuf,
    }

    impl TempFastq {
        fn new(contents: &[u8]) -> Self {
            let timestamp = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "seqsource-fastx-{}-{timestamp}-{}.fastq",
                std::process::id(),
                NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed)
            ));
            fs::write(&path, contents).unwrap();
            Self { path }
        }
    }

    impl Drop for TempFastq {
        fn drop(&mut self) {
            let _ = fs::remove_file(&self.path);
        }
    }

    #[test]
    fn reads_first_header_with_lf_crlf_or_no_line_ending() {
        for (contents, expected) in [
            (b"@A00119:1 \nACGT\n+\n!!!!\n".as_slice(), "@A00119:1 "),
            (b"@A00119:1\r\nACGT\r\n".as_slice(), "@A00119:1"),
            (b"@A00119:1".as_slice(), "@A00119:1"),
        ] {
            let fastq = TempFastq::new(contents);
            let header = read_id(&fastq.path).unwrap();
            assert_eq!(header, expected);
            assert_eq!(identify_instrument(&header), "Illumina NovaSeq 6000");
        }
    }

    #[test]
    fn rejects_empty_files_and_non_header_first_lines() {
        for contents in [b"".as_slice(), b"\n", b"A00119:1\nACGT\n"] {
            let fastq = TempFastq::new(contents);
            assert_eq!(
                read_id(&fastq.path).unwrap_err().kind(),
                ErrorKind::InvalidData
            );
        }
    }
}
