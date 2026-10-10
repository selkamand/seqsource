use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

use flate2::read::MultiGzDecoder;

/// Read the first FASTQ header without its line ending, decoding `.gz` files.
pub fn read_id(fastq: &Path) -> io::Result<String> {
    read_id_classified(fastq).map_err(HeaderReadError::into_io_error)
}

#[derive(Debug)]
pub(crate) enum HeaderReadError {
    Empty,
    InvalidHeader,
    Open(io::Error),
    Read(io::Error),
}

impl HeaderReadError {
    fn into_io_error(self) -> io::Error {
        match self {
            Self::Empty => io::Error::new(io::ErrorKind::InvalidData, "FASTQ file is empty"),
            Self::InvalidHeader => io::Error::new(
                io::ErrorKind::InvalidData,
                "first FASTQ line must start with '@'",
            ),
            Self::Open(error) | Self::Read(error) => error,
        }
    }
}

impl std::fmt::Display for HeaderReadError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Empty => formatter.write_str("FASTQ file is empty"),
            Self::InvalidHeader => formatter.write_str("first FASTQ line must start with '@'"),
            Self::Open(error) | Self::Read(error) => error.fmt(formatter),
        }
    }
}

pub(crate) fn read_id_classified(fastq: &Path) -> Result<String, HeaderReadError> {
    let file = File::open(fastq).map_err(HeaderReadError::Open)?;
    if fastq.extension().is_some_and(|extension| extension == "gz") {
        read_header(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        read_header(BufReader::new(file))
    }
}

fn read_header(mut reader: impl BufRead) -> Result<String, HeaderReadError> {
    let mut header = String::new();
    if reader
        .read_line(&mut header)
        .map_err(HeaderReadError::Read)?
        == 0
    {
        return Err(HeaderReadError::Empty);
    }

    if header.ends_with('\n') {
        header.pop();
        if header.ends_with('\r') {
            header.pop();
        }
    }

    if !header.starts_with('@') {
        return Err(HeaderReadError::InvalidHeader);
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        io::{ErrorKind, Write},
        path::PathBuf,
        sync::atomic::{AtomicU64, Ordering},
        time::{SystemTime, UNIX_EPOCH},
    };

    use flate2::{Compression, write::GzEncoder};

    use super::read_id;
    use crate::core::identify_instrument;

    static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

    struct TempFastq {
        path: PathBuf,
    }

    impl TempFastq {
        fn new(contents: &[u8]) -> Self {
            Self::with_suffix(contents, ".fastq")
        }

        fn with_suffix(contents: &[u8], suffix: &str) -> Self {
            let timestamp = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "seqsource-fastx-{}-{timestamp}-{}{suffix}",
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

    fn gzip(contents: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(contents).unwrap();
        encoder.finish().unwrap()
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

    #[test]
    fn reads_gzipped_fastq_headers() {
        for suffix in [".fastq.gz", ".fq.gz"] {
            let fastq = TempFastq::with_suffix(&gzip(b"@ST-E00185:1\r\nACGT\r\n"), suffix);
            let header = read_id(&fastq.path).unwrap();
            assert_eq!(header, "@ST-E00185:1");
            assert_eq!(identify_instrument(&header), "Illumina HiSeq X");
        }
    }

    #[test]
    fn rejects_gzipped_empty_and_non_header_first_lines() {
        for contents in [b"".as_slice(), b"A00119:1\n"] {
            let fastq = TempFastq::with_suffix(&gzip(contents), ".fastq.gz");
            assert_eq!(
                read_id(&fastq.path).unwrap_err().kind(),
                ErrorKind::InvalidData
            );
        }
    }

    #[test]
    fn returns_an_error_for_invalid_gzip_data() {
        let fastq = TempFastq::with_suffix(b"not gzip data", ".fastq.gz");
        assert!(read_id(&fastq.path).is_err());
    }
}
