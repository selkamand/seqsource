use std::{fmt, fs, io, path::Path};

use crate::{core, fastx};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqAnalysis {
    pub instrument_guess: String,
    pub code: String,
    pub fastq_header: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FailureKind {
    MissingFile,
    EmptyFastq,
    InvalidFastq,
    NonFilePath,
    OtherReadError,
}

#[derive(Debug)]
pub struct FastqCheckError {
    pub kind: FailureKind,
    pub reason: String,
}

impl fmt::Display for FastqCheckError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.reason.fmt(formatter)
    }
}

impl std::error::Error for FastqCheckError {}

/// Check only the first FASTQ header, then classify it using the shared patterns.
pub fn analyse_fastq(path: &Path) -> Result<FastqAnalysis, FastqCheckError> {
    let metadata = fs::metadata(path).map_err(|error| FastqCheckError {
        kind: if error.kind() == io::ErrorKind::NotFound {
            FailureKind::MissingFile
        } else {
            FailureKind::OtherReadError
        },
        reason: format!("could not access FASTQ file {}: {error}", path.display()),
    })?;

    if !metadata.is_file() {
        return Err(FastqCheckError {
            kind: FailureKind::NonFilePath,
            reason: format!("FASTQ path is not a regular file: {}", path.display()),
        });
    }

    let fastq_header = fastx::read_id_classified(path).map_err(|error| {
        let kind = match &error {
            fastx::HeaderReadError::Empty => FailureKind::EmptyFastq,
            fastx::HeaderReadError::InvalidHeader => FailureKind::InvalidFastq,
            fastx::HeaderReadError::Open(source) if source.kind() == io::ErrorKind::NotFound => {
                FailureKind::MissingFile
            }
            fastx::HeaderReadError::Open(_) => FailureKind::OtherReadError,
            fastx::HeaderReadError::Read(source)
                if path.extension().is_some_and(|extension| extension == "gz")
                    || matches!(
                        source.kind(),
                        io::ErrorKind::InvalidData
                            | io::ErrorKind::InvalidInput
                            | io::ErrorKind::UnexpectedEof
                    ) =>
            {
                FailureKind::InvalidFastq
            }
            fastx::HeaderReadError::Read(_) => FailureKind::OtherReadError,
        };
        FastqCheckError {
            kind,
            reason: format!(
                "could not read first FASTQ header from {}: {error}",
                path.display()
            ),
        }
    })?;

    let (instrument_guess, code) = core::identify_instrument_with_codes(&fastq_header);
    Ok(FastqAnalysis {
        instrument_guess,
        code,
        fastq_header,
    })
}
