use std::io::{self, Write};

use crate::analysis::{FastqAnalysis, FastqCheckError};

const ANALYSIS_HEADER: [&str; 7] = [
    "id",
    "fastq",
    "instrument_guess",
    "code",
    "fastq_header",
    "fastq_check_success",
    "failure_reason",
];

pub fn write_analysis_header(writer: &mut impl Write) -> io::Result<()> {
    write_record(writer, &ANALYSIS_HEADER)
}

pub fn write_analysis_row(
    writer: &mut impl Write,
    id: &str,
    fastq: &str,
    check: &Result<FastqAnalysis, FastqCheckError>,
) -> io::Result<()> {
    let (guess, code, header, success, reason) = match check {
        Ok(analysis) => (
            analysis.instrument_guess.as_str(),
            analysis.code.as_str(),
            analysis.fastq_header.as_str(),
            "true",
            "",
        ),
        Err(error) => ("", "", "", "false", error.reason.as_str()),
    };
    write_record(writer, &[id, fastq, guess, code, header, success, reason])
}

/// Write a TSV record, quoting fields containing tabs, quotes, or line breaks.
pub fn write_record(writer: &mut impl Write, fields: &[&str]) -> io::Result<()> {
    for (index, field) in fields.iter().enumerate() {
        if index > 0 {
            writer.write_all(b"\t")?;
        }
        if field
            .chars()
            .any(|character| matches!(character, '\t' | '"' | '\r' | '\n'))
        {
            writer.write_all(b"\"")?;
            writer.write_all(field.replace('"', "\"\"").as_bytes())?;
            writer.write_all(b"\"")?;
        } else {
            writer.write_all(field.as_bytes())?;
        }
    }
    writer.write_all(b"\n")
}
