use std::{
    io::Write,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use serde::Deserialize;

use crate::{
    analysis::{self, FailureKind, FastqAnalysis, FastqCheckError},
    tsv,
};

#[derive(Debug, Deserialize)]
struct ManifestRow {
    id: String,
    fastq: String,
}

/// Each entry contains everything needed for an independent FASTQ check.
#[derive(Debug)]
struct ManifestEntry {
    id: String,
    fastq: String,
    resolved_path: PathBuf,
}

struct EntryOutcome {
    entry: ManifestEntry,
    check: Result<FastqAnalysis, FastqCheckError>,
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct AnalysisSummary {
    pub total: usize,
    pub successful: usize,
    pub missing_file: usize,
    pub empty_fastq: usize,
    pub invalid_fastq: usize,
    pub non_file_path: usize,
    pub other_read_error: usize,
}

impl AnalysisSummary {
    pub fn failed(&self) -> usize {
        self.total - self.successful
    }

    fn record(&mut self, outcome: &EntryOutcome) {
        self.total += 1;
        match &outcome.check {
            Ok(_) => self.successful += 1,
            Err(error) => match error.kind {
                FailureKind::MissingFile => self.missing_file += 1,
                FailureKind::EmptyFastq => self.empty_fastq += 1,
                FailureKind::InvalidFastq => self.invalid_fastq += 1,
                FailureKind::NonFilePath => self.non_file_path += 1,
                FailureKind::OtherReadError => self.other_read_error += 1,
            },
        }
    }
}

/// Validate every record before writing any output.
fn read_manifest(manifest: &Path) -> Result<Vec<ManifestEntry>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(manifest)
        .with_context(|| format!("could not open manifest {}", manifest.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("could not read manifest header from {}", manifest.display()))?;
    for required in ["id", "fastq"] {
        if headers.iter().filter(|header| *header == required).count() != 1 {
            bail!(
                "manifest {} must contain exactly one '{required}' column",
                manifest.display()
            );
        }
    }

    let directory = manifest.parent().unwrap_or_else(|| Path::new("."));
    let mut entries = Vec::new();
    for (index, row) in reader.deserialize::<ManifestRow>().enumerate() {
        let row = row.with_context(|| {
            format!(
                "invalid manifest record {} in {}",
                index + 2,
                manifest.display()
            )
        })?;
        if row.id.trim().is_empty() || row.fastq.trim().is_empty() {
            bail!(
                "manifest record {} in {} must have nonempty id and fastq fields",
                index + 2,
                manifest.display()
            );
        }
        let path = Path::new(&row.fastq);
        let resolved_path = if path.is_absolute() {
            path.to_path_buf()
        } else {
            directory.join(path)
        };
        entries.push(ManifestEntry {
            id: row.id,
            fastq: row.fastq,
            resolved_path,
        });
    }
    Ok(entries)
}

fn analyse_entry(entry: ManifestEntry) -> EntryOutcome {
    let check = analysis::analyse_fastq(&entry.resolved_path);
    EntryOutcome { entry, check }
}

pub fn analyse_from_manifest(manifest: &Path, output: &mut impl Write) -> Result<AnalysisSummary> {
    let entries = read_manifest(manifest)?;
    tsv::write_analysis_header(output)?;

    let mut summary = AnalysisSummary::default();
    for entry in entries {
        let outcome = analyse_entry(entry);
        summary.record(&outcome);
        tsv::write_analysis_row(
            output,
            &outcome.entry.id,
            &outcome.entry.fastq,
            &outcome.check,
        )?;
    }
    Ok(summary)
}
