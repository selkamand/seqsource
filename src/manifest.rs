use std::{
    io::Write,
    num::NonZeroUsize,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail, ensure};
use indicatif::ProgressBar;
use rayon::{ThreadPoolBuilder, prelude::*};
use serde::Deserialize;

use crate::{
    analysis::{self, FailureKind, FastqAnalysis, FastqCheckError},
    tsv,
};

#[cfg(windows)]
const HOME_ENV: &str = "USERPROFILE";
#[cfg(not(windows))]
const HOME_ENV: &str = "HOME";

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
        let resolved_path = resolve_fastq_path(&row.fastq, directory).with_context(|| {
            format!(
                "cannot resolve FASTQ path {:?} in manifest record {} in {}",
                row.fastq,
                index + 2,
                manifest.display()
            )
        })?;
        entries.push(ManifestEntry {
            id: row.id,
            fastq: row.fastq,
            resolved_path,
        });
    }
    Ok(entries)
}

fn resolve_fastq_path(raw_path: &str, directory: &Path) -> Result<PathBuf> {
    let path = Path::new(raw_path);
    if let Ok(suffix) = path.strip_prefix("~") {
        let home = std::env::var_os(HOME_ENV)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("{HOME_ENV} is unset or empty"))?;
        let home = PathBuf::from(home);
        ensure!(home.is_absolute(), "{HOME_ENV} must be an absolute path");
        Ok(home.join(suffix))
    } else if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(directory.join(path))
    }
}

fn analyse_entry(entry: ManifestEntry) -> EntryOutcome {
    let check = analysis::analyse_fastq(&entry.resolved_path);
    EntryOutcome { entry, check }
}

pub fn analyse_from_manifest(
    manifest: &Path,
    threads: NonZeroUsize,
    output: &mut impl Write,
) -> Result<AnalysisSummary> {
    analyse_from_manifest_inner(manifest, threads, output, None)
}

/// Analyse manifest rows while advancing `progress` after each FASTQ check.
/// The total is set after the manifest has been validated.
pub fn analyse_from_manifest_with_progress(
    manifest: &Path,
    threads: NonZeroUsize,
    output: &mut impl Write,
    progress: &ProgressBar,
) -> Result<AnalysisSummary> {
    analyse_from_manifest_inner(manifest, threads, output, Some(progress))
}

fn analyse_from_manifest_inner(
    manifest: &Path,
    threads: NonZeroUsize,
    output: &mut impl Write,
    progress: Option<&ProgressBar>,
) -> Result<AnalysisSummary> {
    let entries = read_manifest(manifest)?;
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .build()
        .context("could not start FASTQ worker threads")?;
    if !entries.is_empty() {
        if let Some(progress) = progress {
            progress.set_length(entries.len() as u64);
            progress.tick();
        }
    }
    // Vec's indexed parallel iterator collects results in manifest order.
    let outcomes: Vec<_> = pool.install(|| {
        entries
            .into_par_iter()
            .map(|entry| {
                let outcome = analyse_entry(entry);
                if let Some(progress) = progress {
                    progress.inc(1);
                }
                outcome
            })
            .collect()
    });

    tsv::write_analysis_header(output)?;

    let mut summary = AnalysisSummary::default();
    for outcome in outcomes {
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
