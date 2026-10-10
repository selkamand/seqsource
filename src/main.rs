use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};
use std::{
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Analyse a single FASTQ file
    Fastq {
        /// FASTQ file path
        fastq: PathBuf,
        /// Print only the instrument guess
        #[arg(long)]
        simple: bool,
        /// Override the ID derived from the FASTQ filename
        #[arg(long, conflicts_with = "simple")]
        id: Option<String>,
    },

    /// Download the first read from an SRA accession
    Download {
        /// SRA accession code
        sra: String,
    },
    /// Analyse FASTQ paths from a TSV manifest
    Manifest {
        /// TSV file with id and fastq columns
        #[arg(value_name = "MANIFEST.TSV")]
        manifest: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Fastq { fastq, simple, id } => {
            let analysis = seqsource::analysis::analyse_fastq(&fastq)?;
            if simple {
                println!("{}", analysis.instrument_guess);
            } else {
                let id = match id {
                    Some(id) => {
                        ensure!(!id.trim().is_empty(), "--id must not be empty");
                        id
                    }
                    None => derive_id(&fastq)?,
                };
                let path = fastq.to_str().context("FASTQ path is not valid UTF-8")?;
                let mut output = BufWriter::new(io::stdout().lock());
                seqsource::tsv::write_analysis_header(&mut output)?;
                seqsource::tsv::write_analysis_row(&mut output, &id, path, &Ok(analysis))?;
                output.flush()?;
            }
        }
        Commands::Download { sra } => {
            eprintln!("Downloading first read for SRA acession: {sra}");
            let filepath = seqsource::sra::download_first_read(&sra)?;
            eprintln!("Successfully downloaded read to {}", filepath.display());
        }
        Commands::Manifest { manifest } => {
            let started = Instant::now();
            let mut output = BufWriter::new(io::stdout().lock());
            let summary = seqsource::manifest::analyse_from_manifest(&manifest, &mut output)?;
            output.flush()?;
            eprintln!(
                "Samples: total={} successful={} failed={}",
                summary.total,
                summary.successful,
                summary.failed()
            );
            eprintln!(
                "Failures: missing_file={} empty_fastq={} invalid_fastq={} non_file_path={} other_read_error={}",
                summary.missing_file,
                summary.empty_fastq,
                summary.invalid_fastq,
                summary.non_file_path,
                summary.other_read_error
            );
            eprintln!("Elapsed: {:.3}s", started.elapsed().as_secs_f64());
        }
    };
    Ok(())
}

fn derive_id(fastq: &Path) -> Result<String> {
    let basename = fastq
        .file_name()
        .and_then(|name| name.to_str())
        .context("cannot derive ID from FASTQ filename; supply --id")?;
    let id = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
        .iter()
        .find_map(|suffix| basename.strip_suffix(suffix))
        .unwrap_or(basename);
    ensure!(
        !id.trim().is_empty(),
        "derived FASTQ ID is empty; supply --id"
    );
    Ok(id.to_string())
}
