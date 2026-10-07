use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};
use std::{fs, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Validate a single FASTQ file
    Fastq {
        /// FASTQ file path
        fastq: PathBuf,
    },

    /// Download the first read from an SRA accession
    Download {
        /// SRA accession code
        sra: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Fastq { fastq } => {
            let metadata = fs::metadata(&fastq)
                .with_context(|| format!("could not access FASTQ file {}", fastq.display()))?;
            ensure!(
                metadata.is_file(),
                "FASTQ path is not a regular file: {}",
                fastq.display()
            );

            let fastq_header = seqsource::fastx::read_id(&fastq)?;
            let instrument = seqsource::core::identify_instrument(&fastq_header);
            println!("{}", instrument);
        }
        Commands::Download { sra } => {
            eprintln!("Downloading first read for SRA acession: {sra}");
            let filepath = seqsource::sra::download_first_read(&sra)?;
            eprintln!("Successfully downloaded read to {}", filepath.display());
        }
    };
    Ok(())
}
