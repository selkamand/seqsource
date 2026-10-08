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
        /// Print prediction, instrument code, and FASTQ header as TSV
        #[arg(long)]
        detailed: bool,
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
        Commands::Fastq { fastq, detailed } => {
            let metadata = fs::metadata(&fastq)
                .with_context(|| format!("could not access FASTQ file {}", fastq.display()))?;
            ensure!(
                metadata.is_file(),
                "FASTQ path is not a regular file: {}",
                fastq.display()
            );

            let fastq_header = seqsource::fastx::read_id(&fastq)?;
            if detailed {
                let (prediction, code) =
                    seqsource::core::identify_instrument_with_codes(&fastq_header);
                println!(
                    "{}\t{}\t{}",
                    tsv_field(&prediction),
                    tsv_field(&code),
                    tsv_field(&fastq_header)
                );
            } else {
                println!("{}", seqsource::core::identify_instrument(&fastq_header));
            }
        }
        Commands::Download { sra } => {
            eprintln!("Downloading first read for SRA acession: {sra}");
            let filepath = seqsource::sra::download_first_read(&sra)?;
            eprintln!("Successfully downloaded read to {}", filepath.display());
        }
    };
    Ok(())
}

fn tsv_field(value: &str) -> String {
    if value
        .chars()
        .any(|character| matches!(character, '\t' | '"' | '\r' | '\n'))
    {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}
