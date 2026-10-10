# SeqSource

> [!WARNING]
> This tool is in early development and is not yet ready for use.

Rapidly infer instrument used to generate whole-genome sequencing dataset based on fastq.

## Install

```
cargo install --git https://github.com/selkamand/seqsource
```

## Quick Start

Run on a single FASTQ file (plain or gzip-compressed):

```
# Run on a single fastq
seqsource fastq path/to/sequence.fastq

# Same as above but manually specify id to use for output
seqsource fastq --id sample_1 path/to/sequence.fastq.gz
```

Or check many files based on a manifest with `id` and `fastq` columns:

```tsv
id	fastq
sample_1	reads/sample_1.fastq.gz
sample_2	reads/sample_2.fastq.gz
```

```sh
seqsource manifest path/to/manifest.tsv
```

## Output

```tsv
id	fastq	instrument_guess	code	fastq_header	fastq_check_success	failure_reason
```
