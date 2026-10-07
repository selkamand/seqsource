# SeqSource

Rapidly infer instrument used to generate whole-genome sequencing dataset based on fastq.

## Install

```
cargo install --git https://github.com/selkamand/seqsource
```

## Quick Start

Run on a single FASTQ file, plain or gzip-compressed

```
seqsource fastq path/to/sequence.fastq
seqsource fastq path/to/sequence.fastq.gz
```

Run on many samples using a manifest tsv with the following columns 

1. id: A unique identifier for that sample
2. fastq: Path to a fastq file

```
seqsource --manifest manifest.tsv
```

## Output

A tsv is printed to stdout with the following columns

1. id
2. instrument_code
3. inferred_instrument
