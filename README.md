# SeqSource

> [!WARNING]
> This tool is in early development and is not yet ready for use.

Rapidly infer instrument used to generate whole-genome sequencing dataset based on fastq.

## Install

```
cargo install --git https://github.com/selkamand/seqsource
```

## Quick Start

Run on a single FASTQ file (plain or gzip-compressed)

```
seqsource fastq path/to/sequence.fastq
seqsource fastq --detailed path/to/sequence.fastq.gz
```

## Output

By default, `seqsource fastq` prints the predicted instrument name on one line:

```
Illumina NovaSeq 6000
```

With `--detailed`, it prints one headerless TSV row with three columns: prediction,
instrument code, and the full first FASTQ header (including `@`):

```
Illumina HiSeq X	ST-E00545	@ST-E00545:419:HLL3MCCXY:8:1101:5041:1713
```

When there is no match, the first two columns are `unknown`. If multiple patterns
match, the prediction lists their names separated by semicolons after `ambiguous:`,
and the codes are separated by colons in the same order. Fields containing tabs or
quotes are quoted, with embedded quotes doubled.
