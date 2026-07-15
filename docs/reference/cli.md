# CLI options and thresholds

Run `mlvamaps call --help` for the complete parser-generated reference.

## General controls

| Option | Default | Purpose |
| --- | --- | --- |
| `-o`, `--output`, `--outdir` | `results` | Output directory. |
| `--sample-id` | Input filename | Sample identifier. |
| `-t`, `--threads` | `32` | Worker/native threads; `0` uses all CPUs. |
| `--quiet` | Off | Suppress live progress. |
| `--max-primer-mismatches` | `3` | Allowed primer errors for paired-primer detection. |
| `--profiles` | None | Known MLVA profile TSV. |

## FASTQ filtering

| Option | Default |
| --- | --- |
| `--min-read-length` | `50` |
| `--max-read-length` | `100000` |
| `--min-qscore` | `0.0` |

## VSEARCH clustering

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-cluster-size` | `2` | Minimum reads in a retained variant. |
| `--cluster-min-identity` | `0.97` | Gap-aware global identity threshold. |
| `--vsearch-bin` | `vsearch` | Executable or path override. |

## Representative mapping and SNP evidence

| Option | Default | Purpose |
| --- | --- | --- |
| `--minibwa-bin` | `minibwa` | Executable or path override. |
| `--no-locus-mapping` | Off | Skip dominant-locus mapping and SNP evidence; per-cluster minibwa alignment remains required. |
| `--min-mapping-quality` | `0` | Minimum accepted primary-alignment MAPQ. |
| `--min-base-quality` | `20` | Minimum base quality used in depth and SNP evidence. |
| `--min-snp-depth` | `3` | Minimum quality-filtered depth at a position. |
| `--min-snp-alternate-reads` | `2` | Minimum reads supporting an alternate base. |
| `--min-snp-frequency` | `0.2` | Minimum alternate allele fraction. |

## Allele calling

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-depth` | `10` | Reads required to avoid `LOW_DEPTH`. |
| `--min-posterior` | `0.75` | Required top repeat-count posterior. |

## Assembly read support

| Option | Purpose |
| --- | --- |
| `--reads FASTQ` | Map accurate reads to extracted products with minibwa. |
| `--bam BAM_OR_SAM` | Measure support from existing assembly alignments. |
| `--alignments BAM_OR_SAM` | Alias for `--bam`. |

`--reads` and `--bam` cannot be used together.

## External executable overrides

- `--amplirust-bin`
- `--vsearch-bin`
- `--minibwa-bin`

These options are useful for testing, containers, and installations whose
executables are not on the default `PATH`.
