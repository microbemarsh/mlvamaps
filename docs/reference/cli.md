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
| `--database` | None | Per-locus reference sequences for MAFFT alignment, RAxML-NG trees, and EPA-ng placement. |
| `--raxml-model` | `GTR+G` | RAxML-NG nucleotide model used for every locus tree. |

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

## Variant mixture estimation

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-mixture-fraction` | `0.01` | Minimum EM-estimated fraction for a variant to affect mixed-locus status and appear separately in the main abundance plot. |

Lower estimates remain in `vntr_mixture_abundance.tsv` as `TRACE` evidence and
are combined into one trace segment in the report.

## Representative mapping and SNP evidence

| Option | Default | Purpose |
| --- | --- | --- |
| `--minimap2-bin` | `minimap2` | Executable or path override. |
| `--no-locus-mapping` | Off | Skip dominant-locus minimap2 mapping and SNP evidence. |
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
| `--reads FASTQ` | Map reads to extracted products with minimap2. |
| `--bam BAM_OR_SAM` | Measure support from existing assembly alignments. |
| `--alignments BAM_OR_SAM` | Alias for `--bam`. |

`--reads` and `--bam` cannot be used together.

## External executable overrides

- `--amplirust-bin`
- `--vsearch-bin`
- `--minimap2-bin`
- `--mafft-bin`
- `--raxml-ng-bin`
- `--epa-ng-bin`

These options are useful for testing, containers, and installations whose
executables are not on the default `PATH`.
