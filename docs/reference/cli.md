# CLI options and thresholds

Run `mlvamaps call --help` for the complete parser-generated reference.

## Completed-results MYOGA export

```bash
mlvamaps export-myoga --results results/ --metadata metadata.tsv \
  --metadata-id shared_identifier -o global_mlva/
```

| Option | Default | Purpose |
| --- | --- | --- |
| `--results DIR` | required | Recursively searched completed MLVAmaps result root. |
| `--metadata CSV_OR_TSV` | required | Metadata joined to recorded sample IDs. |
| `--metadata-id COLUMN` | `shared_identifier` | Metadata join column. |
| `--latitude COLUMN` | `latitude` | Latitude column; standard aliases are recognized for the default. |
| `--longitude COLUMN` | `longitude` | Longitude column; standard aliases are recognized for the default. |
| `--min-callable-fraction FRACTION` | `0` | Retain samples with any exact numeric repeat count; positive values require that fraction of assayed loci. |
| `--min-callable-loci COUNT` | `0` | Additional exact-call count; both callable thresholds apply. |
| `--min-pairwise-loci COUNT` | `1` | Required shared exact-call count for a supported pair. |
| `--min-pairwise-fraction FRACTION` | `0` | No fractional completeness requirement; pairs must still share at least `--min-pairwise-loci` exact calls. |
| `--distance {repeat,categorical}` | `repeat` | Metric used for the square matrix and relatedness tree; both remain in the pairwise table. |
| `--combined-markers` | off | Also recover accepted amplicons, infer per-locus repeat-masked SNP trees, and write a combined SNP/repeat tree. |
| `--loci TSV` | unset | Original rich locus panel used to mask retained amplicons when precomputed masked queries are absent. |
| `--phylogeny-snp-weight FLOAT` | `1.0` | Weight of normalized per-locus SNP-tree distance in the combined marker distance. |
| `--phylogeny-repeat-weight FLOAT` | `1.0` | Weight of normalized per-locus repeat distance in the combined marker distance. |
| `-t`, `--threads` | `32` | MAFFT/RAxML-NG CPU budget for combined-marker inference; `0` uses all CPUs. |
| `--mafft-bin` | `mafft` | MAFFT executable used for retrospective locus alignments. |
| `--raxml-ng-bin` | `raxml-ng` | RAxML-NG executable used when a locus has at least four SNP haplotypes. |
| `--raxml-model` | `DNA` | RAxML-NG model or model-selection set. |
| `-o`, `--output DIR` | required | Export directory. |
| `--force` | off | Replace existing export files. |

See [dataset aggregation and MYOGA export](../workflows/myoga-export.md) for
discovery rules, exact formulas, overlap pruning, and scientific caveats.

The primary input may be one supported FASTA/FASTQ file or a directory of such
files. Directory discovery is non-recursive, and each file is written beneath
`OUTDIR/<sample-id>/`; see [input formats](input-formats.md#input-directories).

## General controls

| Option | Default | Purpose |
| --- | --- | --- |
| `-p`, `--panel` | Required | Primer list or rich locus panel for the selected pipeline. |
| `-i`, `--input` | Required | Input FASTA/FASTQ path or directory; `sr` selects paired/single Illumina input. |
| `-o`, `--output`, `--outdir` | `results` | Output directory. |
| `--sample-id` | Input filename | Sample identifier. |
| `-t`, `--threads` | `32` | CPU budget shared across concurrent EPA-ng locus workers; `0` uses all CPUs. |
| `--quiet` | Off | Suppress live progress. |
| `--max-primer-mismatches` | `2` | Allowed primer errors for paired-primer detection, matching MLVA_finder's default. |
| `--profiles` | None | Known MLVA profile TSV. |
| `--database` | None | Reference-build directory whose fixed trees are reused, or a sequence-only database built on demand. |
| `--reference-metadata` | None | Reference date, coordinates, location, and source TSV/CSV; `reference_metadata.tsv` is auto-detected in database directories. |
| `--raxml-model` | `DNA` | Model-selection set when a sequence-only database requires new locus trees; ignored for reusable trees. |
| `--phylogeny-snp-weight` | `1.0` | Weight of normalized SNP-tree distance in combined marker ranking. |
| `--phylogeny-repeat-weight` | `1.0` | Weight of normalized tandem-repeat distance in combined marker ranking. |

## FASTQ filtering

### Illumina input and QC

| Option | Default | Meaning |
| --- | --- | --- |
| `-i sr` | Required | Select the dedicated Illumina short-read workflow. |
| `--fq1 FASTQ` | None | Mate 1 or single-end Illumina FASTQ. |
| `--fq2 FASTQ` | None | Mate 2; requires `--fq1`. |
| `--short-reads` | Off | Treat `-i DIRECTORY` as paired short reads and match exact `PREFIX_1.fastq.gz`/`PREFIX_2.fastq.gz` names. |
| `--read-technology` | `auto` | Compatibility override; `-i sr` selects Illumina automatically. |
| `--short-min-read-length` | `40` | Minimum post-trim length. |
| `--short-min-mean-quality` | `15` | Minimum mean Phred quality. |
| `--short-trim-quality` | `0` | Conservative 3-prime trim threshold; zero disables trimming. |
| `--short-min-pair-retention` | `0.5` | Fraction of a molecule's mates that must pass; a good mate may remain as an explicit orphan. |
| `--short-min-informative-molecules` | `3` | Boundary-informative molecules required to avoid `LOW_DEPTH`. |
| `--manifest` | None | Failure-isolated batch TSV. |
| `--sample-metadata` | None | CSV/TSV joined by sample ID. |
| `--force` | Off | Rerun already-successful manifest samples. |
| `--keep-intermediates` | Off | Retain per-locus SKESA FASTQs, contigs, and logs. |
| `--skesa-bin PATH` | `skesa` | Required native assembler executable for Illumina local assembly. |

Separate mates are supported. Interleaved data are not guessed or accepted.
Directory discovery is non-recursive, uses the shared prefix as `sample_id`,
and rejects any discovered prefix that is missing mate 1 or mate 2.

| Option | Default | Meaning |
| --- | --- | --- |
| `--taxon-screen-index DEACON_IDX` | None | Retain only reads matching a supplied target-taxon Deacon pangenome index before loading reads into mlvamaps. See [bede/deacon-indexes](https://github.com/bede/deacon-indexes) for build information. |
| `--taxon-screen-abs-threshold` | `1` | Required absolute number of shared minimizers. |
| `--taxon-screen-rel-threshold` | `0` | Required relative proportion of shared minimizers. |
| `--deacon-bin PATH` | `deacon` | Deacon executable. |

The taxon screen is opt-in. When enabled, Deacon receives the same `--threads`
budget as the rest of the pipeline and writes its retained FASTQ and JSON
summary under `taxon_screen/`.

| Option | Default |
| --- | --- |
| `--min-read-length` | `50` |
| `--max-read-length` | `100000` |
| `--min-qscore` | `15.0` |
| `--sample-mode` | `metagenome` |
| `--read-calling-convention` | `assembly` |
| `--fastq-strategy` | `recruit` |
| `--recruitment-preset` | None |
| `--recruitment-min-identity` | `0.9` |
| `--recruitment-min-aligned-bp` | `100` |
| `--recruitment-min-locus-margin` | `10` |
| `--recruitment-database` | None; falls back to `--database` |

`--fastq-strategy recruit` competitively maps retained reads to database or
synthetic locus products, records presence independently from genotype, and
uses primer pairing as a specificity fallback. Use `primer` for the historical
complete-primer-product workflow. `--recruitment-database` uses canonical
products from a reference build without also requesting phylogenetic
placement.

## Deprecated clustering compatibility

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-cluster-size` | `1` | Ignored compatibility option; mapping groups retain low-depth evidence. |
| `--cluster-min-identity` | `0.97` | Ignored compatibility option; FASTQ sequence clustering is no longer used. |
| `--vsearch-bin` | `vsearch` | Ignored compatibility option retained for older command lines. |

## Variant mixture estimation

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-mixture-fraction` | `0.01` | Minimum EM-estimated fraction for a variant to affect mixed-locus status and appear separately in the main abundance plot. |
| `--min-secondary-reads` | `2` | Reads required to promote an abundance-supported secondary from candidate to confirmed. |

Lower estimates remain in `vntr_mixture_abundance.tsv` as `TRACE` evidence and
are combined into one trace segment in the report.

## Representative mapping and SNP evidence

| Option | Default | Purpose |
| --- | --- | --- |
| `--minimap2-bin` | `minimap2` | Executable or path override. |
| `--no-locus-mapping` | Off | Skip downstream dominant-representative mapping and SNP evidence; default locus recruitment still runs unless `--fastq-strategy primer` is selected. |
| `--min-mapping-quality` | `0` | Minimum accepted primary-alignment MAPQ. |
| `--min-base-quality` | `20` | Minimum base quality used in depth and SNP evidence. |
| `--min-snp-depth` | `3` | Minimum quality-filtered depth at a position. |
| `--min-snp-alternate-reads` | `2` | Minimum reads supporting an alternate base. |
| `--min-snp-frequency` | `0.2` | Minimum alternate allele fraction. |

## Allele calling

Illumina mode uses transparent evidence classes rather than the long-read
posterior as its gate. Exact calls require two directly observed repeat
boundaries. Insert size may narrow an interval but never supplies an interval
midpoint as an exact call.

| Option | Default | Purpose |
| --- | --- | --- |
| `--min-depth` | `1` | Informative reads required to avoid `LOW_DEPTH`. One repeat-informative read is sufficient for a provisional call. |
| `--min-posterior` | `0.75` | Required top repeat-count probability for FASTQ and assembly calls. |
| `--repeat-range-tolerance` | `1.0` | Number of repeats allowed beyond either expected locus bound before assigning `OUT_OF_RANGE`. |
| `--max-confidence-depth` | `25` | Maximum effective dominant-cluster evidence used to sharpen FASTQ allele confidence. |
| `--debug-disagreements` | Off | Write read-level mapping-versus-anchor-measurement evidence and a locus-level disagreement summary. |

`--read-calling-convention assembly` makes primer-spanning reads use the same
product-size calibration and rounding as assembly products. The alternative
`probabilistic` convention evaluates the unrounded read measurement directly
on the half-unit grid.

The default `--sample-mode metagenome` flags any meaningful secondary allele
as `MULTIPLE_VARIANTS`; candidates below `--min-secondary-reads` remain visible
without changing the signature. Use `--sample-mode isolate` for cultured material to
retain the historical 80% dominance rule. In both modes the dominant allele
uses the assembly-equivalent convention, so it can be compared directly with
a later cultured assembly while posterior probability and dominant fraction
retain the original detection uncertainty.

## Assembly read support

| Option | Default | Purpose |
| --- | --- | --- |
| `--algorithm {legacy,novel}` | `legacy` | Assembly allele caller. `legacy` reproduces MLVA_finder product selection and rounding; `novel` uses the depth-aware probabilistic half-unit model. |
| `--assembly-round-tolerance FRACTION` | `0.25` | Integer tolerance for the legacy algorithm and compatibility CSVs. |
| `--reads FASTQ` | None | Map reads to extracted products with minimap2. |
| `--bam BAM_OR_SAM` | None | Measure support from existing assembly alignments. |
| `--alignments BAM_OR_SAM` | None | Alias for `--bam`. |

`--reads` and `--bam` cannot be used together.

## Validation command

```bash
mlvamaps validate \
  --truth assembly_results/calls.tsv \
  --long-read long_results/calls.tsv \
  --illumina short_results/calls.tsv \
  -o validation
```

The detailed table classifies exact matches, covered intervals, incorrect
calls, unresolved loci, false positives, and false negatives. The summary
reports exact-call accuracy, callable fraction, interval coverage, profile
recovery, and the false exact-call rate.

## External executable overrides

- `--amplirust-bin`
- `--skesa-bin`
- `--minimap2-bin`
- `--mafft-bin`
- `--raxml-ng-bin`
- `--epa-ng-bin`

These options are useful for testing, containers, and installations whose
executables are not on the default `PATH`.

## Reference database builder

`mlvamaps build-reference -i DIR -p PANEL.csv --metadata
metadata.csv -o OUT` extracts every locus from each metadata-linked assembly,
writes the per-locus database, and runs MAFFT plus RAxML-NG. The `-p` option
auto-detects minimal primer lists and rich locus panels.

| Option | Default | Purpose |
| --- | --- | --- |
| `--multiple-products` | `exclude` | Exclude, choose the best, or fail on equally good multiple products. |
| `--max-primer-mismatches` | `2` | Maximum Amplirust primer errors. |
| `--min-references-per-tree` | `3` | Minimum references required to infer a locus tree. |
| `-t`, `--threads` | `32` | Parallel assembly-extraction workers and maximum MAFFT/RAxML-NG threads; `0` uses all CPUs. RAxML-NG retries low-pattern loci with fewer threads automatically. |
| `--quiet` | off | Suppress per-assembly and per-locus progress updates. |
| `--raxml-model` | `DNA` | Model-selection set used independently for each reference tree. |

RAxML-NG 2.x selects a nucleotide model independently for each locus from the
default `DNA` model set. Pass an explicit model such as `--raxml-model GTR+G`
for reproducibility or compatibility with older RAxML-NG releases. Model
choice cannot resolve references whose retained SNP sequences are exactly
identical.
