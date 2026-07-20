# Output file reference

## Files returned for every call

| File | Meaning |
| --- | --- |
| `calls.tsv` | Compact per-locus result shared by FASTQ and assembly modes. |
| `mlva_fingerprint.tsv` | Wide sample-by-locus repeat-count fingerprint. |
| `mlva_fingerprint_probabilistic.tsv` | Long-form calls with confidence values. |
| `profile_matches.tsv` | Closest profile rows; header-only without a profile database. |
| `novelty_scores.tsv` | Nearest profile, score, and interpretation. |
| `report.html` | Self-contained human-readable report. |
| `locus_repeat_counts.tsv` | Exact individual-locus repeat counts in a compact long-form table. |

## Optional phylogenetic placement outputs

When `--database` is supplied, `phylogeny/` contains one directory per reference
locus. Each directory has the raw references, MAFFT reference alignment,
RAxML-NG reference search files and best-tree Newick, and—when a query locus
was callable—the query FASTA, MAFFT `--add --keeplength` alignment, query-only
aligned FASTA, and EPA-ng `epa_result.jplace`. The default RAxML-NG model is
`GTR+G` and can be changed with `--raxml-model`. EPA-ng consumes the optimized
RAxML-NG `.bestModel` file and places the query without changing the reference
topology. Summed matching distances combine the placement pendant/distal lengths
with the RAxML-NG reference-tree branch lengths.

| File | Meaning |
| --- | --- |
| `phylogeny/locus_status.tsv` | Whether each database locus received a query placement. |
| `phylogeny/locus_phylogenetic_distances.tsv` | Query-to-reference patristic distance plus EPA-ng edge, likelihood weight, pendant length, and distal length. |
| `phylogeny/phylogenetic_matches.tsv` | References ranked by the sum of distances across all placed loci. |

## Compact call columns

`calls.tsv` contains:

```text
sample_id
locus_id
present
repeat_count
repeat_count_raw
product_size_bp
read_depth
mean_coverage
status
evidence
```

Fields that do not apply to an input mode remain blank. A FASTQ call normally
has read depth but no assembly product size. An assembly-only call has product
size but no read depth unless support data are supplied.

## FASTQ outputs

| File | Meaning |
| --- | --- |
| `qc_summary.tsv` | Input and read-filtering totals. |
| `filtered_reads.fastq.gz` | Reads retained after length and quality filtering. |
| `filtered_reads.fasta` | Lossless sequence projection used by Amplirust. |
| `read_locus_assignments.tsv` | Primer-supported locus, orientation, score, and assignment QC. |
| `read_repeat_features.tsv` | Repeat coordinates, counts, patterns, motifs, and quality per read. |
| `vntr_asv_table.tsv` | Retained VSEARCH variants and representative statistics. |
| `vntr_asv_memberships.tsv` | Cluster membership and per-read substitution/indel evidence. |
| `vntr_asv_representatives.fasta` | Observed representative repeat sequences. |
| `vntr_mixture_abundance.tsv` | Emu-inspired EM abundance, estimated read count, and meaningful/trace classification for every retained variant. |
| `locus_mapping_references.fasta` | Dominant observed representative amplicons. |
| `locus_read_alignments.sam` | minimap2 locus-relative mappings. |
| `locus_mapping_summary.tsv` | Mapping rate, depth, coverage, and SNP totals by locus. |
| `locus_snps.tsv` | Filtered representative-relative SNP evidence. |
| `read_level_allele_predictions.tsv` | Per-read repeat-count probabilities and evidence weights. |
| `allele_calls.tsv` | Bayesian locus calls, posterior detail, meaningful-variant count, and dominant mixture fraction. |
| `amplirust/` | Native primer-pairing evidence. |
| `vsearch/` | Native unique sequences, centroids, and memberships. |
| `minimap2/` | Dominant-locus mapping FASTQ and reference diagnostics. |

When `--no-locus-mapping` is used, the mapping summary and SNP tables are
header-only and the dominant-locus reference/SAM diagnostics are not produced.
Parasail global alignments used for membership edit metrics still run because
they are independent of locus-wide read mapping.

## Assembly outputs

| File | Meaning |
| --- | --- |
| `assembly_amplicons.tsv` | Accepted product coordinates, orientation, size, and primer mismatches. |
| `assembly_amplicons.fasta` | Extracted assembly primer products. |
| `read_support.tsv` | Optional mapped reads and mean coverage per product. |
| `read_support.sam` | minimap2 alignments when `--reads` is used. |
| `amplirust/` | Native in-silico PCR evidence. |

## FASTQ call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Sufficient depth and a decisive in-range posterior. |
| `LOW_DEPTH` | Fewer reads than `--min-depth`. |
| `AMBIGUOUS` | Weak top posterior or insufficient separation from the second call. |
| `OUT_OF_RANGE` | Best repeat count lies outside the configured locus range. |
| `MULTIPLE_VARIANTS` | Several EM-meaningful variants remain and the dominant estimated fraction is below 0.8. |
| `LOCUS_DROPOUT` | No retained read evidence produced a prediction. |

## Assembly call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Product found and repeat count calculated. |
| `PRESENT_COUNT_UNKNOWN` | Product found but panel metadata cannot convert its size to repeat count. |
| `NOT_FOUND` | No accepted paired-primer product found. |
