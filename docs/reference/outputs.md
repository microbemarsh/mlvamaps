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
| `locus_mapping_references.fasta` | Dominant observed representative amplicons. |
| `locus_read_alignments.sam` | minibwa locus-relative mappings. |
| `locus_mapping_summary.tsv` | Mapping rate, depth, coverage, and SNP totals by locus. |
| `locus_snps.tsv` | Filtered representative-relative SNP evidence. |
| `read_level_allele_predictions.tsv` | Per-read repeat-count probabilities and evidence weights. |
| `allele_calls.tsv` | Bayesian locus calls and posterior detail. |
| `amplirust/` | Native primer-pairing evidence. |
| `vsearch/` | Native unique sequences, centroids, and memberships. |
| `minibwa/` | Per-cluster membership alignments plus dominant-locus mapping FASTQ, references, and indexes. |

When `--no-locus-mapping` is used, the mapping summary and SNP tables are
header-only and the dominant-locus reference/SAM diagnostics are not produced.
Per-cluster minibwa alignments used for membership edit metrics are still run.

## Assembly outputs

| File | Meaning |
| --- | --- |
| `assembly_amplicons.tsv` | Accepted product coordinates, orientation, size, and primer mismatches. |
| `assembly_amplicons.fasta` | Extracted assembly primer products. |
| `read_support.tsv` | Optional mapped reads and mean coverage per product. |
| `read_support.sam` | minibwa alignments when `--reads` is used. |
| `minibwa/` | Assembly-product reference and minibwa indexes when `--reads` is used. |
| `amplirust/` | Native in-silico PCR evidence. |

## FASTQ call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Sufficient depth and a decisive in-range posterior. |
| `LOW_DEPTH` | Fewer reads than `--min-depth`. |
| `AMBIGUOUS` | Weak top posterior or insufficient separation from the second call. |
| `OUT_OF_RANGE` | Best repeat count lies outside the configured locus range. |
| `MULTIPLE_VARIANTS` | Several retained variants exist without a sufficiently dominant cluster. |
| `LOCUS_DROPOUT` | No retained read evidence produced a prediction. |

## Assembly call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Product found and repeat count calculated. |
| `PRESENT_COUNT_UNKNOWN` | Product found but panel metadata cannot convert its size to repeat count. |
| `NOT_FOUND` | No accepted paired-primer product found. |
