# Output file reference

## Files returned for every call

| File | Meaning |
| --- | --- |
| `calls.tsv` | Compact per-locus result shared by FASTQ and assembly modes. |
| `mlva_fingerprint.tsv` | Wide sample-by-locus repeat-count fingerprint. |
| `mlva_fingerprint_probabilistic.tsv` | Long-form calls with confidence values. |
| `profile_matches.tsv` | Closest profile rows; header-only without a profile database. |
| `report.html` | Self-contained interpretation report with conditional sample findings, locus-quality flags, gel evidence, profile matches, closest reference genomes, and expandable technical tables. |
| `locus_repeat_counts.tsv` | Exact individual-locus repeat counts in a compact long-form table. |
| `allele_probability_distribution.tsv` | Ranked integer/half-unit allele probabilities, selected state, and inference method. |

## Reference builder outputs

| File | Meaning |
| --- | --- |
| `database/LOCUS.fasta` | Unmasked extracted reference amplicons in the format accepted by `--database`. |
| `database/reference_metadata.tsv` | Metadata normalized to a `reference_id` key. |
| `database/reference_sequence_index.tsv` | Canonical amplicon, repeat-masked SNP, and complete marker SHA-256 keys used by the default exact-match fast path. |
| `database/reference_assemblies.tsv` | Reference-ID, source-assembly path, and canonical whole-genome SHA-256 used for assembly-query tie breaking. |
| `reference_build_manifest.tsv` | Per-reference/locus product counts, selected product, primer errors, and exclusion status. |
| `phylogeny/LOCUS.tree` | Portable Newick SNP tree for the locus. |
| `phylogeny/LOCUS/references.aligned.fasta` | Repeat-masked MAFFT alignment used for tree inference. |
| `phylogeny/LOCUS/reference.mlvamaps.raxml.log` | Full RAxML-NG output for every attempted thread count. |
| `phylogeny/reference_tree_status.tsv` | Tree completion or insufficient-reference status for every panel locus. |
| `phylogeny/reference_marker_components.tsv` | Retained repeat measurements and the masking method for each reference marker. |
| `phylogeny/reference_haplotype_groups.tsv` | Mapping from every reference ID to the representative repeat-masked SNP haplotype used as a tree tip. Identical SNP haplotypes are collapsed when at least two distinct haplotypes remain. |
| `myoga_metadata.csv` | Metadata with a `genome_id` column matching Newick tip labels. |

## Optional phylogenetic placement outputs

When `--database` is supplied, `phylogeny/` contains one directory per reference
locus. Each directory has the raw references, fixed MAFFT reference alignment,
RAxML-NG tree/model artifacts, and—when a query locus
was callable—the query FASTA, MAFFT `--add --keeplength` alignment, query-only
aligned FASTA, and EPA-ng `epa_result.jplace`. Artifacts from a reference-build
database are copied without rerunning RAxML-NG. For sequence-only databases
that lack reusable trees, RAxML-NG selects a
model independently for each locus from its `DNA` set; this can be changed with
`--raxml-model`. EPA-ng consumes the optimized
RAxML-NG `.bestModel` file and places the query without changing the reference
topology. Callable loci are placed concurrently within the `--threads` CPU
budget. Post-placement patristic matrices, cross-locus SNP/repeat aggregation,
and neighbor joining use NumPy's compiled C/BLAS kernels rather than Python
reference-pair loops. Summed matching distances combine the placement
pendant/distal lengths
with the RAxML-NG reference-tree branch lengths. The output retains both the
highest-likelihood-weight placement distance and the expected distance across
all candidate placements, weighted by normalized likelihood weight ratios.
Ranking uses the likelihood-weighted sum and reports its gap to the next
reference.

Before launching alignment or placement, MLVAMaps checks the persistent sequence
index. A query matching the same reference or tied reference group at every
configured panel locus is reported with zero distance and status
`EXACT_AMPLICON_MATCH` or `EXACT_MARKER_MATCH`; MAFFT, RAxML-NG, and EPA-ng are
skipped for that query. Every primer set must be callable in the query and
represented in the index. A missing or dropped-out locus disables the exact
fast path, and the query continues through the alignment and EPA-ng workflow.

| File | Meaning |
| --- | --- |
| `phylogeny/locus_status.tsv` | Whether each database locus received a query placement. |
| `phylogeny/locus_phylogenetic_distances.tsv` | Best-placement and likelihood-weighted query-to-reference patristic distances, placement entropy, EPA-ng edge, likelihood weight, pendant length, and distal length. |
| `phylogeny/phylogenetic_matches.tsv` | Complete references ranked by likelihood-weighted summed distance across all placed loci, with raw best-placement sums and rank gaps. |
| `phylogeny/marker_components.tsv` | Query and reference repeat counts, repeat-unit haplotypes, masking coordinates, and SNP-sequence lengths. |
| `phylogeny/locus_marker_distances.tsv` | Per-locus EPA/tree SNP distance, direct aligned SNP divergence, exact-match status, hybrid normalized SNP distance, and explicit repeat-count distance for every reference. |
| `phylogeny/combined_marker_matches.tsv` | References ranked by the configurable weighted sum of identity-aware hybrid SNP and repeat distances. Exact assembly-query ties are secondarily ranked by canonical genome identity, MUMmer4 `dnadiff` SNPs, indel bases, and one-to-one aligned fraction without changing the marker distance. |
| `phylogeny/whole_genome_dnadiff.tsv` | Interpreted whole-genome identity, SNP, indel, and alignment-coverage results for tied exact references. Individual MUMmer reports are retained under `phylogeny/dnadiff/`. |
| `phylogeny/closest_reference_bands.tsv` | Exact per-locus amplicon sizes and repeat calls for the top combined-marker reference (or top SNP-tree match when repeat-aware ranking is unavailable), used for the reference lane in the generated gel. |
| `phylogeny/combined_markers.tree` | MYOGA-compatible Newick neighbor-joining tree inferred from the combined normalized SNP-plus-repeat distance matrix. Tip labels are reference IDs plus the query sample ID. |

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
primary_read_depth
mean_coverage
allele_confidence
second_best_repeat_count
second_best_probability
inference_method
dominant_variant_fraction
num_candidate_variants
num_confirmed_secondary_variants
secondary_alleles
allele_distribution
status
evidence
```

Fields that do not apply to an input mode remain blank. A FASTQ call based on a
complete dominant local product includes its `product_size_bp` and unrounded
`repeat_count_raw`; a provisional partial-read call may leave product size
blank. An assembly-only call has product size but no read depth unless support
data are supplied.

## FASTQ outputs

| File | Meaning |
| --- | --- |
| `qc_summary.tsv` | Input and read-filtering totals. |
| `filtered_reads.fastq.gz` | Reads retained after length and quality filtering. |
| `filtered_reads.fasta` | Lossless sequence projection used by Amplirust. |
| `locus_recruited_reads.tsv` | Competitive per-read locus mapping, candidate allele, alignment quality, and presence/genotype evidence class. |
| `locus_presence.tsv` | Per-locus mapped, full-product, and repeat-informative read counts with presence status. |
| `local_locus_products.fasta` | Modal-length majority products reconstructed from the dominant read cluster and used for assembly-equivalent primary repeat counts. |
| `recruitment/locus_recruitment_references.fasta` | Database-derived or synthetic competitive locus/allele reference bank. |
| `recruitment/read_recruitment.sam` | Raw competitive long-read mappings used for presence and provisional genotype evidence. |
| `read_locus_assignments.tsv` | Primer-supported locus, orientation, score, and assignment QC. |
| `read_repeat_features.tsv` | Repeat coordinates, counts, patterns, motifs, and quality per read. |
| `vntr_asv_table.tsv` | Retained VSEARCH variants and representative statistics. |
| `vntr_asv_memberships.tsv` | Cluster membership and per-read substitution/indel evidence. |
| `vntr_asv_representatives.fasta` | Observed representative repeat sequences. |
| `vntr_mixture_abundance.tsv` | Emu-inspired EM abundance, estimated read count, and dominant/confirmed-secondary/candidate/trace evidence tier for every retained variant. |
| `locus_mapping_references.fasta` | Dominant observed representative amplicons. |
| `locus_read_alignments.sam` | minimap2 locus-relative mappings. |
| `locus_mapping_summary.tsv` | Mapping rate, depth, coverage, and SNP totals by locus. |
| `locus_snps.tsv` | Filtered representative-relative SNP evidence. |
| `read_level_allele_predictions.tsv` | Per-read repeat-count probabilities, unrounded measurement, uncertainty, and evidence weights. |
| `allele_calls.tsv` | Assembly-equivalent dominant-product calls plus capped read confidence, product size, raw repeat measurement, measurement source, total and primary depth, candidate/confirmed-secondary counts, and dominant mixture fraction. |
| `amplirust/` | Native primer-pairing evidence. |
| `vsearch/` | Native unique sequences, centroids, and memberships. |
| `minimap2/` | Dominant-locus mapping FASTQ and reference diagnostics. |

Recruitment presence statuses:

| Status | Meaning |
| --- | --- |
| `PRESENT_GENOTYPED` | At least two complete recruited products support genotyping. |
| `PRESENT_PROVISIONAL` | One complete product or repeat-boundary-spanning partial evidence is available. |
| `PRESENT_UNTYPED` | Locus-specific mapping establishes presence but does not resolve both repeat boundaries. |
| `NO_EVIDENCE` | No mapping passed recruitment thresholds. |

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
| `legacy_output.csv` | Historical row-oriented locus details, including zero-based primer positions and mismatch display. |
| `legacy_mlva_analysis.csv` | Historical wide repeat-count layout. |
| `legacy_predicted_pcr_sizes.csv` | Historical wide product-size layout. |
| `legacy_primer_mismatches.txt` | Historical primer mismatch summary. |

## FASTQ call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Sufficient depth and a decisive in-range posterior. |
| `LOW_DEPTH` | Fewer dominant-cluster reads than `--min-depth`. |
| `AMBIGUOUS` | Weak top posterior or insufficient separation from the second call. |
| `OUT_OF_RANGE` | Best repeat count lies outside the configured locus range. |
| `MULTIPLE_VARIANTS` | At least one confirmed secondary remains in metagenome mode; isolate mode additionally requires dominant fraction below 0.8. Candidate and trace variants do not force this status. |
| `LOCUS_DROPOUT` | No retained read evidence produced a prediction. |

## Assembly call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Product found and repeat count calculated. |
| `AMBIGUOUS` | Product found, but the leading allele is below the configured posterior threshold or leads the runner-up by less than 0.2. |
| `PRESENT_COUNT_UNKNOWN` | Product found but panel metadata cannot convert its size to repeat count. |
| `NOT_FOUND` | No accepted paired-primer product found. |
