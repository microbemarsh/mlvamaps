# Output file reference

All generated FASTA and FASTQ artifacts are gzip-compressed by default and use
a matching `.gz` suffix. Input files are never modified.

## Files returned for every call

| File | Meaning |
| --- | --- |
| `calls.tsv` | Compact per-locus result shared by FASTQ and assembly modes. |
| `mlva_fingerprint.tsv` | Wide sample-by-locus repeat-count fingerprint. |
| `mlva_fingerprint_probabilistic.tsv` | Long-form calls with confidence values. |
| `profile_matches.tsv` | Complete machine-readable match table corresponding to the HTML results. Rows are labeled `mlva_profile` for direct repeat-profile comparisons and `sequence_reference` for combined marker references. Includes every ranked match, metadata, distances, comparison counts, and source-specific SNP/repeat fields. Header-only only when neither `--profiles` nor `--database` supplies reference matches. |
| `profile_match_loci.tsv` | Long-form profile comparison with one row per profile and locus, including query/profile alleles, absolute difference, match status, and profile-allele probability. |
| `report.html` | Self-contained interpretation report with sample findings, locus-quality flags, FASTQ SPOARS/assembly-PCR concordance, gel evidence, profile matches, closest reference genomes, and technical tables. |
| `locus_repeat_counts.tsv` | Exact individual-locus repeat counts in a compact long-form table. |
| `allele_probability_distribution.tsv` | Ranked integer/half-unit allele probabilities, selected state, and inference method. |

Multi-taxon databases additionally produce the following under `phylogeny/`:

| File | Meaning |
| --- | --- |
| `taxonomic_identification.tsv` | Backward-compatible summary columns plus assignment rank/status, categorical confidence, closest and runner-up distances, margins, locus support, bootstrap stability, and input mode. |
| `taxonomic_identification_evidence.tsv` | Ranked candidate taxa with distance, similarity (not probability), and bootstrap winner fraction. |
| `taxonomic_identification_loci.tsv` | Per-locus recovery, taxonomic weight, favored taxon, support/conflict state, and available FASTQ depth/consensus evidence. |
| `taxonomic_identification.json` | Complete versioned machine-readable assignment, candidate, and locus evidence object. |

Illumina calls additionally always write `sample_summary.tsv`,
`myoga_samples.csv`, and `myoga_loci.csv`.

## Dataset-level MYOGA export outputs

`mlvamaps export-myoga` reads completed calls and writes the following without
rerunning sample analysis:

| File | Meaning |
| --- | --- |
| `myoga_metadata.tsv` | Original metadata restricted/padded to final tree samples, with canonical `sample_id`, validated `latitude`, and validated `longitude`. |
| `mlva_profiles.tsv` | Final sample-by-locus exact repeat-count matrix; unresolved loci are empty. |
| `mlva_calls_long.tsv` | Auditable final sample-by-locus calls retaining available `calls.tsv` evidence fields. |
| `mlva_pairwise_distances.tsv` | Both categorical and repeat-count distance components over shared exact calls, including unsupported pre-tree pairs. |
| `mlva_distance_matrix.tsv` | Symmetric selected-distance matrix for final tree samples, with a zero diagonal. |
| `mlva_nj.tree` | Deterministic neighbor-joining MLVA relatedness tree; absent when no sample passes filtering. |
| `samples_used.tsv` | Final tree sample paths, callability, metadata, and coordinate status. |
| `samples_excluded.tsv` | Tree exclusions and non-fatal geography exclusions with reason codes and details. |
| `export_summary.tsv`, `export_summary.txt` | Machine- and human-readable discovery, filtering, overlap, and output totals. |

With `export-myoga --combined-markers`, the export additionally writes:

| File | Meaning |
| --- | --- |
| `combined_marker_sequence_status.tsv` | Per-sample/locus recovery source, masking method, and explicit reason when no safe sequence was usable. |
| `locus_tree_status.tsv` | Per-locus sample and SNP-haplotype counts, inference method, scale, alignment, and tree path. |
| `locus_trees/LOCUS/samples.tree` | Sample-tip locus tree derived from the repeat-masked SNP tree distance; identical SNP haplotypes are restored as zero-distance tips. |
| `locus_snp_distances.tsv` | Per-locus pairwise patristic SNP and repeat distances, locus scales, normalized components, and weighted locus distance. |
| `combined_marker_pairwise_distances.tsv` | Component means, shared-locus counts, weights, combined distances, and overlap status for each sample pair. |
| `combined_marker_distance_matrix.tsv` | Complete combined-marker matrix after deterministic overlap pruning. |
| `combined_marker_nj.tree` | Neighbor-joining tree from the combined marker matrix. |
| `combined_marker_metadata.tsv` | Metadata restricted to exactly the combined-tree tips. |

See the [MYOGA export workflow](../workflows/myoga-export.md) for formulas and
filter behavior.

## Reference builder outputs

| File | Meaning |
| --- | --- |
| `database/LOCUS.fasta.gz` | Gzip-compressed unmasked reference amplicons accepted by `--database`. |
| `database/reference_metadata.tsv` | Metadata normalized to a `reference_id` key. |
| `database/reference_sequence_index.tsv` | Canonical amplicon, repeat-masked SNP, and complete marker SHA-256 keys used by the default exact-match fast path. |
| `database/taxon_locus_discrimination.tsv` | Build-time per-locus taxonomic weight, normalized information gain, reference coverage, and supporting counts. |
| `database/reference_assemblies.tsv` | Reference-ID, source-assembly path, and canonical whole-genome SHA-256 used for assembly-query tie breaking. |
| `reference_build_manifest.tsv` | Per-reference/locus product counts, selected product, primer errors, and exclusion status. |
| `reference_locus_amplifiability.tsv` | Per-locus retained amplicon and genome counts, amplifiable percentage, and `NO_AMPLICONS`, `INSUFFICIENT_REFERENCES`, or `BUILT` tree status. |
| `phylogeny/LOCUS.tree` | Portable Newick SNP tree for the locus. |
| `phylogeny/LOCUS/references.aligned.fasta.gz` | Repeat-masked MAFFT alignment used for tree inference. |
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

Before launching alignment or placement, mlvamaps checks the persistent sequence
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
| `phylogeny/taxonomic_identification.tsv` | Automatic sample-level nearest-taxon result (`SUPPORTED`, `AMBIGUOUS`, or `INSUFFICIENT_EVIDENCE`) when reference metadata contain `taxon_id`. |
| `phylogeny/taxonomic_identification_evidence.tsv` | Ranked per-taxon distance, coverage-adjusted score, reference count, and informative-locus evidence. |
| `phylogeny/taxon_assignment.tsv` | Optional calibrated `POSITIVE`, `NEGATIVE`, or `INDETERMINATE` target-taxon result, prediction set, compatibility p-values, bootstrap support, and QC. Written when `--target-taxon-id` and `--taxon-calibration` are supplied. |
| `phylogeny/taxon_assignment_candidates.tsv` | Repeat, SNP, and joint distances, conformal nonconformity/p-values, acceptance state, and nearest references for every labeled taxon. |
| `phylogeny/taxon_assignment_loci.tsv` | Per-locus target-versus-best-alternative marker distances, margin, placement uncertainty, and favored state. |

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

### Illumina-specific evidence

| File | Meaning |
| --- | --- |
| `short_read_qc_summary.tsv` | Input/retained reads and pairs, orphans, and empirical insert-size values when estimable. |
| `short_read_recruitment_summary.tsv` | Unique, ambiguous, discordant, and orphan pair counts per locus. |
| `short_read_assembly_summary.tsv` | Per-locus SKESA status, contig sizes, depth, and failure reason. |
| `filtered_reads_1.fastq.gz`, `filtered_reads_2.fastq.gz` | Quality-filtered mates, written with fast gzip compression for downstream native recruitment. |
| `filtered_orphan_reads.fastq.gz` | Retained single mates whose partner failed QC; empty when no orphans are present. |
| `sample_summary.tsv` | One normalized sample row for batch aggregation. |
| `myoga_samples.csv` | MYOGA metadata; `genome_id` equals `sample_id` and generated sample tree-tip IDs. |
| `myoga_loci.csv` | Long-form exact calls, intervals, evidence, and confidence. |
| `batch_status.tsv` | Success, failure, or resume status for every manifest sample. |

Illumina `calls.tsv` preserves the compact columns and appends read technology,
evidence class, boundary support, informative molecules, local assembly,
repeat interval, confidence explanation, failure/warning, and mixture-support
fields. Empty `repeat_count` plus populated `repeat_count_min` and
`repeat_count_max` is an interval. Empty exact and interval fields with
`PRESENCE_ONLY` means detected but not sized.

When Illumina mode receives `--database`, complete primer-bounded locus products
also produce the standard `phylogeny/` reference-placement outputs. Their
ranked `sequence_reference` rows are appended to `profile_matches.tsv` and
rendered in `report.html`. Unresolved loci are omitted from sequence placement.

| File | Meaning |
| --- | --- |
| `qc_summary.tsv` | Input and read-filtering totals. |
| `taxon_screen/taxon_screened_reads.fastq.gz` | Optional reads retained by the target-taxon Deacon pangenome screen before downstream analysis. |
| `taxon_screen/taxon_screen_summary.json` | Native Deacon input/output read and base totals, thresholds, and throughput. |
| `filtered_reads.fastq.gz` | Reads retained after length and quality filtering. |
| `filtered_reads.fasta.gz` | Lossless sequence projection used by native primer pairing. |
| `locus_recruited_reads.tsv` | Competitive per-read locus mapping, candidate allele, alignment quality, and presence/genotype evidence class. |
| `locus_presence.tsv` | Per-locus mapped, full-product, and repeat-informative read counts with presence status. |
| `local_locus_products.fasta.gz` | SPOARS POA consensus contigs reconstructed from the dominant cluster and passed through assembly-mode PCR calling. |
| `local_assembly_concordance.tsv` | Per-locus raw read-length range and mode, POA consensus length, PCR product size, raw/final repeat calls, support, and fallback status. |
| `local_assembly_pcr/` | Native Sassy primer matches and extracted products from the per-locus POA contigs. |
| `recruitment/locus_recruitment_references.fasta.gz` | Database-derived or synthetic competitive locus/allele reference bank. |
| `recruitment/read_recruitment.sam` | Raw competitive long-read mappings used for presence and provisional genotype evidence. |
| `read_locus_assignments.tsv` | Primer-supported locus, orientation, score, and assignment QC. |
| `read_repeat_features.tsv` | Repeat coordinates, counts, patterns, motifs, and quality per read. |
| `mapped_variant_table.tsv` | Competitive mapping-derived repeat-product groups and support statistics. |
| `mapped_read_memberships.tsv` | Per-read mapping group plus substitution/indel evidence. |
| `mapped_variant_representatives.fasta.gz` | Diagnostic repeat representatives for mapped product groups. |
| `vntr_mixture_abundance.tsv` | EM abundance, estimated read count, and dominant/confirmed-secondary/candidate/trace evidence tier for every mapped product group. |
| `locus_mapping_references.fasta.gz` | Dominant SPOARS assembly-PCR products used for read support and SNP mapping. |
| `locus_read_alignments.sam` | minimap2 locus-relative mappings. |
| `locus_mapping_summary.tsv` | Mapping rate, depth, coverage, and SNP totals by locus. |
| `locus_snps.tsv` | Filtered representative-relative SNP evidence. |
| `read_level_allele_predictions.tsv` | Per-read repeat-count probabilities, unrounded measurement, uncertainty, and evidence weights. |
| `read_locus_disagreement_audit.tsv` | Optional (`--debug-disagreements`) CIGAR, extraction, anchor, repeat, mapping-reference allele, and measured-read allele evidence per recruited read. |
| `locus_disagreement_summary.tsv` | Optional (`--debug-disagreements`) mapping/measurement disagreement counts and combined-read versus consensus call agreement. |
| `allele_calls.tsv` | Assembly-equivalent dominant-product calls plus capped read confidence, product size, raw repeat measurement, measurement source, total and primary depth, candidate/confirmed-secondary counts, and dominant mixture fraction. |
| `in_silico_pcr/` | Native primer-pairing evidence from filtered reads. |
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
| `assembly_amplicons.fasta.gz` | Extracted assembly primer products. |
| `read_support.tsv` | Optional mapped reads and mean coverage per product. |
| `read_support.sam` | minimap2 alignments when `--reads` is used. |
| `amplirust/` | Native in-silico PCR evidence. |
| `legacy_output.csv` | Historical row-oriented locus details, including zero-based primer positions and mismatch display. |
| `legacy_mlva_analysis.csv` | Historical wide repeat-count layout. |
| `legacy_predicted_pcr_sizes.csv` | Historical wide product-size layout. |
| `legacy_primer_mismatches.txt` | Historical primer mismatch summary. |

For a directory input containing assemblies, mlvamaps also writes
`MLVA_analysis_<input-directory>.csv` at the top level of the output directory.
It combines all per-assembly `legacy_mlva_analysis.csv` rows in deterministic
filename order and assigns the historical zero-padded `key` values.

## FASTQ call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Sufficient depth and a decisive in-range posterior. |
| `LOW_DEPTH` | Fewer dominant-cluster reads than `--min-depth`. |
| `AMBIGUOUS` | Weak top posterior or insufficient separation from the second call. |
| `OUT_OF_RANGE` | Best repeat count exceeds the configured review range by more than `--repeat-range-tolerance`. The observed allele is retained rather than clipped. |
| `MULTIPLE_VARIANTS` | At least one confirmed secondary remains in metagenome mode; isolate mode additionally requires dominant fraction below 0.8. Candidate and trace variants do not force this status. |
| `LOCUS_DROPOUT` | No retained read evidence produced a prediction. |

`allele_calls.tsv` also contains a more explicit `evidence_status`: `CONFIDENT`,
`PROVISIONAL_LOW_DEPTH`, `SINGLE_MOLECULE_PROVISIONAL`, `AMBIGUOUS`, or
`NO_INFORMATIVE_READS`. The legacy `call_status` values remain unchanged for
backwards compatibility.

## Assembly call statuses

| Status | Meaning |
| --- | --- |
| `PASS` | Product found and repeat count calculated. |
| `AMBIGUOUS` | Product found, but the leading allele is below the configured posterior threshold or leads the runner-up by less than 0.2. |
| `PRESENT_COUNT_UNKNOWN` | Product found but panel metadata cannot convert its size to repeat count. |
| `NOT_FOUND` | No accepted paired-primer product found. |
