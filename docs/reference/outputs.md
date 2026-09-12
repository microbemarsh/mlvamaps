# Output file reference

All generated FASTA and FASTQ artifacts are gzip-compressed by default and use
a matching `.gz` suffix. Input files are never modified.

## Single-sample and batch layouts

A single input writes files directly to the requested output directory:

```text
OUTDIR/
├── report.html
├── calls.tsv
├── mlva_fingerprint.tsv
└── ...
```

A directory or manifest run separates each sample from batch-wide aggregates:

```text
OUTDIR/
├── SAMPLE_A/
│   ├── report.html
│   ├── calls.tsv
│   └── ...
├── SAMPLE_B/
│   └── ...
└── batch_summary/
    ├── batch_status.tsv
    ├── calls.tsv
    ├── mlva_fingerprint.tsv
    └── ...
```

Files below a sample directory contain only that sample. Files below
`batch_summary/` combine available results across the run; failed samples are
recorded in `batch_summary/batch_status.tsv`. The `batch_summary` name is
reserved and cannot be used as a sample ID in a batch.

## Files returned for every call

| File | Meaning |
| --- | --- |
| `calls.tsv` | Compact per-locus result shared by FASTQ and assembly modes. |
| `mlva_fingerprint.tsv` | Wide sample-by-locus repeat-count fingerprint. |
| `mlva_fingerprint_probabilistic.tsv` | Long-form calls with confidence values. |
| `profile_matches.tsv` | Complete machine-readable match table corresponding to the HTML results. Rows are labeled `mlva_profile` for direct repeat-profile comparisons and `mapping_reference` for alignment-based references. Includes every ranked match, metadata, distances, comparison counts, and source-specific support fields. Header-only only when neither `--profiles` nor `--database` supplies reference matches. |
| `profile_match_loci.tsv` | Long-form profile comparison with one row per profile and locus, including query/profile alleles, absolute difference, match status, and profile-allele probability. |
| `report.html` | Self-contained interpretation report with sample findings, locus-quality flags, FASTQ SPOARS/assembly-PCR concordance, gel evidence, profile matches, closest reference genomes, and technical tables. |
| `locus_repeat_counts.tsv` | Exact individual-locus repeat counts in a compact long-form table. |
| `allele_probability_distribution.tsv` | Ranked integer/half-unit allele probabilities, selected state, and inference method. |
| `common_locus_calls.tsv` | Technology-neutral FASTQ calls with probability, margin, molecule support, and explicit called/low-coverage/unresolved/ambiguous/mixed/not-found status. |
| `molecule_candidate_evidence.tsv` | Optional detailed competitive-alignment and VNTR evidence for each molecule/candidate state; intended for validation rather than the default report. |
| `candidate_mapping/candidate_alignments.bam` | Optional (`--keep-intermediates`) compressed competitive candidate alignments streamed from minimap2 through htslib. It is removed during normal operation; no text candidate SAM is produced. |

## Dataset-level MYOGA export outputs

`mlvamaps export-myoga` reads completed calls and writes the following without
rerunning sample analysis:

| File | Meaning |
| --- | --- |
| `myoga_metadata.tsv` | Original metadata restricted/padded to final tree samples, with canonical `sample_id`, validated `latitude`, and validated `longitude`. |
| `mlva_profiles.tsv` | Final sample-by-locus exact repeat-count matrix; unresolved loci are empty. |
| `mlva_calls_long.tsv` | Auditable final sample-by-locus calls retaining available `calls.tsv` evidence fields. |
| `mlva_pairwise_distances.tsv` | Both categorical and repeat-count distance components over shared exact calls, including unsupported pre-tree pairs. |
| `mlva_distance_matrix.tsv` | Symmetric selected-distance matrix for all samples passing per-sample thresholds, with a zero diagonal and blank cells where no shared exact call permits a distance. |
| `mlva_nj.tree` | Deterministic neighbor-joining MLVA relatedness tree; absent when no sample passes filtering. |
| `samples_used.tsv` | Final tree sample paths, callability, metadata, and coordinate status. |
| `samples_excluded.tsv` | Tree exclusions and non-fatal geography exclusions with reason codes and details. |
| `export_summary.tsv`, `export_summary.txt` | Machine- and human-readable discovery, filtering, overlap, and output totals. |

With `--combined-markers`, `combined_marker_distance_matrix.tsv` and
`combined_marker_metadata.tsv` retain all discovered valid samples. Pairwise
cells are blank only when no shared recovered marker permits a distance; the
neighbor-joining tree remains limited to a finite complete submatrix.

## Reference builder outputs

| File | Meaning |
| --- | --- |
| `database/LOCUS.fasta.gz` | Gzip-compressed unmasked reference amplicons accepted by `--database`. |
| `database/reference_metadata.tsv` | Metadata normalized to a `reference_id` key. |
| `database/reference_assemblies.tsv` | Reference-ID, source-assembly path, and canonical whole-genome SHA-256 retained for source provenance. |
| `database/competitive_mapping/candidate_contexts.fasta` | Unique observed or synthetic repeat-state sequence hypotheses used for competitive FASTQ mapping. |
| `database/competitive_mapping/candidate_metadata.tsv` | Stable candidate IDs, repeat coordinates/states, background IDs, taxon data, and observed/synthetic status. |
| `database/competitive_mapping/candidate_provenance.tsv` | Candidate-to-reference, taxon, and background relationships retained after sequence deduplication. |
| `database/competitive_mapping/short.mmi` | minimap2 index using `-k 21 -w 11` for Illumina reads. |
| `database/competitive_mapping/long.mmi` | minimap2 index using `-k 11 -w 5` for accurate long/amplicon reads. |
| `database/deacon/reference_genomes.fasta` | Combined real genomes used for broad target-group recruitment; synthetic alleles are excluded. |
| `database/deacon/target_recruitment.idx` | Deacon minimizer index over the real reference genomes. |
| `manifest.json` | Completed schema-2.0 database manifest with checksums, versions, parameters, and asset paths. |
| `reference_build_manifest.tsv` | Per-reference/locus product counts, selected product, primer errors, and exclusion status. |
| `reference_locus_amplifiability.tsv` | Per-locus retained amplicon and genome counts, amplifiable percentage, and `NO_AMPLICONS`, `INSUFFICIENT_REFERENCES`, or `BUILT` tree status. |
| `myoga_metadata.csv` | Metadata with a `genome_id` column matching Newick tip labels. |

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
| `short_read_mapping_evidence.tsv` | Per-locus state, candidate scores, molecule support, boundary evidence, MAPQ, and context provenance. |
| `short_read_run_metadata.json` | minimap2 and mlvamaps versions, resolved parameters, database source, and insert-size estimate. |
| `filtered_reads_1.fastq.gz`, `filtered_reads_2.fastq.gz` | Quality-filtered mates, written with fast gzip compression for downstream native recruitment. |
| `filtered_orphan_reads.fastq.gz` | Retained single mates whose partner failed QC; empty when no orphans are present. |
| `sample_summary.tsv` | One normalized sample row for batch aggregation. |
| `myoga_samples.csv` | MYOGA metadata; `genome_id` equals `sample_id` and generated sample tree-tip IDs. |
| `myoga_loci.csv` | Long-form exact calls, intervals, evidence, and confidence. |
| `batch_summary/batch_status.tsv` | Success, failure, or resume status for every directory or manifest sample. |

Illumina `calls.tsv` preserves the compact columns and appends read technology,
evidence class, boundary support, informative molecules, local assembly,
repeat interval, confidence explanation, failure/warning, and mixture-support
fields. Empty `repeat_count` plus populated `repeat_count_min` and
`repeat_count_max` is an interval. Empty exact and interval fields with
`PRESENCE_ONLY` means detected but not sized.

When Illumina mode receives `--database`, complete primer-bounded locus products
also produce the standard `classification/` reference-support outputs. Their
ranked `sequence_reference` rows are appended to `profile_matches.tsv` and
rendered in `report.html`. Informative alignments can contribute to reference
support even when the repeat count is unresolved.

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
| `in_silico_pcr/` | Sassy-backed primer-search and paired-product evidence (`primers.csv`, `matches.tsv`, and `products.fasta.gz`). |
| `legacy_output.csv` | Historical row-oriented locus details, including zero-based primer positions and mismatch display. |
| `legacy_mlva_analysis.csv` | Historical wide repeat-count layout. |
| `legacy_predicted_pcr_sizes.csv` | Historical wide product-size layout. |
| `legacy_primer_mismatches.txt` | Historical primer mismatch summary. |

For a directory input containing assemblies, mlvamaps also writes
`batch_summary/MLVA_analysis_<input-directory>.csv`.
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
| `NOT_FOUND` | No product was eligible for the historical assembly repeat-count rule. Candidate products can still appear in `assembly_amplicons.tsv` when repeat calibration is unavailable or the raw repeat count is 100 or greater. |

### Coverage-gated missing-locus evidence

When read coverage passes the configured gate, confidently undetected loci
contribute a reference-specific log-likelihood penalty. This heuristic does not
establish biological absence. Assembly inputs do not use this gate. The
classification table records the settings, affected loci and penalty; these
penalties never enter repeat-profile trees. See [the model](../concepts/mapping-classification.md).

## Mapping classification (default with a sequence database)

| File under `classification/` | Meaning |
| --- | --- |
| `mapping_reference_matches.tsv` | Ranked reference groups with joint log likelihood, model weight, EM fraction, equivalent IDs and missing-locus penalty. |
| `molecule_reference_likelihoods.tsv.gz` | Original per-molecule/reference log scores before the sample-level missing-locus penalty. |
| `observed_reference_vntrs.fasta` and `observed_reference_metadata.tsv` | Actual observed mapping targets and their reference/locus membership; no synthetic alleles. |
| `assembly_reference_evidence.tsv` and `query_amplicons.fasta` | Assembly-only reference-relative mismatch counts, CIGARs, repeat differences and observed query products. |
| `closest_reference_bands.tsv` | Amplicon sizes and repeat counts for the leading reference, used in the report gel. |
| `classification.json` | Error-rate estimates, coverage gate, penalties, grouped results and EM objective/convergence diagnostics. |
| `taxonomic_identification.tsv` and `taxonomic_identification_evidence.tsv` | Closest-reference identification and ranked reference-group evidence, with equivalent reference IDs, optional taxon annotations and model support. Species support is not pooled. |
| `mlva_profiles.tree` | Newick neighbor-joining tree of observed repeat profiles; no classification likelihoods or absence penalties enter it. |
| `mlva_profile_distances.tsv` and `mlva_profile_tree_metadata.tsv` | Distance matrix and metadata whose IDs match the tree tips. |
| `mlva_profile_tree_status.tsv` | Locus set, distance definition, and explicit reason when no profile tree is available. |

Mapping results appear as `match_type=mapping_reference` in `profile_matches.tsv`.
The distance is negative joint log likelihood, not the legacy normalized
SNP/repeat distance. In mixed mode rank follows the fitted component fraction,
so distance alone need not be monotonic in rank. `equivalent_references` lists
IDs the evidence cannot distinguish; `reference_id` is only a representative.
Use the `mapping_reference_matches` result key; legacy combined-marker
keys and phylogenetic typing outputs are removed.

See [mapping classification and profile trees](../concepts/mapping-classification.md)
for the likelihood model, defaults, interpretation, and multi-sample MYOGA export.
