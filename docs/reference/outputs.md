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
| `locus_repeat_counts.tsv` | Conventional individual-locus repeat counts in a compact long-form table; numeric counts can include estimates. |
| `allele_probability_distribution.tsv` | Ranked integer/half-unit allele probabilities, selected state, and inference method. |
| `common_locus_calls.tsv` | Technology-neutral FASTQ calls with probability, margin, informative molecule support, total presence support, and explicit called/low-coverage/unresolved/ambiguous/mixed/not-found status. |
| `molecule_candidate_evidence.tsv` | Optional detailed competitive-alignment and VNTR evidence for each molecule/candidate state; intended for validation rather than the default report. |
| `candidate_mapping/candidate_alignments.bam` | Optional (`--keep-intermediates`) compressed competitive candidate alignments streamed from minimap2 through htslib. It is removed during normal operation; no text candidate SAM is produced. |
| `short_read_diagnostics.tsv` | Paired-end calls only: unique recruited and rescued fragments, direct count, original candidate estimate, graph attempt/outcome/length, final count and call method. Dedicated audit output; primary call columns are unchanged. |
| `short_read_microassembly.fasta` | Paired-end calls with `--keep-intermediates`: successful read-supported local products, including an empty file when none were recovered. |
| `genotype_evidence.tsv` | FASTQ uncertainty status, independent support, anchors, strands, candidate support, observation envelope, error-model source, and computational fallback. |
| `sequence_alleles.tsv` | FASTQ complete products identically observed in at least three independent molecules, with SHA-256 identifiers, repeat sequence, pattern, interruptions, and partial terminal bases. |
| `library_insert_distribution.tsv` | Illumina effective insert mean/SD and calibration source (`mapped_unique_flanks`, `user_supplied`, or `unavailable`), support, and selection-bias limitation. |
| `overlap_fragment_calibration.tsv` | Illumina bounded exact-overlap calibration diagnostics, retained molecule counts, component mean/SD, and heterogeneity or censoring rejection. An accepted component has `whole_library=False` and `unobserved_tail_excluded=False`; it does not replace whole-library calibration. |
| `joint_pair_diagnostics.tsv` | Illumina paired-read rescue outcome per locus, including molecule support, competition, calibration source, and reasons for abstaining. |
| `locus_group_calls.tsv` | Illumina shared alleles and original molecule support for indistinguishable locus groups; separate from individual locus assignments and fingerprints. |
| `fastq_performance.tsv` | FASTQ mapping, recruitment/extraction, and inference/reconstruction timings, mapped-molecule throughput, and inference-tier fractions. |
| `variant_refinement_skips.tsv` | Long-read diagnostic variant alignments skipped because refinement would exceed resource limits. |

## Dataset-level MYOGA export outputs

`mlvamaps export-myoga` reads completed calls and writes the following without
rerunning sample analysis:

| File | Meaning |
| --- | --- |
| `myoga_metadata.tsv` | Original metadata restricted/padded to final tree samples, with canonical `sample_id`, validated `latitude`, and validated `longitude`. |
| `mlva_profiles.tsv` | Complete sample-by-locus numerical predictions; join `mlva_calls_long.tsv` for evidence and prior-only flags. Distance matrices separately exclude prior-only predictions. |
| `mlva_calls_long.tsv` | Auditable final sample-by-locus calls retaining available `calls.tsv` evidence fields. |
| `mlva_pairwise_distances.tsv` | Both categorical and repeat-count distance components over shared exact calls, including unsupported pre-tree pairs. |
| `mlva_distance_matrix.tsv` | Symmetric selected-distance matrix for all samples passing per-sample thresholds, with a zero diagonal and blank cells where no shared exact call permits a distance. |
| `mlva_nj.tree` | Deterministic neighbor-joining MLVA relatedness tree; absent when no sample passes filtering. |
| `samples_used.tsv` | Final tree sample paths, callability, metadata, and coordinate status. |
| `samples_excluded.tsv` | Tree exclusions and non-fatal geography exclusions with reason codes and details. |
| `export_summary.tsv`, `export_summary.txt` | Machine- and human-readable discovery, filtering, overlap, and output totals. |

With `--combined-markers`, `combined_marker_distance_matrix.tsv` and
`combined_marker_metadata.tsv` retain all discovered valid samples. Pairwise
cells are blank only when neither a shared exact repeat call nor a shared
recovered SNP sequence permits a distance; the
neighbor-joining tree remains limited to a finite complete submatrix.

`alignment_likelihood_distance_matrix.tsv` compares every discovered sample's
alignment-likelihood/EM composition using Hellinger distance. The matching
pairwise table, metadata, and neighbor-joining tree are also written. This is a
symmetric distance between fitted reference-component profiles, not a genomic
SNP distance; cells remain blank when a sample has no usable `classification.json`.

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

Database-backed FASTQ calls retain direct measurements and also report
low-confidence `ESTIMATED` counts when direct sizing is unresolved. The fallback
uses repeat-spanning CIGAR indels with sufficient aligned support at both flanks,
or calibrated opposite-flank pairs. Conflicting measurements, partial-repeat bounds,
and reference-conditioned alignments retain numerical estimates and distributions. Estimates
are included in calls and fingerprints, but not promoted to measured molecule
evidence for reference classification. Model support compares available alleles;
it is not a calibrated probability of correctness.

`estimation_method`, `estimation_reason`, and `estimate_support` are appended to
`common_locus_calls.tsv`, `calls.tsv`, `locus_repeat_counts.tsv`, and long-read
`allele_calls.tsv`. They identify the method, uncertainty, and supporting molecule
count. `molecule_candidate_evidence.tsv` also records `aligned_flank_bases` and
`spanning_repeat_indel_delta` for auditing alignment estimates. Estimates do not
inflate measured `molecule_support`.

`presence_support` counts detected independent molecules, while
`molecule_support` counts informative molecules. `primary_molecule_support` and
`secondary_molecule_support` count observations whose rounded measurement matches
each reported allele. `calls.tsv` uses total presence for `read_depth` and the
observed primary count for `primary_read_depth`; fitted allele fractions do not
create read counts. Mates and alternative mappings share their molecule identity.
Confidence thresholds never erase these observed depths or the selected repeat
prediction.

`common_locus_calls.tsv` also includes `repeat_count_raw`, `inference_method`,
`probability_kind`, `allele_fractions`, `em_iterations`, `em_converged`, and
`consensus_status`. `best_probability` and `candidate_distribution` describe
single-allele model support; `allele_fractions` and dominant/secondary fractions
describe the mixture model and are not probabilities of a correct call.
The model supports have not been calibrated against real-sample error rates.

`molecule_candidate_evidence.tsv` preserves `raw_repeat_count` and `repeat_sigma`
in repeat units, plus optional `pair_repeat_count` and `pair_repeat_sigma`.
Paired selection limits and `fragment_*` geometry fields retain the inputs to the
conditional length likelihood, including `fragment_length_class` values
`boundary_pair` and `repeat_mate`. Methods `paired_fragment_length`,
`short_read_joint_length_likelihood`, and `consensus_supported` identify estimates
from external or conservatively estimated library calibration and supported long-read consensus,
respectively. Neither adds independently measured classification evidence.

Method `joint_pair_read_likelihood` (`conditional_joint_pair_likelihood` in the
common inference table) identifies the additional joint paired-read placement
rescue. It can use an external/unique-flank calibration or a separately guarded
overlap component. Its point estimate remains conditional on the pure-repeat
and fragment-component assumptions: it is `ESTIMATED`, leaves the upper length
unresolved, and creates neither a reconstructed query sequence nor measured
classification evidence.

`report.html` includes an always-visible table for every panel locus, with the
count, status, support, and evidence or uncertainty reason. `≈` marks estimates.
FASTQ counts are mandatory for every valid configured locus with an eligible prior,
including no-data loci. Presence statuses such as `NOT_FOUND`, `LOCUS_DROPOUT`, and
`PRESENT_COUNT_UNKNOWN` describe observed evidence and can coexist with a numeric
prediction. Zero is emitted only when it is the selected allele. A locus without
informative evidence or an eligible prior fails with a configuration error.

### Complete FASTQ predictions and confidence

The shared finalizer runs after reconstruction and inference. Its priority is:
joint direct measurements, conflicting measured candidates, calibrated paired
lengths, repeat-spanning alignments, partial-repeat bounds, reference-conditioned
alignments, then a locus-specific prior. It preserves the selected number when
probability, margin, coverage, mixture, or computational thresholds fail.

Priors give equal mass to **distinct observed reference alleles**, irrespective
of duplicate references. If none exist, a panel with product-size calibration
and a repeat unit supplies its positive nominal count; otherwise an explicit
non-default repeat range supplies a uniform discrete prior. The default `0..100`
range is never informative. The legacy nominal default of zero is also not an
informative prior; an observed zero allele or explicit zero range remains valid.
Ties select the smaller count. Large priors use a
bounded grid and record the limit. Priors describe a catalog or panel assumption,
not population allele frequencies, presence, or evidence of a sample genotype.

The primary and genotype tables report the following fields. Confidence metadata
also travels with compact counts, allele calls, probabilistic fingerprints and
MYOGA locus exports. Existing aliases such as `called_repeat_count` and
`primary_repeat_count_raw` remain in `allele_calls.tsv`.

| Field | Meaning |
| --- | --- |
| `repeat_count`, `repeat_count_raw` | Selected reported and underlying estimates; raw values survive uncertainty thresholds. |
| `estimate_source` | Measurement, reconstruction, paired length, spanning alignment, partial-bound likelihood, reference-conditioned evidence, or the explicit prior source. |
| `genotype_status` | Evidence classification, independent of numerical completion. |
| `model_support` | Relative support among tested candidates; not calibrated accuracy. |
| `data_evidence_confidence` | Uncalibrated evidence-strength score: selected support × independent depth/(depth + minimum depth), discounted by source, measurement sigma, ambiguity and resource limits. Exactly zero for prior-only predictions. |
| `second_best_repeat_count`, `second_best_support` | Strongest alternative and its model support, when represented. |
| `repeat_interval_low`, `repeat_interval_high` | Observation envelope or likelihood support interval; not a calibrated credible interval. The upper field is blank when the upper length remains unresolved. |
| `candidate_grid_low`, `candidate_grid_high` | Computational search endpoints, when reported; not biological bounds. |
| `upper_unbounded`, `finite_length_support` | Whether available length evidence leaves the upper length unresolved or supports a finite interval in the fitted model. |
| `paired_span_support`, `partial_bound_support`, `repeat_rich_mate_support` | Unique molecules supplying paired geometry, anchored bounds, and repeat-only mate observations. Categories can overlap and should not be summed as independent depth. |
| `length_model_status` | Conditional length fit, bounds-only result, unsupported sample mode, or computational fallback status. |
| `length_null_fraction` | Selected zero or 1% background assumption in the conditional fragment model; uncalibrated. |
| `pair_rescue_support` | Distinct original molecules contributing to the joint paired-read placement rescue. |
| `pair_rescue_boundary_support` | Informative junction molecules supporting the rescue's boundary constraints; a subset of its molecule support. |
| `pair_rescue_calibration_source` | Fragment calibration used by the rescue, including the explicitly selected `exact_overlap_component`. |
| `pair_rescue_status` | Rescue outcome retained with the call; successful fits remain `estimated`, with no finite upper-length claim. |
| `independent_support` | Distinct contributing molecules; mates and duplicate references add no observations. Zero for prior-only predictions. |
| `prior_only` | True when no informative sample repeat evidence contributed. Presence can still be ambiguous or flank-supported. |
| `computational_fallback` | Candidate/matrix/refinement limit and bounded fallback reason. |
| `confidence_reason` | Explanation of evidence and limitations. |

Complete numerical fingerprints agree with the primary per-locus table. Profile
matching, repeat-profile trees, MYOGA distances and SNP-plus-repeat distances
exclude prior-only values. Reference classification continues to use original
alignments and measured molecule counts; finalized predictions create no new
classification evidence. Legacy assembly rows retain their existing behavior.

## FASTQ outputs

### Illumina-specific evidence

| File | Meaning |
| --- | --- |
| `short_read_qc_summary.tsv` | Input/retained reads and pairs, orphans, and legacy user-supplied insert-size settings. Use `library_insert_distribution.tsv` for the effective calibration. |
| `short_read_recruitment_summary.tsv` | Unique, ambiguous, discordant, and orphan pair counts per locus. |
| `short_read_mapping_evidence.tsv` | Per-locus state, candidate scores, molecule support, boundary evidence, MAPQ, and context provenance. |
| `short_read_run_metadata.json` | minimap2 and mlvamaps versions, configured parameters, database source, and legacy external insert-size settings. Automatic calibration is recorded in `library_insert_distribution.tsv`. |
| `filtered_reads_1.fastq.gz`, `filtered_reads_2.fastq.gz` | Quality-filtered mates, written with fast gzip compression for downstream native recruitment. |
| `filtered_orphan_reads.fastq.gz` | Retained single mates whose partner failed QC; empty when no orphans are present. |
| `sample_summary.tsv` | One normalized sample row for batch aggregation. |
| `myoga_samples.csv` | MYOGA metadata; `genome_id` equals `sample_id` and generated sample tree-tip IDs. |
| `myoga_loci.csv` | Long-form numerical predictions, intervals, evidence, confidence, and prior-only flags; also produced by long-read FASTQ calls. |
| `batch_summary/batch_status.tsv` | Success, failure, or resume status for every directory or manifest sample. |

Illumina `calls.tsv` preserves the compact columns and appends read technology,
evidence class, boundary support, informative molecules, local assembly,
repeat interval, confidence explanation, failure/warning, and mixture-support
fields. Current valid configured FASTQ loci retain a numeric `repeat_count`;
`repeat_interval_low`/`high` describe measured or candidate envelopes separately.
Older results can have empty counts and populated `repeat_count_min`/`max`, or
empty sizing fields with `PRESENCE_ONLY`. Those historical formats remain readable.

`locus_group_calls.tsv` reports `member_loci`, a stable `group_id`, shared repeat
alleles, independent molecule support, and `group_status`. Its `locus_assignment`
remains `AMBIGUOUS_LOCUS`. These rows do not appear as additional loci in the
fingerprint or as classification observations. `SHARED_ALLELE` establishes a
pooled allele measurement; `MULTIPLE_ALLELES` retains distinct group alleles.
Group fractions cannot identify which member carries an allele or establish
strain fractions. `LOWER_BOUND` and `UNRESOLVED` rows retain absent full sizing
without filling a group count from the reference prior.

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

The compatible tables below retain legacy statuses. **`PASS` does not mean
`EXACT`.** Consult `genotype_status` in `common_locus_calls.tsv` or the
supplementary `genotype_evidence.tsv`:

| Evidence status | Meaning |
| --- | --- |
| `EXACT` | A resolved call with at least `max(3, minimum_molecules)` independent complete products, agreeing conventional counts, measurement sigma at most 0.25 repeats, no quality flags or ambiguous bases, and no derived consensus votes. These are evidence criteria, not a calibrated error guarantee. |
| `ESTIMATED` | A numeric count failing the exact-evidence criteria, including low coverage, uncertain measurements, reconstruction, and paired-length estimates. |
| `LOWER_BOUND` | A uniquely anchored terminal motif run supplies `repeat_lower_bound`; the numeric prediction depends on a conditioned prior and an explicit expansion-tail assumption. |
| `UNRESOLVED` | Evidence or computation cannot resolve a genotype confidently; a bounded numeric prediction is still retained. |
| `PRIOR_ONLY` | No informative sample repeat evidence; numerical prediction comes entirely from the explicit locus prior and data-evidence confidence is zero. |
| `AMBIGUOUS_LOCUS` | Only competing locus assignments support the locus; those molecules do not vote for a repeat count. |
| `POSSIBLE_MIXTURE` | Supported secondary-allele evidence; this does not confirm a strain mixture. |

`repeat_interval_low`/`high` enclose the measured observations plus/minus three
measurement sigmas; they are **not calibrated confidence intervals**.
`log_likelihood_gap` compares the leading raw allele states. Legacy fields named
`probability` describe model support, not validated probabilities of correctness.
Inference tiers describe direct collapsed evidence (1), targeted motif refinement
(2), and joint heterogeneous/paired inference (3); zero means no usable sizing
observations. The highest applicable tier is reported.

Sequence allele IDs hash the observed complete product, including flanks.
Distinct products can have the same repeat count. These supplementary products
are not automatically added to combined-marker distances; the existing exporter
continues to mask repeat regions in its SNP component. TSV replacement is atomic
per file; the output directory is not a transactional snapshot of an entire run.
See the [validation and compatibility report](../workflows/fastq-evidence-validation.md).

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
