# Automatic MLVA taxon assignment

When `mlvamaps call` receives a reference database whose
`reference_metadata.tsv` contains `taxon_id`, it automatically ranks all
annotated taxa. No calibration artifact or requested target is needed. Use
`--no-taxon-identification` to disable this behavior, or
`--taxon-identification` to require metadata and fail clearly when it is absent.

For loci recovered in the query and represented in every candidate taxon,
mlvamaps reuses `locus_marker_distances.tsv`. At reference-build time it writes
`database/taxon_locus_discrimination.tsv`. Each locus weight is normalized
mutual information between its repeat/SNP marker signature and taxon label,
multiplied by reference coverage and the fraction of signatures observed more
than once. Thus shared markers receive zero or little weight and singleton
isolate signatures are discounted.

Runtime taxon distance is the weighted mean over recovered discriminatory loci,
followed by the mean of the nearest `--taxon-k` complete references (default 3).
The transparent display similarity is:

```text
taxon_score = locus_recovery_fraction / (1 + taxon_distance)
```

The similarity is not a probability. Species assignment additionally requires
compatible distance, enough discriminatory loci, runner-up separation,
non-conflicting evidence, and configured bootstrap stability. FASTQ/Illumina
calls require 90% locus recovery, three discriminatory loci, and 1.5 times the
assembly margin. Missing loci are never interpreted as taxon absence.

When species criteria fail, metadata columns `species_group`/`taxon_group`/`group`,
`genus`, `family`, `order`, `class`, and `phylum` are considered in that order.
The most specific value shared by the closest competitors is reported; otherwise
the assignment is unresolved.

Outputs are summary, candidate, and locus-evidence TSVs plus
`phylogeny/taxonomic_identification.json`.

## Optional calibrated target validation

`mlvamaps` can test whether a sample's MLVA repeat counts and repeat-masked
marker sequences are compatible with a requested target taxon. The method uses
the existing MAFFT, EPA-ng, direct aligned-SNP, and tandem-repeat evidence. It
does not require Kraken2, Deacon, SPUMONI, or whole-genome taxonomic
classification.

## Interpretation

The result is one of:

- `POSITIVE`: the target is the only taxon in the calibrated joint prediction
  set, repeat and SNP channels both support it, locus-bootstrap support passes,
  and assignment QC passes.
- `NEGATIVE`: the calibrated test excludes the target and supports at least one
  labeled alternative with adequate data.
- `INDETERMINATE`: multiple taxa remain compatible, no taxon is compatible,
  repeat and SNP evidence disagree, bootstrap support is insufficient, or
  assignment QC fails.

Conformal p-values are **reference-cohort compatibility values**. They are not
posterior probabilities that the organism is present. EPA-ng likelihood weight
ratios are retained as placement-QC evidence and are not interpreted as species
probabilities.

## Reference metadata

The sequence database must contain both target and relevant near-neighbor
references. Its `reference_metadata.tsv` requires these columns for assignment:

```text
reference_id	taxon_id	taxon_name
BA_001	1392	Bacillus anthracis
BC_001	1396	Bacillus cereus
```

`taxon_id` is the stable class label. `taxon_name` is optional display text.
All records assigned the same `taxon_id` must use a consistent name.

A target-only database cannot establish specificity. For closely related taxa,
include diverse lineages and known confounders. Identical MLVA marker
haplotypes occurring in multiple taxa should produce an indeterminate result.

## Calibration

Calibration is label-conditional and leave-one-reference-out. The input distance
table requires one row for every query reference, candidate reference, and
locus:

```text
query_reference_id	reference_id	locus_id	normalized_repeat_distance	normalized_snp_distance
BA_001	BA_001	L1	0	0
BA_001	BA_002	L1	0.2	0.1
BA_001	BC_001	L1	1.5	1.2
```

The calibrator always removes rows where `query_reference_id == reference_id`,
preventing a reference from calibrating against itself. Distances should be
generated with the same repeat masking, normalization, marker weights, and
reference trees used for runtime placement.

Build the artifact with:

```bash
mlvamaps calibrate-taxa \
  --reference-distances reference_leave_one_out_distances.tsv \
  --reference-metadata reference_build/database/reference_metadata.tsv \
  --sequence-index reference_build/database/reference_sequence_index.tsv \
  --k 3 \
  --alpha 0.05 \
  --minimum-loci 3 \
  --output reference_build/database/taxon_calibration
```

Outputs:

- `taxon_calibration.json`: versioned class counts, marker weights, thresholds,
  panel/database signatures, and per-class conformal score distributions.
- `taxon_calibration_scores.tsv`: auditable leave-one-reference-out scores.

The runtime rejects calibration artifacts whose panel or sequence-database
signature differs from the active reference build.

## Calling

Run the normal MLVA workflow with a database, target, and calibration artifact:

```bash
mlvamaps call \
  --panel panel.tsv \
  --input sample.fastq.gz \
  --database reference_build \
  --target-taxon-id 1392 \
  --taxon-calibration reference_build/database/taxon_calibration/taxon_calibration.json \
  --output results
```

The same options work for long-read FASTQ, Illumina, and assembly workflows.
The older target-specific conformal assignment remains available as an optional
advanced validation workflow and is disabled unless both `--target-taxon-id`
and `--taxon-calibration` are supplied. `calibrate-taxa` is retained for this
purpose; it is not required for normal automatic identification.

Important optional controls:

- `--taxon-alpha`: override conformal prediction-set alpha.
- `--taxon-min-loci`: override the artifact's minimum callable loci.
- `--taxon-min-locus-fraction`: required callable panel fraction.
- `--taxon-bootstrap-replicates`: deterministic locus bootstrap count (default 200).
- `--taxon-min-bootstrap-support`: winner support required for an automatic
  species assignment and target support required for `POSITIVE` (default 0.90).
- `--taxon-max-placement-entropy`: optional EPA-ng uncertainty ceiling.
- `--taxon-min-placement-lwr`: optional median best-placement LWR floor.

## Runtime outputs

Under `phylogeny/`:

- `taxon_assignment.tsv`: final decision, reason, prediction set, target and
  alternative compatibility values, bootstrap support, placement QC, and locus
  counts.
- `taxon_assignment_candidates.tsv`: one row per candidate taxon with repeat,
  SNP, and joint distances, nonconformity scores, conformal p-values, and
  accepted/excluded states.
- `taxon_assignment_loci.tsv`: per-locus target-versus-alternative distances,
  margin, placement uncertainty, and interpretation.

The HTML report includes the final assignment while explicitly distinguishing
compatibility p-values from positivity probabilities.

## Validation requirements

Before operational use, freeze the panel, reference cohort, and calibration
artifact, then evaluate an independent sample-level test cohort. Report
sensitivity, specificity, predictive values under stated prevalence
assumptions, and the indeterminate rate. Include target isolates, each important
near neighbor, low-depth samples, mixtures, and lineages excluded from marker
selection. Threshold selection must occur before final held-out evaluation.