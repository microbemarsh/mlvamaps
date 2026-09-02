# Dataset aggregation and MYOGA export

`mlvamaps export-myoga` converts already completed MLVAmaps result directories
into a filtered sample-by-locus dataset, pairwise MLVA distances, and a
MYOGA-ready neighbor-joining tree. It reads existing result files and never
reruns MLVA calling.

```bash
mlvamaps export-myoga \
  --results results/ \
  --metadata sramic_curated_metadata.tsv \
  --metadata-id shared_identifier \
  --latitude latitude \
  --longitude longitude \
  --min-callable-fraction 0 \
  -o global_mlva/
```

CSV and TSV metadata are accepted. A tab in the header takes precedence over
commas elsewhere in a record. `--metadata-id` is matched exactly to the
`sample_id` recorded in `calls.tsv`; the directory basename is not substituted
for a recorded ID. Latitude and longitude aliases used by `--sample-metadata`
are recognized when the default coordinate names are requested.

## Discovery and filtering

The exporter recursively finds `calls.tsv`, prefers per-sample files over the
combined copy in `batch_summary/`, and consults `batch_status.tsv` and
`sample_summary.tsv` when available. Legacy batch roots that store aggregate
files directly at their top level remain supported. Failed, incomplete, duplicate, and
malformed results are recorded in `samples_excluded.tsv` rather than stopping
other samples. Directories containing recognizable MLVAmaps outputs but no
`calls.tsv` are reported as `MISSING_CALLS_FILE`.

A locus is callable only when its final `repeat_count` in `calls.tsv` is a
finite numeric value. `present=yes`, a repeat interval, partial repeat evidence,
or presence-only evidence does not create an exact allele. Missing calls remain
empty in `mlva_profiles.tsv`; they are never converted to zero.

By default, the exporter retains any sample with at least one finite numeric
VNTR `repeat_count`; this is equivalent to `--min-callable-fraction 0` and
`--min-callable-loci 0`. Samples with no exact VNTR calls are still excluded
because no repeat-count distance can be calculated for them.

Use a positive `--min-callable-fraction` when a completeness filter is desired.
Its denominator is the number of loci actually assayed for that sample, as
recorded by its `calls.tsv`; it is not the union of loci from unrelated panels
in the same export. For example, a threshold of `0.8` requires 20 of 25 assayed
loci or 12 of 14 assayed loci. When `--min-callable-loci` is also set, a sample
must satisfy both criteria, equivalently the larger required locus count.

Panel order is recovered from the sample with the most locus rows, with a
deterministic sample-ID tie break. Loci found only in other results are appended
in natural alphanumeric order. This also applies when a batch-level `calls.tsv`
contains a different number of locus rows for each sample. The wide profile uses
the union of loci and leaves loci absent from a sample's input rows empty. Sample
IDs are naturally sorted. Identical MLVA profiles are retained as independent
observations.

## Pairwise MLVA distances

For samples `i` and `j`, `Cij` is the set of loci with finite exact repeat
counts in both samples. The export reports both:

```text
categorical_differences = count(repeat_i != repeat_j over Cij)
categorical_distance    = categorical_differences / |Cij|

repeat_distance_raw = sum(abs(repeat_i - repeat_j) over Cij)
repeat_distance     = repeat_distance_raw / |Cij|
```

The default tree metric is `repeat`, the mean absolute repeat-count difference.
Use `--distance categorical` for the proportion of differing shared alleles.
The existing combined-marker normalization in MLVAmaps scales query/reference
repeat differences by dispersion in a supplied reference database. That scale
cannot be recovered reliably from ordinary completed sample results, so this
export does not invent or silently reuse a dataset-derived replacement in its
default repeat-only mode. The explicitly requested combined-marker mode below
uses a documented dataset-derived scale.

A pair is supported by default when it shares at least one exact call
(`--min-pairwise-loci 1`); the default `--min-pairwise-fraction 0` adds no
fractional completeness requirement. The fractional denominator, when a
positive threshold is requested, is the intersection of the two samples'
assayed loci, so mixing compatible panels does not count loci absent from one
panel as failed calls. Unsupported pairs remain in
`mlva_pairwise_distances.tsv` with `comparison_status=insufficient_overlap` and
empty normalized distances. The exporter never substitutes zero or a maximum
distance. Before tree construction it repeatedly removes the sample with the
most unsupported relationships, then the fewest callable loci and natural
sample order as deterministic tie breaks, until the retained distance matrix is
complete. These removals are reported as `INSUFFICIENT_PAIRWISE_OVERLAP`.

## Optional per-locus SNP and combined-marker trees

Use the opt-in mode when retained sequence artifacts are available:

```bash
mlvamaps export-myoga \
  --results results/ \
  --metadata metadata.tsv \
  --combined-markers \
  --loci panel.tsv \
  --phylogeny-snp-weight 1 \
  --phylogeny-repeat-weight 1 \
  -t 32 \
  -o global_mlva/
```

`--loci` is unnecessary when every usable locus already has the repeat-masked
`phylogeny/LOCUS/query.fasta.gz` written by a prior `--database` analysis. The
exporter prefers that exact precomputed SNP sequence. Otherwise it uses the
selected `calls.tsv:evidence` product from `assembly_amplicons.fasta.gz`, or a
unique compatible product from `local_assembly_pcr/matches.tsv` or
`local_locus_products.fasta.gz`, and masks it with the original rich panel.
The exact assembly evidence identifier is authoritative; fallback products
must agree with the final called product size and be sequence-unique. A
sequence whose VNTR tract cannot be bounded is rejected rather than
allowing repeat-length gaps into the SNP model. Every decision is recorded in
`combined_marker_sequence_status.tsv`.

Exact duplicate SNP sequences are collapsed before inference. A locus with at
least four SNP haplotypes is aligned with MAFFT and inferred with RAxML-NG;
two or three haplotypes use their aligned pairwise SNP distances, and an
invariant locus has SNP distance zero. The resulting haplotype distance is
expanded back to all sample tips, so duplicate samples remain present at zero
SNP distance in `locus_trees/LOCUS/samples.tree`.

For locus `l`, the retrospective export uses:

```text
snp_scale_l = median positive patristic distance among SNP haplotypes
snp_component_ijl = patristic_distance_ijl / snp_scale_l

repeat_scale_l = max(population standard deviation of accepted repeat counts, 0.5)
repeat_component_ijl = abs(repeat_il - repeat_jl) / repeat_scale_l

combined_distance_ij = mean over shared accepted loci of
    (snp_weight * snp_component_ijl + repeat_weight * repeat_component_ijl)
```

The scale falls back to `1` when a locus has no positive SNP distance. A mean,
not a sum, prevents pairs with more shared loci from becoming artificially
more distant. The same pairwise overlap thresholds used by the repeat-only
export apply. Unsupported relationships are left missing, and samples are
removed deterministically until the combined matrix is complete. Use
`combined_marker_nj.tree` with `combined_marker_metadata.tsv`; their sample IDs
match exactly.

## Metadata and geography

`myoga_metadata.tsv` has exactly one row per final tree tip and a canonical
`sample_id` equal to that tip label. Unmatched samples remain in the MLVA tree
with blank metadata. Original metadata columns are retained; input columns that
would conflict with canonical `sample_id`, `latitude`, or `longitude` are
preserved with a `metadata_` prefix.

Coordinates must be finite and satisfy `-90 <= latitude <= 90` and
`-180 <= longitude <= 180`. Missing metadata, missing coordinates, and invalid
coordinates are recorded in `samples_excluded.tsv` with `scope=geography`, but
do not remove an otherwise usable MLVA sample. Invalid source coordinate text
is preserved in its `metadata_` column while the canonical coordinate is left
blank. Tree exclusions have `scope=tree`.

## Outputs

```text
global_mlva/
├── myoga_metadata.tsv
├── mlva_profiles.tsv
├── mlva_calls_long.tsv
├── mlva_pairwise_distances.tsv
├── mlva_distance_matrix.tsv
├── mlva_nj.tree
├── samples_used.tsv
├── samples_excluded.tsv
├── export_summary.tsv
└── export_summary.txt
```

With `--combined-markers`, the directory also contains:

```text
locus_trees/LOCUS/
    samples.repeat_masked.aligned.fasta
    haplotypes.raxml.tree        # loci with >=3 SNP haplotypes
    samples.tree
combined_marker_sequence_status.tsv
locus_tree_status.tsv
locus_snp_distances.tsv
combined_marker_pairwise_distances.tsv
combined_marker_distance_matrix.tsv
combined_marker_nj.tree
combined_marker_metadata.tsv
```

`mlva_calls_long.tsv` and `mlva_profiles.tsv` contain the final tree samples.
The pairwise long table also retains threshold-passing samples subsequently
removed for insufficient overlap, making those decisions auditable. The square
matrix contains only final samples in identical row and column order.

`mlva_nj.tree` is built with MLVAmaps' deterministic NumPy neighbor-joining
implementation. One- and two-sample exports produce valid simple Newick trees;
no tree is written for zero retained samples. Negative reconstruction branch
lengths are clamped to zero. Supply `mlva_nj.tree` and `myoga_metadata.tsv` to
MYOGA. Their tip/sample-ID sets are identical, including samples whose
coordinates are blank.

This is an **MLVA relatedness tree**, not a whole-genome phylogeny or a
nucleotide-substitution model. Repeat homoplasy, locus-specific mutation rates,
missing-data overlap, and the small number of MLVA loci limit evolutionary
interpretation.

The combined tree is likewise a repeat-aware multilocus marker relatedness
tree, not a whole-genome phylogeny. Its dataset-derived normalization makes it
appropriate for relationships within this export, but distances should not be
compared numerically across independently normalized exports.

Existing batch `myoga_samples.csv` and `myoga_loci.csv` outputs remain
unchanged for backward compatibility. The retrospective exporter is separate
and can be run after any completed batch. Existing export files are protected;
pass `--force` to replace them deterministically.
