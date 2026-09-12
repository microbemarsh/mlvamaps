# Dataset aggregation and MYOGA export

`mlvamaps export-myoga` converts already completed `mlvamaps` result directories
into a filtered sample-by-locus dataset, pairwise MLVA distances, and a
MYOGA-ready neighbor-joining tree. It reads existing result files and never
reruns MLVA calling.

Add `--combined-markers` to restore the retrospective repeat-masked SNP plus
repeat-distance workflow. It uses eight threads by default (`-t/--threads`),
retains every discovered valid sample in the combined distance matrix and
metadata, and leaves a distance blank only when two samples have no recoverable
shared marker. Supply the rich locus panel with `--loci` when retained amplicons
must be masked because precomputed masked query sequences are unavailable. For
current assembly results, the exporter reads
`classification/query_amplicons.fasta`; partial successful calls
(`run_status=success_partial`) are valid export inputs and are not discarded.

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
Directories containing recognizable `mlvamaps` outputs but no
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
A pair is supported by default when it shares at least one exact call
(`--min-pairwise-loci 1`); the default `--min-pairwise-fraction 0` adds no
fractional completeness requirement. The fractional denominator, when a
positive threshold is requested, is the intersection of the two samples'
assayed loci, so mixing compatible panels does not count loci absent from one
panel as failed calls. Unsupported pairs remain in
`mlva_pairwise_distances.tsv` with
`comparison_status=insufficient_overlap`. Their calculated distance is retained
when the pair has at least one shared exact call; it is empty only when no
distance can be calculated. The complete square distance matrix retains every
sample that passes the per-sample callable threshold, including blank cells for
pairs with no shared exact calls. The exporter never substitutes zero or a
maximum distance. Before tree construction it repeatedly removes the sample
with the
most unsupported relationships, then the fewest callable loci and natural
sample order as deterministic tie breaks, until the retained distance matrix is
complete. These removals are reported as `INSUFFICIENT_PAIRWISE_OVERLAP`.

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

`mlva_calls_long.tsv` and `mlva_profiles.tsv` contain the final tree samples.
The pairwise long table also retains threshold-passing samples subsequently
removed for insufficient overlap, making those decisions auditable. The square
matrix contains all threshold-passing samples in identical row and column order;
blank cells mean that a pair has no shared exact call.

`mlva_nj.tree` is built with `mlvamaps`' deterministic NumPy neighbor-joining
implementation. One- and two-sample exports produce valid simple Newick trees;
no tree is written for zero retained samples. Negative reconstruction branch
lengths are clamped to zero. Supply `mlva_nj.tree` and `myoga_metadata.tsv` to
MYOGA. Their tip/sample-ID sets are identical, including samples whose
coordinates are blank.

This is an **MLVA relatedness tree**, not a whole-genome phylogeny or a
nucleotide-substitution model. Repeat homoplasy, locus-specific mutation rates,
missing-data overlap, and the small number of MLVA loci limit evolutionary
interpretation.

Existing batch `myoga_samples.csv` and `myoga_loci.csv` outputs remain
unchanged for backward compatibility. The retrospective exporter is separate
and can be run after any completed batch. Existing export files are protected;
pass `--force` to replace them deterministically.
