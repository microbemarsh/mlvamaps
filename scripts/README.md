# Maintenance scripts

Scripts in this directory are narrow maintenance and validation utilities rather than part
of the public `mlvamaps` command-line interface.

## `benchmark_fastq_inference.py`

Run clean and difficult short-read fixtures against the current checkout without
a baseline or server data:

```bash
python scripts/benchmark_fastq_inference.py --challenge --layouts paired shotgun \
  --workdir /path/to/empty-benchmark-directory --copies 2 --repetitions 3 \
  --output fastq-benchmark.json
```

This writes `benchmark.json` in the work directory and copies that report to
`--output` when supplied. Without `--workdir`, `--challenge` creates and prints a
retained temporary directory's report path. Standalone retained runs also support
ablations and candidate/cell limits. Plain timing runs without these options keep
their original compact output. Challenge estimates remain in the accuracy
denominator; executing the benchmark does not guarantee exact agreement.

Run both source versions in fresh processes on exactly the same retained FASTQs:

```bash
python scripts/benchmark_fastq_inference.py \
  --baseline-source /path/to/baseline-checkout \
  --workdir /path/to/empty-comparison-directory \
  --repetitions 3 --copies 2 --challenge \
  --ablations no_error_model no_motif_refinement no_partial_bounds
```

Use the same Python environment with the development dependencies, minimap2,
and Sassy on `PATH`. Preparation is excluded from timing. The harness alternates
baseline/updated order and runs workers serially. It refuses to overwrite a
nonempty work directory and retains FASTQs, SHA-256 hashes, source hashes,
arguments, outputs, BAMs, logs, measurements, truth comparisons, and every
numeric/status disagreement in `comparison.json`. A failed worker makes the
benchmark fail. Output hashes check repeatability across trials.

Add `--ablations no_short_assembly` to compare the updated mapping-only workflow
on the same files. This invokes `--short-no-assembly` and disables both local
assembly and consensus reconstruction; pair-likelihood rescue requires isolate
mode, the benchmark default.

The clean fixture exercises three layouts and 11 known alleles. With
`--challenge`, `fastq_challenge_fixture.py` adds fixed-seed paired and long-read
fixtures with expansions, contractions, errors, long arrays, ambiguous loci,
mixtures, and missing/low coverage. Their truth is specified before calling;
neither version defines truth. `--layouts paired shotgun --challenge` includes
only the short-read challenge; the simulated long-read challenge requires
selecting `long`. `--copies` changes coverage by copying reads with
new IDs; these are not independent sequencing-error replicates. All loci stay in
the summaries. Stratified metrics and accuracy/call-rate curves retain missing
calls; model-support thresholds are not calibrated probabilities.

The original timing-only invocation remains available. See the
[validation report](../docs/workflows/fastq-evidence-validation.md) for measured
results, limitations, profile reproduction, and the retained disagreements.

## `convert_geonome_metadata.py`

Converts a Geonome metadata table, reference directory, or manifest into the
`reference_metadata.tsv` schema accepted by mlvamaps. It is retained because it
is documented in `docs/reference/input-formats.md` and has dedicated tests.

Run it from the repository root:

```bash
python scripts/convert_geonome_metadata.py /path/to/geonome/metadata.tsv \
  --output reference_metadata.tsv
```

Do not add one-off organism-specific converters here. Prefer generic import
logic in the package, or keep project-specific transformations with the source
dataset so their provenance remains explicit.

## `compare_short_read_results.py`

Compare matching `(sample_id, locus_id)` calls from files or result trees, for
either short or long reads:

```bash
python scripts/compare_short_read_results.py compare \
  --assembly assembly_results --reads updated_results \
  --baseline previous_results > concordance.json
```

The JSON includes paired-call `n`, exact/within-one/within-two agreement, MAE,
median absolute error, tie-corrected Spearman correlation, and per-locus results.
Missing calls are counted separately, never treated as zero; undefined
correlations are null. Duplicate sample/locus rows are rejected. `--short-reads`
remains an alias for `--reads`.

The comparison also reports:

- `sample_matching`: shared and unmatched sample IDs, without silently renaming samples;
- `by_read_status` and `by_estimation_method`: agreement for `PASS`, `LOW_DEPTH`,
  `ESTIMATED`, and other available statuses, and for each estimate method;
- `per_sample`: sample-specific error and missingness summaries;
- `mean_signed_error`: FASTQ minus assembly, distinguishing expansion from contraction bias;
- `disagreements`: ranked numeric discrepancies and missing calls, with source
  paths, statuses, method, and supporting depth.

When sibling diagnostic tables exist, the audit joins `genotype_evidence.tsv`,
`short_read_diagnostics.tsv`, `query_locus_reconstruction.tsv`,
`library_insert_distribution.tsv`, and `sample_summary.tsv` by sample/locus.
Disagreements retain those rows and their source paths. `by_evidence_observation`
groups results by recorded evidence, assembly failure reason, reconstruction
status, and insert-calibration availability/reason. These are diagnostic
observations, not proven causes of a discrepancy. Older results without those
files remain comparable.

`exact_agreement_all_assembly_numeric` includes every locus with an assembly
count, counting missing FASTQ calls as failures. The paired-count metrics also
include `greater_than_2` and `max_absolute_error`. `presence` compares explicit
`yes`/`no` fields separately, including cases where equal numeric predictions
have different presence evidence; `NOT_FOUND` alone does not prove absence.

The all-call result retains estimates. Status strata describe the recorded call
status, not independently validated measurement quality, and must be interpreted
alongside their callable sample/locus counts. Assembly calls are a comparator,
not guaranteed truth. This utility measures paired locus-call agreement; its
Spearman statistic is not a distance-matrix Mantel test.

For a FASTQ-versus-assembly concordance investigation, first inspect sample IDs,
per-locus errors, coverage, and estimate frequency on matched samples with the
same panel and rounding convention. A higher correlation after excluding
estimates is not itself a calling-accuracy improvement: report how many calls
were excluded. Inspect the existing `local_assembly_concordance.tsv` alongside
discordant long-read calls to compare consensus sizing with per-read sizing.

For comparisons with cgMLST or Parsnp, keep sample order and sample membership
identical across matrices. Evaluate repeat-only, within-locus SNP-only, and
combined distances separately using consistent callable loci. cgMLST allele IDs
are categorical labels, so their numeric differences are not allele distances.
The combined exporter normalizes repeats and SNP distances using the available
cohort; differences in missingness and per-locus variability can change their
effective contribution across runs. Inspect component availability and keep
the normalization, masks, and component weights fixed in the evaluation before
attributing a change to genotyping. Do not select weights to maximize correlation
on the same cohort used to report performance.

For ONT reruns, specify `--read-technology ont` (or `ont-hq` for accurate ONT)
and the appropriate sample mode. The CLI now passes that choice to mapping;
older versions accepted the option without applying it. `--sample-mode isolate`
sets isolate reference interpretation; it does not by itself correct noisy
repeat measurements. No server-data improvement is established by the synthetic
regression tests alone.

Time two checkouts on identical paired reads with the same Python environment:

```bash
python scripts/compare_short_read_results.py benchmark \
  --baseline-source /path/to/previous/checkout --source . \
  --panel panel.tsv --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \
  --sample-id sample --assembly assembly_results/sample \
  --sample-mode isolate --assembly-round-tolerance 0.25 \
  --outdir benchmark_results --threads 2 --trials 3 > runtime.json
```

Use an empty output directory. `--database` is also accepted. Trials alternate
baseline and updated versions in fresh processes. Runtime includes QC, mapping,
inference and output writing but excludes imports. Peak RSS is the maximum of
parent and largest child, not their concurrent total. The utility also reports
the fraction of updated loci attempting micro-assembly. Assembly calls are only
read for comparison after inference. Run this on the server's matched dataset
to assess real-data agreement and throughput.
The benchmark passes the same explicit sample mode, rounding, and optional
`--short-insert-mean`/`--short-insert-sd` to both versions and retains available
`fastq_performance.tsv` stage timings. Supply independently calibrated insert
settings together; the comparison script does not estimate them from assembly
repeat counts.
