# Maintenance scripts

Scripts in this directory are narrow maintenance and validation utilities rather than part
of the public `mlvamaps` command-line interface.

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
  --outdir benchmark_results --threads 2 --trials 3 > runtime.json
```

Use an empty output directory. `--database` is also accepted. Trials alternate
baseline and updated versions in fresh processes. Runtime includes QC, mapping,
inference and output writing but excludes imports. Peak RSS is the maximum of
parent and largest child, not their concurrent total. The utility also reports
the fraction of updated loci attempting micro-assembly. Assembly calls are only
read for comparison after inference. Run this on the server's matched dataset
to assess real-data agreement and throughput.
