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

Compare matching `(sample_id, locus_id)` calls from files or result trees:

```bash
python scripts/compare_short_read_results.py compare \
  --assembly assembly_results --short-reads updated_results \
  --baseline previous_results > concordance.json
```

The JSON includes paired-call `n`, exact/within-one/within-two agreement, MAE,
median absolute error, tie-corrected Spearman correlation, and per-locus results.
Missing calls are counted separately, never treated as zero; undefined
correlations are null. Duplicate sample/locus rows are rejected.

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
