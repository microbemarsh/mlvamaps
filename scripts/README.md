# Maintenance scripts

Scripts in this directory are narrow data-conversion utilities rather than part
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
## Cross-mode validation

`benchmark_cross_mode.py` compares existing assembly/SR/LR runs using a TSV
manifest with `sample_id`, `mode`, and `outdir` columns. Relative directories are
resolved against the manifest directory.

```bash
python scripts/benchmark_cross_mode.py runs.tsv --output comparison.json
```

It writes JSON metrics and `comparison.loci.tsv`: exact/within-one repeat
concordance, mean absolute difference, callable loci, full profiles, masked SNP
concordance and mixture component/fraction concordance. Missing calls never
become zero. Optional `--distance-matrix MODE=matrix.tsv` arguments compute tied-rank
Spearman correlation over shared finite off-diagonal distances. SNP comparisons
exclude unknown bases; component identity requires matching full masked sequence.
