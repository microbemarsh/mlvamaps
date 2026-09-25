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

The locus table includes calls missing from either mode, with `comparison` equal
to `exact`, `discordant`, `missing_in_a`, or `missing_in_b`. The JSON reports
`exact_repeat_recovery_by_mode`: exact matches divided by all callable loci in
that mode, including those uncalled in the other mode. For an assembly/SR
comparison, its `assembly` value measures recovery of assembly calls and should
be reviewed alongside shared-call concordance and `missing_call_loci_by_mode`.
Loci uncalled in both modes do not contribute to these metrics. Assembly calls
are comparators, not an independently verified ground truth.

For a server-only dataset, keep the reads and assemblies there and compare
existing outputs with a manifest such as:

```tsv
sample_id	mode	outdir
isolate1	assembly	/path/to/isolate1/assembly
isolate1	sr_before	/path/to/isolate1/sr_before
isolate1	sr_after	/path/to/isolate1/sr_after
```

Use identical reads, panels and settings for the two SR runs, and distinct output
directories. Repeat at the coverage levels relevant to the experiment (for
example 1x, 2x, 5x, 10x and full depth), using the same subsampled molecules for
before/after runs and several subsampling seeds. Use a separate manifest per
depth and seed. Inspect discordant loci and missing calls as well as aggregate
metrics. Compare wall time under the same CPU allocation using
`short_read_run_metadata.json` → `performance.stage_seconds.total`; recruitment
and inference timings are also in `reconstruction_metadata.json`. These metadata
do not report peak RAM; measure that with the server's job accounting if needed.

## Short-read recruitment performance

`benchmark_repeat_fitting.py` compares exhaustive native candidate alignment
against equivalent-score reuse on synthetic flanking reads. It checks exact
equality of all likelihoods, posteriors, mixture fractions and call fields:

```bash
python scripts/benchmark_repeat_fitting.py --molecules 500 --unit 60 \
  --maximum 100 --output repeat_fitting_timing.json
```

This measures repeat fitting only, separately from recruitment and I/O.

`benchmark_fastq_io.py` compares native decompression and compression backends
on the same subset, including QC and retained-read hash checks. Install the
`fastq` Python extra or conda-forge's `rapidgzip` and `python-isal` first:

```bash
python scripts/benchmark_fastq_io.py --reads1 R1.fastq.gz --reads2 R2.fastq.gz \
  --pairs 300000 --threads 8 --scratch-dir "$SLURM_TMPDIR" --output fastq_io_timing.json
```

This isolates input/QC/output costs; sample-wide timings are recorded in
`short_read_run_metadata.json`. QC overlaps recruitment in the default engine.

`benchmark_short_read_recruitment.py` times a bounded subset of paired or
single-end FASTQs at several worker counts and checks evidence equivalence.
Run it in the installed mlvamaps environment and within your CPU allocation:

```bash
python scripts/benchmark_short_read_recruitment.py \
  --reads1 R1.fastq.gz --reads2 R2.fastq.gz --primers panel.tsv \
  --database reference_build --pairs 10000 --threads 1 4 8 \
  --output recruitment_timing.json
```

This isolates recruitment, including process startup and audit hashing. It does
not measure reference classification or end-to-end sample runtime. Production
calls also record stage timings in `reconstruction_metadata.json`.
