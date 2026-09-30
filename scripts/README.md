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
Loci uncalled in both modes are retained in call-rate denominators and labeled
`missing_in_both`. The table also joins short-read depth, length and evidence
features. Assembly calls
are comparators, not an independently verified ground truth.

To compare directly with assembly + MLVA_finder, add a `mlva_finder` mode row.
For that mode, `outdir` points to MLVA_finder's detailed `*_output.csv`, which
preserves full locus names. `source_sample` selects the exact `strain` field
(usually the assembly filename); `sample_id` pairs it with the SR run. Use the
same primer panel and record matching rounding, primer-error and binning settings.
Example TSV:

```tsv
sample_id	mode	outdir	source_sample
sample_A	mlva_finder	mlva/assemblies_output.csv	sample_A.fasta
sample_A	sr	sr/sample_A
```

```bash
python scripts/benchmark_cross_mode.py runs.tsv --output comparison.json --require-exact
```

The reports are written before the command exits. `--require-exact` returns 1
unless every reported locus has matching called repeats **and calibrated PCR
size in base pairs**, with no missing sample runs. A locus unresolved in both
outputs also fails. Matching rounded alleles cannot hide a one-base size error.
Wrong, missing and provisional SR calls all prevent complete agreement.
Missing values remain distinct from zero. `strict_repeat_match` and
`strict_size_match` record the separate checks; `strict_match` combines them.
The locus table includes `product_size_bp_a/b`, `size_comparison`, and
`size_absolute_difference_bp`. Size checks use `calls.tsv.product_size_bp`,
not the physical sequence length or a depth diagnostic. Older outputs without
product size must be regenerated before strict validation can pass.
The metric
`exact_repeat_recovery_by_mode.mlva_finder` includes all MLVA_finder calls in its
denominator, so dropping difficult loci cannot improve it to 100%.
Use the detailed CSV rather than the shortened `MLVA_analysis_*.csv` table to
avoid ambiguous locus-name matching. Numeric primary-count agreement alone
does not establish agreement on mixture composition or SNPs.

## Synthetic repeat-sizing datasets

```bash
python scripts/make_repeat_sizing_inputs.py --output /tmp/synthetic_sizing
python -m pytest -q tests/test_mlva_finder_concordance.py
```

The generator writes deterministic artificial assemblies, three-column primer
panels, read FASTQs, reference loci, and `truth.tsv` under three sample folders:

- `spanning`: 26 loci, repeat-unit lengths 2–96 bp, zero/whole/partial repeats,
  strict rounding boundaries, forward/reverse primer indels, flank indels,
  reverse orientation, and an 80-copy allele. Reads cover complete products.
- `shotgun`: 11 loci with 440–536 bp products, paired 150-bp reads from 250-bp
  fragments at staggered positions. No fragment contains a complete product.
- `unresolved`: an 80-copy repeat with boundary pairs but no spanning sequence
  or independently calibrated insert distribution. This negative control must
  retain a zero-confidence prior estimate constrained by the observed lower
  bound, and fail strict concordance.

Every reference contains five repeats; all sample alleles differ. `truth.tsv`
distinguishes physical repeat-tract length, actual sequence length, calibrated
PCR size, and unrounded MLVA allele. A flank indel can change the MLVA allele
without changing the physical tandem-repeat count.

The regression runs the public FASTQ caller and assembly caller, comparing all
37 measurable loci against frozen output from the pinned upstream MLVA_finder
script, with default rounding and with rounding disabled. It needs minimap2,
Sassy and the package dependencies, but no network access. See the
[oracle provenance](../tests/data/mlva_finder_oracle/README.md) for regeneration.
These controlled error-free reads test sizing and reconstruction, not accuracy
on noisy, unevenly covered real libraries.

## Server comparisons

For a server-only dataset, keep the reads and assemblies there and compare
existing outputs with a manifest such as:

```tsv
sample_id	mode	outdir
isolate1	assembly	/path/to/isolate1/assembly
isolate1	sr_before	/path/to/isolate1/sr_before
isolate1	sr_after	/path/to/isolate1/sr_after
```

Use identical reads, panels and settings for the two SR runs, and distinct output
directories. Short reads use one reference-backed competitive caller; no engine
selection is required:

```bash
mlvamaps call -i sr --fq1 R1.fastq.gz --fq2 R2.fastq.gz -p primers.tsv --database references \
  -t 8 -o sr_competitive
```

Repeat at the coverage levels relevant to the experiment (for
example 1x, 2x, 5x, 10x and full depth), using the same subsampled molecules for
before/after runs and several subsampling seeds. Use a separate manifest per
depth and seed. Inspect discordant loci and missing calls as well as aggregate
metrics. Compare wall time under the same CPU allocation using
`short_read_run_metadata.json` → `performance.stage_seconds.total`; recruitment
and inference timings are also in `reconstruction_metadata.json`. These metadata
do not report peak RAM; measure that with the server's job accounting if needed.

## Short-read recruitment performance

`benchmark_reconstruction.py` generates six synthetic loci with off-target
background reads, runs the complete short-read path and records wall time and
peak resident memory. Direct and unresolved scenarios exercise separate evidence paths within the
same competitive caller. Each JSON result records its runtime, call correctness,
peak process memory and stage timings:

```bash
python scripts/benchmark_reconstruction.py --molecules 10000 --threads 1 --output bench/direct
python scripts/benchmark_reconstruction.py --molecules 10000 --scenario unresolved --output bench/unresolved
```

These synthetic measurements do not establish performance on real libraries.

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
  --pairs 10000 --threads 1 4 8 \
  --output recruitment_timing.json
```

This isolates recruitment, including process startup and audit hashing. It does
not measure reference classification or end-to-end sample runtime. Production
calls also record stage timings in `reconstruction_metadata.json`.

## Sub-1× validation

```bash
python scripts/benchmark_low_coverage.py --output bench/sub1x --seeds 10
python scripts/benchmark_cross_mode.py bench/sub1x/runs.tsv --output bench/sub1x/comparison.json
```

The first script samples paired 150-bp reads at 0.125×, 0.25×, 0.5× and 0.75×,
with 320±30-bp fragments and a 0.1% substitution rate. Identical reads are run
with rich and calibrated primer-only panels, without a database or supplied
insert-size parameters. Truth is measured by assembly mode on the source
contigs with the same panel used for the corresponding FASTQ run. Every locus is counted, including missing and uncertain calls; the
output separates definitive call errors/recovery from provisional estimates.
These small synthetic controls do not establish superiority over Shovill.

For that comparison, run Shovill plus assembly mode and this caller on identical
subsampled server FASTQs, with independent seeds and equal CPU allocations.
Record Shovill's full command and version: its documented defaults include
contig-coverage filtering, so evaluate both defaults and any explicitly chosen
low-coverage settings ([Shovill documentation](https://github.com/tseemann/shovill)).
Compare both outputs against an independent high-quality truth where available,
using the manifest workflow above; disagreement alone does not establish which
method is correct. Include missing loci and total elapsed time for assembly
plus typing, rather than timing typing alone.
