# MLVA Seer

MLVA Seer calls VNTR/MLVA loci from a primer table plus either sequencing reads
or an assembly. The main output is a small `calls.tsv` table that says whether
each VNTR was found and, when possible, the repeat count.

## Install

```bash
conda env create -f environment.yml
conda activate mlva-seer
```

For development from this repository:

```bash
python3 -m pip install -e ".[dev]"
pytest
```

## Quick Start

Call directly from reads:

```bash
mlva-seer call primers.tsv sample.fastq.gz
```

Call from an assembly:

```bash
mlva-seer call primers.tsv assembly.fasta
```

Call from an assembly and use reads to add minimap2 depth support:

```bash
mlva-seer call primers.tsv assembly.fasta --reads sample.fastq.gz
```

By default results go in `results/`. The file most users want first is:

```text
results/calls.tsv
```

It includes:

```text
sample_id
locus_id
present
repeat_count
product_size_bp
read_depth
mean_coverage
status
evidence
```

Use `--outdir my_results` to choose a different output folder and `--sample-id`
to set the sample name.

## Primer File

Use a CSV or TSV with these columns:

```text
locus_id
forward_primer
reverse_primer
repeat_unit_length_bp
expected_product_size_bp
nominal_repeat_units
```

The first three columns are required. The repeat/product/unit columns are
strongly recommended because they let MLVA Seer convert product sizes into
repeat counts.

Example:

```text
locus_id	forward_primer	reverse_primer	repeat_unit_length_bp	expected_product_size_bp	nominal_repeat_units
vrrA_12bp_314bp_10U	CACAACTACCACCGATGGCACA	GCGCGTTTCGTTTGATTCATAC	12	314	10
```

Example primer files are in `examples/seer_lab_Ba/`.

## Outputs

For FASTQ input, MLVA Seer also writes detailed evidence tables, filtered reads,
an MLVA fingerprint, optional profile matches, and `report.html`.

For assembly input, MLVA Seer writes extracted primer products to
`assembly_amplicons.tsv` and `assembly_amplicons.fasta`. If `--reads` is supplied,
it maps reads to those products with minimap2 and writes `read_support.tsv`.

## More Commands

Simulate reads for testing:

```bash
mlva-seer simulate \
  --loci examples/mlva_loci.example.tsv \
  --sample-id SIM1 \
  --depth 500 \
  --outdir simulated
```

Export amplirust products from assemblies:

```bash
mlva-seer extract-amplicons \
  --input assembly.fasta \
  --primers examples/seer_lab_Ba/mlva_seer_primers.example.tsv
```

`mlva-nanopore` remains available as a backwards-compatible command name.
