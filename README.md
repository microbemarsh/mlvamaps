# mlva-seer

mlva-seer calls VNTR/MLVA loci from a primer table plus either sequencing reads
or an assembly. The main output is a small `calls.tsv` table that says whether
each VNTR was found and, when possible, the repeat count.

<p align="center">
  <img width="500" height="400" alt="image" src="https://github.com/user-attachments/assets/9d6b1032-d4c2-40c8-9bfd-7c93e4f496dc" />
</p>

## Install

The recommended install method is through conda. Please make sure to have a working conda installation first. If you don't, please select the appropriate version found [here](https://github.com/conda-forge/miniforge).

A bioconda release will happen in the future, but until then, install using this:

```bash
# clone repo
git clone https://github.com/microbemarsh/mlva_seer.git

# cd to repo
cd mlva_seer

# make conda environment
conda env create -f environment.yml

# activate conda environment
conda activate mlva-seer

# install mlva-seer without invoking the PEP 517/pyproject build path
python setup.py install
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

Or use an existing BAM/SAM aligned to the assembly:

```bash
mlva-seer call primers.tsv assembly.fasta --bam assembly_reads.bam
```

For read-depth mapping, minimap2 uses its default alignment settings unless you
choose a preset with `--minimap2-preset`, for example `--minimap2-preset sr`.

To compare a sample against known MLVA types, add a profile database:

```bash
mlva-seer call primers.tsv assembly.fasta --profiles mlva_profiles.tsv
```

This repository includes a converted UF B. anthracis profile database at:

```text
data/uf_ba_mlva_profiles.tsv
```

To rebuild it from the original UF table:

```bash
python scripts/convert_uf_ba_vntrs.py /path/to/uf_ba_vntrs.tsv
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

Use `-o my_results` (or `--output my_results`) to choose a different output
folder and `--sample-id`
to set the sample name.

MLVA Seer prints live progress updates while it runs. Use `--quiet` to suppress
them in scripts. It uses 32 worker threads by default; pass a different number
with `-t`/`--threads`, or explicitly use `--threads 0` to use all available CPU
cores.

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
an MLVA fingerprint, optional profile matches, and `report.html`. VNTR reads
are clustered by Savont from quality-bearing, primer-trimmed amplicons. Degenerate
primer pairing and read-to-locus assignment are performed by Amplirust using a
lossless FASTA projection of the filtered reads; the original FASTQ qualities
remain attached to the assigned reads. MLVA Seer
writes one FASTQ per locus and submits them together as a single pooled Savont
run, so Savont owns the complete `-t`/`--threads` allocation without nested
Python chunking. Only Savont's retained cluster assignments are consumed.
MLVA Seer independently calculates support, frequency, repeat count, and edit
evidence. Each cluster is represented by its most-supported observed repeat
sequence, using mean read quality and sequence as deterministic tie breakers.
Savont's synthesized consensus is never used as the reported ASV sequence, so
the selected representative's actual indel state is preserved.
`vntr_asv_representatives.fasta` contains those observed sequences.
`vntr_asv_table.tsv` contains per-ASV read and indel counts, while
`vntr_asv_memberships.tsv` retains each clustered read's raw/aligned repeat
sequence and its insertions, deletions, and substitutions relative to the
selected representative. Raw Savont outputs are retained under `results/savont/`
for diagnostics but are not used for downstream abundance or sequence analysis.
Savont 0.6.1 or newer is required. ASVs require two supporting reads by
default; change this with `--savont-min-cluster-size`.

For assembly input, MLVA Seer writes extracted primer products to
`assembly_amplicons.tsv` and `assembly_amplicons.fasta`. If `--reads` is supplied,
it maps reads to those products with minimap2 and writes `read_support.tsv`.
If `--bam` is supplied, read depth is estimated from the assembly-aligned
BAM/SAM. The assembly `report.html` includes a generated gel electrophoresis
image; band position follows product size and band intensity follows read-depth
support when FASTQ or BAM/SAM evidence is available. Without depth evidence,
present loci are drawn with a default band intensity. When profile matching is
available, the closest reference profile is drawn as a neighboring gel lane, and
band labels are shown on hover to keep the gel readable.

If `--profiles` is supplied, both FASTQ and assembly runs write
`mlva_fingerprint.tsv`, `profile_matches.tsv`, and closest-profile summaries in
`report.html`. The profile table should include `profile_id`, optional
`strain_id`, and one column per VNTR locus.

## More Commands

Simulate reads for testing:

```bash
mlva-seer simulate \
  --loci examples/mlva_loci.example.tsv \
  --sample-id SIM1 \
  --depth 500 \
  -o simulated
```

Export amplirust products from assemblies:

```bash
mlva-seer extract-amplicons \
  --input assembly.fasta \
  --primers examples/seer_lab_Ba/mlva_seer_primers.example.tsv
```

## Motivation and recognition

mlva-seer was designed to speed up the current MLVA typing process in the Seer Lab at the University of Florida. It was heavily influenced by [MLVA_finder](https://github.com/i2bc/MLVA_finder) and [amplirust](https://github.com/erdikilic/amplirust), the latter being a dependency. 
