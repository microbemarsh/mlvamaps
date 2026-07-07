# MLVA Seer

MLVA Seer was designed to meet the needs of the SEER Lab at the University of Florida.
mlva-seer is a software built for Nanopore amplicon-based MLVA/VNTR
genotyping. It is designed to modernize older assembly-based in-silico PCR
workflows, such as `i2bc/MLVA_finder`, by calling VNTR repeat copy number
directly from fastq reads while retaining interpretable evidence
tables for microbiology and outbreak surveillance teams.

This first implementation performs primer-based read assignment, flank-aware repeat
array extraction, motif-aware repeat parsing, simple VNTR repeat variant
clustering, read-level probabilistic repeat-count prediction, Bayesian-style
per-locus aggregation, MLVA fingerprint construction, known-profile matching,
novelty scoring, and HTML reporting.

For any assembly/contig/reference in-silico PCR step, MLVA Seer delegates to
[amplirust](https://github.com/erdikilic/amplirust) rather than maintaining a
separate primer-search implementation. amplirust provides fast primer matching,
FASTA/GenBank input support, circular genome handling, reverse-complement
search, mismatch/edit-distance settings, product FASTA output, and detailed TSV
statistics.

## Install

For conda (preferred), including the amplirust executable used by the
assembly/reference extraction command:

```bash
conda env create -f environment.yml
conda activate mlva-seer
```

From the repository root:

```bash
python3 -m pip install -e .
```

For tests:

```bash
python3 -m pip install -e ".[dev]"
pytest
```

## Primer Input

The simplest required analysis inputs are:

```text
--input assembly.fasta, contigs.fasta, or reads.fastq.gz
--primers primers.tsv
```

The recommended primer table is tab-delimited:

```text
locus_id
forward_primer
reverse_primer
repeat_unit_length_bp
expected_product_size_bp
nominal_repeat_units
```

`examples/mlva_seer_primers.example.tsv` is a cleaned, headered copy derived
from the legacy MLVA_finder primer file:

```text
vrrA_12bp_314bp_10U  CACAACTACCACCGATGGCACA  GCGCGTTTCGTTTGATTCATAC  12  314  10
```

The tool also accepts the raw legacy whitespace-delimited file shape:

```text
locus_id forward_primer reverse_primer
```

Forward primer is listed first, followed by reverse primer. Reverse primers
should be supplied as the primer oligo sequence in 5' to 3' orientation; the
caller searches for its reverse complement in the oriented amplicon read.

## Full Loci File

For richer assays, `mlva_loci.tsv` is tab-delimited:

```text
locus_id
chrom_or_contig
start
end
forward_primer
reverse_primer
left_flank_sequence
right_flank_sequence
repeat_motif
expected_min_repeats
expected_max_repeats
expected_amplicon_min_bp
expected_amplicon_max_bp
pool_id
```

`left_flank_sequence`, `right_flank_sequence`, and `pool_id` may be blank.

Example files are available in `examples/`, including the recommended
`examples/mlva_seer_primers.example.tsv`, the raw
`examples/insilicoMLVAprimers_all.raw.example.csv`, and the richer
`examples/mlva_loci.example.tsv`.

## Call From FASTQ

```bash
mlva-nanopore call \
  --input sample.fastq.gz \
  --primers examples/mlva_seer_primers.example.tsv
```

Outputs:

```text
results/
  qc_summary.tsv
  filtered_reads.fastq.gz
  read_locus_assignments.tsv
  read_repeat_features.tsv
  vntr_asv_table.tsv
  vntr_asv_consensus.fasta
  read_level_allele_predictions.tsv
  allele_calls.tsv
  mlva_fingerprint.tsv
  mlva_fingerprint_probabilistic.tsv
  profile_matches.tsv
  novelty_scores.tsv
  report.html
```

## Simulate Reads

```bash
mlva-nanopore simulate \
  --loci mlva_loci.tsv \
  --profile mlva_profiles.tsv \
  --profile-id P1 \
  --sample-id SIM1 \
  --depth 500 \
  --error-rate 0.03 \
  --outdir simulated/
```

This writes `SIM1.fastq.gz` and `truth_profile.tsv`, which can be passed into
`mlva-nanopore call`.

## Extract Amplicons From Assemblies

Install amplirust first if not using mlva-seer conda install.

```bash
conda install bioconda::amplirust
```

Then use the primer table to build an amplirust primer CSV and extract
products from FASTA or GenBank inputs:

```bash
mlva-nanopore extract-amplicons \
  --input assembly.fasta \
  --primers examples/mlva_seer_primers.example.tsv
```

Outputs:

```text
assembly_amplicons/
  amplirust_primers.csv
  amplirust_products.fasta
  amplirust_stats.tsv
```

Use `--no-search-rc` to disable reverse-complement searching, `--trim-primers`
to remove primer sequences from extracted products, and `--threads` to control
parallelism.

## Profile Database

The optional `mlva_profiles.tsv` should include `profile_id`, optional
`strain_id`, one column per locus, and optional metadata:

```text
profile_id	strain_id	VNTR_01	VNTR_02	VNTR_03	metadata
P1	        strain_A	5	    4	    8	    outbreak reference
```

The matcher reports exact or nearest known profiles using missing-locus-tolerant
Manhattan distance.

## Current Scope

This is a working research prototype, not yet a validated clinical or regulatory
assay. The current caller is direct-FASTQ and primer/flank driven, while
assembly/reference product extraction is routed through amplirust. Planned
extensions include feeding amplirust products into the same VNTR parser,
minimap2-backed read-to-amplicon assignment, polished consensus calling,
trained read-level models, mixture detection, calibration plots, and
organism-specific binning rules that mirror published MLVA allele conventions.
