# Illumina short-read workflow

Illumina mode is a separate evidence model for shotgun or amplicon reads that
often do not span a complete VNTR. It does not reinterpret each mate as a short
long read. Pair identity is retained from FASTQ validation through recruitment,
merging, local assembly, boundary evaluation, and support summaries.

## Command line

```bash
mlvamaps call -p panel.tsv -i sr \
  --fq1 SRR000001_1.fastq.gz \
  --fq2 SRR000001_2.fastq.gz \
  --profiles profiles.tsv \
  --database reference_build \
  --sample-metadata metadata.tsv \
  --sample-id SRR000001 \
  -o results/SRR000001 -t 8
```

For single-end data, omit `--fq2`. Mates are not inferred from filenames.
Interleaved FASTQ is not supported. Compressed files are read without whole-file
decompression. Pair files must have equal record counts and matching normalized
IDs at every record.

## QC

The defaults require 40 post-trim bases and mean Q15. Three-prime trimming is
disabled unless `--short-trim-quality` is set. With the default
`--short-min-pair-retention 0.5`, a good mate remains as an orphan if the other
mate fails. IDs and mate association are not rewritten.

`short_read_qc_summary.tsv` reports input, rejected, retained, and orphan
counts. When enough exact opposite-orientation mappings exist, it also records
the empirical fragment-span median, median absolute deviation, and pair count.

## Competitive recruitment

Recruitment references follow this order:

1. per-locus products from `--database`;
2. rich panel products synthesized from primers, flanks, motif, and expected
   repeat range;
3. unique primer seeds when a complete template is unavailable.

Canonical locus-specific k-mers recruit all loci competitively. A confident
mate can rescue its unseeded mate. Equal-scoring assignments are ambiguous;
confident mates assigned to different loci are discordant. Neither category is
counted as unique support for multiple loci. This conservative seed model is
deterministic and prevents shared repeat-motif k-mers from driving assignment.

## Read merging and local assembly

Overlapping mates are merged only with at least 20 aligned bases and at most 3%
mismatches. A candidate overlap explainable entirely by the tandem motif is
rejected because it cannot establish the fragment span. Conflicts in accepted
overlaps use the higher-quality base. Every locus is assembled
independently with the native SKESA short-read assembler. There is no Python
assembly fallback: Illumina mode exits with a clear error when SKESA is missing
or fails. Install through the supplied Conda environment or select a compatible
executable with `--skesa-bin PATH`.

Up to four locus assemblies run concurrently. SKESA cores are divided from the
global `--threads` budget, so nested jobs do not oversubscribe the requested
CPU count. The local jobs use a small hash table appropriate for recruited
locus reads and report contigs as short as 40 bases. Multiple contigs remain
ambiguous and cannot create an exact call. `short_read_assembly_summary.tsv`
reports the backend, status, contig count and length, depth estimate, and
failure reason. `--keep-intermediates` retains each locus FASTQ, contig FASTA,
and SKESA log under `short_read_assembly_intermediates/`.

## Exact, interval, and presence evidence

Each contig, merged pair, and original read is evaluated with the same panel
anchors and assembly-calibrated repeat convention used elsewhere in mlvamaps.
Evidence is classified as:

- `COMPLETE_ASSEMBLED_PRODUCT`
- `BOUNDARY_SPANNING_READ_PAIR`
- `BOUNDARY_SPANNING_SINGLE_READ`
- `PARTIAL_REPEAT_EVIDENCE`
- `PRESENCE_ONLY`
- `AMBIGUOUS_ASSEMBLY`
- `MULTIPLE_ALLELES`
- `LOW_DEPTH`
- `NOT_FOUND`

An exact `repeat_count` requires a contig, merged pair, or original read that
directly resolves both repeat boundaries. The two boundaries may be the rich
panel flanks or, when flanks are absent, the product primers. A read inside the
repeat or covering one boundary cannot create an exact value.

Opposite boundaries on separate mates may produce `repeat_count_min` and
`repeat_count_max`. Empirical insert size can narrow that interval when at
least two concordant spans support it. The midpoint is never copied into
`repeat_count`. Without an adequate insert estimate, the panel's expected range
is retained and the reason explains why.

## Mixtures and confidence

Allele support counts only molecules with discriminating boundary evidence.
Repeat-internal locus reads remain in `uninformative_locus_reads`. Multiple
defensible alleles are preserved with primary/secondary support, informative
molecules, fractions, and `mixture_status`. Fractions are left empty when no
allele-discriminating molecule exists.

Confidence reasons are textual and auditable. A primer/flank-bounded local
assembly with adequate molecule support is high confidence; direct molecule
evidence without depth is provisional; conflicts lower confidence; interval or
presence-only rows have no falsely precise probability. The HTML report uses a
dedicated Illumina table and labels unresolved rows explicitly.

## Metadata and MYOGA

Metadata can use `sample_id`, `run_accession`, `sra_run`, or `accession` as its
join key. Aliases for BioSample, dates, coordinates, location, country, host,
isolation source, and study are normalized without removing original columns.

`myoga_samples.csv` uses `genome_id = sample_id`. MYOGA recognizes
`genome_id`, `latitude`, `longitude`, `location`, and `collection_date`
directly. If a Newick tree is generated from the same sample IDs, its tip names
must remain exactly those values. In MYOGA, load the Newick/Parsnp tree and
`myoga_samples.csv`, then select `genome_id` as the ID column if it is not
selected automatically. `myoga_loci.csv` is a long-form companion for external
filtering; MYOGA does not require it.

## Manifest batches and HPC

```text
sample_id\treads1\treads2\tmetadata_id
SRR000001\t/path/SRR000001_1.fastq.gz\t/path/SRR000001_2.fastq.gz\tSAMN000001
SRR000002\t/path/SRR000002.fastq.gz\t.\tSAMN000002
```

```bash
mlvamaps call -p panel.tsv -i sr --manifest samples.tsv \
  --sample-metadata metadata.tsv \
  --profiles profiles.tsv --database reference_build \
  -o results -t 32
```

Samples run independently and write to `results/<sample_id>/`. One malformed or
missing sample is recorded in `batch_status.tsv` without stopping the rest.
Successful sample directories resume by default; use `--force` to recompute.
Combined standard, sample-summary, and MYOGA tables are written at the batch
root. Samples are processed sequentially so a process never retains all batch
FASTQs together. Within a sample, FASTQ/QC/recruitment run in bounded chunks;
only uniquely locus-recruited molecules are retained for assembly.

For Slurm arrays, split the manifest by row while preserving its header and run
one manifest shard per task into separate output roots. Merge the resulting TSV
or CSV files after all tasks complete. mlvamaps does not submit scheduler jobs.

## Synthetic worked example

Generate the tiny offline example:

```bash
python examples/make_illumina_example.py examples/illumina_demo
```

Then run:

```bash
mlvamaps call -p examples/illumina_demo/panel.tsv -i sr \
  --fq1 examples/illumina_demo/SRR_DEMO_1.fastq.gz \
  --fq2 examples/illumina_demo/SRR_DEMO_2.fastq.gz \
  --sample-id SRR_DEMO \
  --sample-metadata examples/illumina_demo/metadata.tsv \
  -o examples/illumina_demo/results
```

Open `results/report.html`, inspect the exact-versus-unresolved evidence, and
load `results/myoga_samples.csv` into MYOGA. Validate against an assembly call:

```bash
mlvamaps call -p examples/illumina_demo/panel.tsv \
  -i examples/illumina_demo/truth.fasta.gz -o examples/illumina_demo/truth
mlvamaps validate \
  --truth examples/illumina_demo/truth/calls.tsv \
  --illumina examples/illumina_demo/results/calls.tsv \
  -o examples/illumina_demo/validation
```

## Troubleshooting

- **Different FASTQ counts or IDs:** regenerate mates together; do not sort one
  file independently.
- **Many ambiguous pairs:** provide a reference build with longer, divergent
  locus flanks and review similar loci in the panel.
- **Presence-only locus:** this is expected when neither reads nor the local
  graph resolve both boundaries. Do not replace the blank call with the
  expected-range midpoint.
- **SKESA not found:** create the supplied Conda environment, install the
  `skesa` package, or pass `--skesa-bin PATH`.
- **No contigs:** review retained depth, read length, and off-target repetitive
  evidence. Exact merged reads may still call a locus when SKESA conservatively
  breaks an ambiguous repeat assembly.
- **MYOGA row does not attach to a tip:** make `genome_id` exactly equal to the
  Newick label, including suffixes and case.
