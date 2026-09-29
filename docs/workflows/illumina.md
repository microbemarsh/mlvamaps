# Illumina paired-end FASTQ

Short reads use **competitive locus recruitment and local VNTR reconstruction**.
There is one standard short-read pathway.

```bash
mlvamaps call -i sr --fq1 sample_R1.fastq.gz --fq2 sample_R2.fastq.gz \
  -p primers.tsv -o results -t 8
mlvamaps call -p panel.tsv -i reads/ --short-reads -o results -t 8
```

No reference database is needed. Recruitment and reconstruction use only the
supplied panel and observed reads, even when `--database` is supplied for later
reference classification. Mate 2 is optional; recruited pairs retain both mates.

Legacy three-column primer panels are supported, including names such as
`vrrA_12bp_314bp_10U`. The name supplies repeat-unit length, nominal product
length and nominal repeat count for the shared MLVA length calibration; those
values may also be supplied as panel columns. Complete primer-bounded reads,
uniquely overlapping pairs and uniquely reconstructed products can be called
without knowing the motif or internal flank sequences.
Without product-length calibration, recovered sequences are retained but their
numeric repeat counts remain blank; repeat-unit size alone is insufficient.

Primer-only loci remain unresolved when a complete product cannot be recovered.
They do not use synthetic repeat candidates or treat the entire sequence
between primers as a repeat tract. When the motif is unknown, overlaps must
span at least two configured repeat units and pass an observed periodicity
check. A rich panel with concrete motifs and repeat-boundary flanks additionally
supports partial boundary evidence and the likelihood fallback below.

1. Stream QC and competitive recruitment against a combined flank index.
2. Classify spans, boundaries, flank pairs and anchored repeat-rich reads.
3. Recover observed products and uniquely overlapping pairs.
4. Reconstruct unresolved loci with a bounded overlap graph and SPOARS consensus.
5. Measure recovered contigs with the same Sassy PCR and length calibration as assembly.
6. Use haploid candidate likelihoods only for unresolved rich-panel loci, then
   interpret products through the shared repeat/SNP layer.

The graph retains contiguous high-quality segments rather than discarding a
whole read for one bad base. Overlaps tolerate up to 2% substitutions but must
have a unique non-repeat offset; ambiguous consensus bases remain `N`.
Sassy runs on the small recovered-contig collection after recruitment.
`locus_reconstruction/reconstructed_contigs.fasta.gz` and
`locus_reconstruction/reconstruction_pcr.tsv` retain the sequences and PCR
measurements shown in the report's local assembly section.

The hierarchy is `DIRECT` → `RECONSTRUCTED` → `INFERRED`. Low-confidence
observed products and informative likelihood candidates retain their best
repeat estimate in `calls.tsv`, the fingerprint, and the report, with
`AMBIGUOUS` status and alternatives/intervals. Equally supported lengths use
the smallest candidate; boundary-only estimates may represent only a lower
bound. A flat likelihood or an uncalibrated primer-only product still cannot
supply a numeric count. Substantial incompatible alleles are `MIXED`.
See the [method definitions and limitations](../concepts/locus-reconstruction.md).

Useful advanced controls:

```text
--short-repeat-fraction 0.7
--short-confidence-threshold 0.8
--short-min-spanning-pairs 2
--short-max-candidate-repeat-count 100
--sr-recruitment-audit compact
--insert-mean 400 --insert-sd 35
```

Supply insert statistics only when measured independently from the library.
Without reliable insert statistics, boundary-only observations cannot resolve
an exact allele. Repeat-rich read abundance is not treated as a copy-number
measurement. `--threads` controls recruitment, I/O and subsequent locus work
without overlapping full allocations.

`calls.tsv` preserves the primary schema and adds call method and confidence.
`short_read_repeat_evidence.tsv` contains evidence counts and failure-mode
features; `reconstructed_locus_variants.tsv` retains recovered sequence.
Candidate likelihood files are empty for loci resolved without inference.
Run metadata includes stage timings and insert statistics.

## Validation commands

```bash
python -m pytest -q
python scripts/benchmark_reconstruction.py --output bench/direct --molecules 10000 --threads 1
python scripts/benchmark_reconstruction.py --output bench/unresolved --molecules 10000 --scenario unresolved --threads 1
python scripts/benchmark_short_read_recruitment.py --help
python scripts/benchmark_cross_mode.py manifest.tsv --output comparison.json
```

The cross-mode manifest has `sample_id`, `mode` (`assembly`, `sr`, `lr`) and
`outdir`. Run assembly and FASTQ independently with the same panel. Assembly is
a validation target and is never supplied as training labels or expected calls.
The comparison reports exact/±1 concordance, dropout, incorrect-call rate,
per-locus call rates, SNP concordance and available short-read QC features.
It also reports `best_estimates_including_ambiguous` concordance separately
from confident-call metrics. The companion `.loci.tsv` includes both best
estimates, statuses, and their absolute differences.
Synthetic performance fixtures are not a substitute for paired real assemblies
and Illumina libraries. Resource measurements are described in the
[refactor validation report](../reference/short-read-validation.md).

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
  --profiles profiles.tsv \
  -o results -t 32
```

Samples run independently and write to `results/<sample_id>/`. One malformed or
missing sample is recorded in `results/batch_summary/batch_status.tsv` without
stopping the rest.
Successful sample directories resume by default; use `--force` to recompute.
Combined standard, sample-summary, and MYOGA tables are written at the batch
root's clearly scoped `batch_summary/` directory.

Independent samples run concurrently under one global `--threads` budget. The
number of active samples is the minimum of sample count, available threads, and
a memory-aware cap of four. Each sample receives an equal integer share of the
budget, and the allocated total never exceeds `--threads`. Set a lower cap for
large samples or memory-constrained nodes:

```bash
export MLVAMAPS_MAX_CONCURRENT_SAMPLES=2
mlvamaps call -p panel.tsv -i short_read_directory/ --short-reads \
  -o results -t 32
```

Completion order does not change combined output order. Within each sample,
FASTQ/QC streams in bounded chunks and process workers recruit molecules
against panel anchors. Progress reports sample/worker allocation and stage
transitions unless `--quiet` is selected.

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
load `results/myoga_samples.csv` into MYOGA. An assembly can be called
independently for scientific concordance assessment:

```bash
mlvamaps call -p examples/illumina_demo/panel.tsv \
  -i examples/illumina_demo/truth.fasta.gz -o examples/illumina_demo/truth
```

## Troubleshooting

- **Different FASTQ counts or IDs:** regenerate mates together; do not sort one
  file independently.
- **Many ambiguous pairs:** review similar primers in the panel. Where known,
  provide longer, divergent locus flanks in a rich panel.
- **Presence-only locus:** examine the evidence reason. Informative partial
  reads retain a provisional estimate or lower bound; flank-only reads without
  usable fragment-length information cannot estimate repeat number. An
  incomplete primer-only product also needs reconstruction or additional
  repeat-boundary information. Expected-range midpoints are not measurements.
- **Database predates the context schema:** omit `--database` for allele calling;
  rebuild it only if reference classification is wanted.
- **MYOGA row does not attach to a tip:** make `genome_id` exactly equal to the
  Newick label, including suffixes and case.
