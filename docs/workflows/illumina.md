# Illumina paired-end FASTQ

SR FASTQ calling requires **`--database`**. Assembly and default accurate
long-read calling remain database independent. Taxonomic classification always
requires a database and runs separately from repeat measurement.

```bash
mlvamaps call -i sr --fq1 sample_R1.fastq.gz --fq2 sample_R2.fastq.gz \
  --database references -p primers.tsv -o results -t 8
mlvamaps call -p panel.tsv -i reads/ --short-reads \
  --database references -o results -t 8
```

A reference build can supply its saved panel when `-p` is omitted. Mate 2 is
optional. QC retains usable orphan mates, and recruitment preserves both mates
when only one maps. Minimap2 is required. There is no short-read engine selector.

The calling sequence is:

1. Filter and validate FASTQs.
2. Competitively map reads against the database's candidate contexts for all
   panel loci. Duplicate references and allele expansions do not multiply
   molecule support. Cross-locus ties remain unassigned.
3. Classify reads against supported reference flanks. Measure complete products
   and reads spanning both repeat boundaries; reconstruct unique overlaps where
   useful. Unobserved bases remain `N`.
4. For unresolved lengths, use inward mate geometry with empirical or supplied
   insert statistics. A single boundary supplies a lower bound, never an exact
   count. Equally supported lengths remain unresolved.
5. Estimate remaining lengths from whole-input k-mer coverage only when an
   observed, primer-connected graph and usable flank depth exist.
6. Apply shared assembly calibration and Sassy PCR to recovered complete
   products, then write repeat/SNP results and database classification.

The database identifies locus context; its stored repeat count is not a sample
call. Reads can establish alleles outside the reference candidate grid, including
half-repeat alleles. Reference-relative lengths assume the unobserved flanks
have the reference lengths; validate these assumptions on real libraries.
Up to eight equally supported reference backgrounds are evaluated. Conflicting
lengths and larger background ties remain unresolved.

Legacy three-column panels are supported. Calibrated locus names or panel
columns supply repeat-unit and product-length calibration. Reference contexts
provide the internal flanks missing from primer-only panels. A database lacking
a locus cannot supply a call for that locus.

There are no minimum supporting-molecule or depth cutoffs. Confidence, support,
intervals and failure reasons remain explicit. The hierarchy is `DIRECT` →
`RECONSTRUCTED` → `INFERRED` → `KMER_DEPTH`, with `AMBIGUOUS`, `MIXED` and no-call
outcomes. `KMER_DEPTH` is a provisional `ESTIMATED` result, with confidence zero
and a coverage-sensitivity interval; it does not supply reconstructed sequence.

`reference_calling/summary.tsv` records outcomes and reference IDs, alongside
candidate contexts and provenance. BAMs are retained with `--keep-intermediates`.
Calls include `reference_assisted` provenance. `short_read_repeat_evidence.tsv`
records evidence counts, product/VNTR lengths and uncertainty;
`reconstructed_locus_variants.tsv` retains sample-derived sequence.
`short_read_run_metadata.json` identifies the caller as
`competitive_reference_likelihood`.

## Performance and CPU allocation

Reference mapping uses the sample's `--threads` budget and writes BAM through
HTSlib. Normal runs reuse the QC FASTQs. A combined replay retains mapped pairs
and orphan mates; input is not replayed separately for each locus. Iterator-only
API inputs use temporary FASTQs. Reference fitting currently runs per locus in
the main process. Memory scales with recruited evidence and the reference
candidate count. Optional rapidgzip/ISA-L acceleration remains available.

Use `--force` when comparing previously completed server runs. Historical
sample-only benchmarks exercise lower-level recovery helpers, which remain
available for synthetic validation; they do not measure the database-backed
SR command's total runtime.

## Validation commands

```bash
python -m pytest -q
python scripts/benchmark_reconstruction.py --output bench/direct --molecules 10000 --threads 1
python scripts/benchmark_reconstruction.py --output bench/unresolved --molecules 10000 --scenario unresolved --threads 1
python scripts/benchmark_reconstruction.py --output bench/pileup --molecules 12000 --scenario pileup --threads 4
python scripts/benchmark_reconstruction.py --output bench/depth --molecules 12000 --scenario depth --threads 2
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
mlvamaps call -p panel.tsv -i sr --manifest samples.tsv --database references \
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
mlvamaps call -p panel.tsv -i short_read_directory/ --short-reads --database references \
  -o results -t 32
```

Completion order does not change combined output order. Within each sample,
FASTQ/QC streams in bounded chunks and process workers recruit molecules
against reference targets. Progress reports sample/worker allocation and stage
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
mlvamaps call -p examples/illumina_demo/panel.tsv -i sr --database references \
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
- **Database predates the context schema:** rebuild it with the current
  `mlvamaps build-reference`; SR FASTQ calling requires current reference assets.
- **MYOGA row does not attach to a tip:** make `genome_id` exactly equal to the
  Newick label, including suffixes and case.
