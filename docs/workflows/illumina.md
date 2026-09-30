# Illumina paired-end FASTQ

Short reads use a **single database-free competitive caller**. Supply the
primer panel and FASTQs; there is no short-read engine selector. Native flank
alignments compete across loci, and physical repeat lengths compete using
sample evidence. Supplying `--database` adds downstream reference classification
without changing the genotype algorithm.

There are no minimum depth, supporting-molecule, spanning-pair or graph-edge
count requirements. Single observations remain usable; confidence and intervals
describe the evidence, while unobserved lengths remain unresolved.

```bash
mlvamaps call -i sr --fq1 sample_R1.fastq.gz --fq2 sample_R2.fastq.gz \
  -p primers.tsv -o results -t 8
mlvamaps call -p panel.tsv -i reads/ --short-reads -o results -t 8
```

Mate 2 is optional; competitive recruitment retains both mates and QC-surviving
orphans. Direct spans, unique pair overlaps and local stitching feed the same
assembly-calibrated repeat/SNP interpreter, including novel and half-repeat
alleles and primer-indel correction. Sample-derived arms permit candidate
inference without reference alleles or cached minimap2 indexes. Minimap2 is
used only for optional downstream database classification.

A boundary-only read gives a lower bound, with equal likelihood for all longer
compatible alleles. Repeated copies cannot manufacture a unique length. Inward
pairs contribute a robust fragment likelihood; learned flank-phase offsets are
removed before comparing physical lengths. The caller scores length vectors
without building or aligning every candidate sequence. Local product recovery,
SNP evidence and the read-depth fallback are automatic parts of this pathway.

Legacy three-column primer panels are supported, including names such as
`vrrA_12bp_314bp_10U`. The name supplies repeat-unit length, nominal product
length and nominal repeat count for the shared MLVA length calibration; those
values may also be supplied as panel columns. Complete primer-bounded reads,
uniquely overlapping pairs and uniquely reconstructed products can be called
without knowing the motif or internal flank sequences.
Without product-length calibration, recovered sequences are retained but their
numeric repeat counts remain blank; repeat-unit size alone is insufficient.

Primer-only loci can learn repeat motifs and both boundary flanks from recruited
sample reads. Primer-anchored consensus tolerates isolated sequencing errors;
once both boundaries are learned, partial reads can support likelihood estimates.
Without those boundaries or a complete product, a primer-connected k-mer graph
may still support a provisional read-depth length estimate.
The caller does not treat the entire sequence between primers as a repeat tract.
When the motif is unknown, overlaps must
span at least two configured repeat units and pass an observed periodicity
check. A rich panel with concrete motifs and repeat-boundary flanks additionally
supports partial boundary evidence and the likelihood fallback below.

1. Stream QC and competitive recruitment against a combined flank index.
2. Classify spans, boundaries, flank pairs and anchored repeat-rich reads.
3. Recover observed products and uniquely overlapping pairs.
4. Pile up recruited reads, consolidate substitution errors, and reconstruct
   unresolved loci through unique overlaps.
5. Measure recovered contigs with the same Sassy PCR and length calibration as assembly.
6. Use haploid candidate likelihoods for unresolved loci with panel or sample-learned boundaries.
7. Estimate remaining unresolved lengths from whole-input k-mer coverage when a
   primer-connected graph is available. Interpret recovered products through
   the shared repeat/SNP layer.

The graph retains contiguous high-quality segments rather than discarding a
whole read for one bad base. Overlaps tolerate up to 2% substitutions but must
have a unique non-repeat offset; ambiguous consensus bases remain `N`.
Read and path redundancy is consolidated before applying reconstruction limits,
so thousands of recruited molecules can contribute to a primer-bounded consensus.
Different-length quality-trimmed fragments can join a containing read through a
unique ungapped placement, retaining their original coordinates and base votes.
They contribute no support outside their observed bases. The 2,048-node cap
applies after this consolidation; node/path failures set `limit_reached=true`.
Consensus retains the length established by the overlap coordinates.
Separate recovered components with the same product length and primer ends
also retain their repeat count; disputed interior bases are masked as `N`.
Sassy runs on the small recovered-contig collection after recruitment.
`locus_reconstruction/reconstructed_contigs.fasta.gz` and
`locus_reconstruction/reconstruction_pcr.tsv` retain the sequences and PCR
measurements shown in the report's local assembly section.

The hierarchy is `DIRECT` → `RECONSTRUCTED` → `INFERRED` → `KMER_DEPTH`. Low-confidence
observed products and informative likelihood candidates retain their best
repeat estimate in `calls.tsv`, the fingerprint, and the report, with
`AMBIGUOUS` status and alternatives/intervals. Equally supported lengths use
the smallest candidate; boundary-only estimates may represent only a lower
bound. A flat likelihood or an uncalibrated primer-only product still cannot
supply a numeric count on its own. `KMER_DEPTH` uses `ESTIMATED` status and
uncalibrated confidence (`0`); its interval describes coverage sensitivity and
retains any wider likelihood uncertainty. It supplies no reconstructed sequence
or SNP claims. Substantial incompatible alleles are `MIXED` and retain priority.
See the [method definitions and limitations](../concepts/locus-reconstruction.md).

Useful advanced controls:

```text
--short-repeat-fraction 0.7
--short-confidence-threshold 0.8
--short-max-candidate-repeat-count 100
--sr-recruitment-audit compact
--insert-mean 400 --insert-sd 35
```

Supply insert statistics only when measured independently from the library.
Without reliable insert statistics, boundary-only observations cannot resolve
an exact allele. Depth estimates normalize whole-input graph k-mer abundance
against flank coverage; recruited repeat-rich read counts alone do not determine
copy number. `--threads` controls recruitment, I/O and subsequent locus work
without overlapping full allocations.

Pileups use HTSlib through the installed [`pysam` Python API](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignmentFile.pileup) at validated read
coordinates, without realigning repeat lengths. Overlapping mates contribute
one molecule vote; conflicts become `N`, and consensus bases require 70%
agreement. The read-depth limit is set to the complete input pool, avoiding
HTSlib's default depth truncation. Overlap mismatch counts use NumPy.
Independent A/C/G/T-only pools without shared molecules are batched into one
native `pysam.samtools.consensus` call per consolidation. Pools containing `N`
or overlapping molecular observations use `AlignmentFile.pileup`, preserving
conflict votes and exact molecule weights. Identical observation patterns share
one pileup entry with their full multiplicity; reads are not subsampled.
Independent loci use process workers for reconstruction preparation and final
recovery, capped by the sample's thread allocation and number of loci. Completed
preparations are reused when recruitment has not changed the reads or template.
Sample-anchor learning and indexing share each prepared sequence pool, and
unchanged loci reuse their learned templates between rescue rounds.
When a template changes, reads are reclassified and obsolete boundary
coordinates and length measurements are discarded.
Pileups use temporary uncompressed BAMs in the system temporary directory;
memory scales with the active pileup depth and retained read evidence. Stage
timings, pileup backend and worker allocation are recorded in
`reconstruction_metadata.json`. Sample-anchor learning and indexing currently
run in the main process; progress logs show each learning phase and its locus
molecule count, including after each rescue scan.

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
