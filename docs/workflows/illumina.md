# Illumina short-read workflow

The default `--sr-engine repeat-likelihood` recruits molecules by non-repeat
flanks, then combines E (enclosing), S (spanning pair), F (flanking), and anchored
FRR evidence. Repeat states are synthetic and can be novel or outside the panel's
expected range. Haploid likelihood and the existing microbial EM model retain
distinct components, including supported SNP haplotypes with the same count.

```bash
mlvamaps call --loci panel.tsv -i sr --fq1 reads_1.fastq.gz --fq2 reads_2.fastq.gz \
  --insert-mean 400 --insert-sd 35 --sample-mode metagenome -o results
```

Insert overrides must describe independently measured library statistics. They
are optional when repeat-independent flank pairs can estimate them or enclosing
reads determine the repeat directly. Single-end data can use E and F evidence;
FRRs without a locus anchor are excluded. Unidentifiable lengths remain blank.

Read-supported non-repeat sequence is reconstructed after repeat inference;
uncovered bases become N. Products enter the same genotype and SNP-mask path as
assembly and LR products. See [formulas, diagnostics, CLI semantics and
limitations](../concepts/locus-reconstruction.md). For regression comparison,
`--sr-engine competitive` retains the [legacy workflow](competitive-fastq.md).

## Performance and CPU allocation

`--threads` now applies to short-read recruitment as well as later locus fitting
and reference mapping. Recruitment uses ordered 256-pair chunks and process
workers to avoid Python's interpreter lock. Each sample uses its existing CPU
allocation: directory/manifest batches first divide the global budget among
samples. Recruitment workers finish before the locus fitting pool starts. Small
inputs use fewer workers, and single-chunk inputs run without process startup.
Only two chunks per worker may be in flight, including completed results waiting
for earlier chunks; the full FASTQ is never submitted at once.

The motif test retains the same cyclic Hamming comparison and 90% threshold,
but runs its base/phase comparisons in NumPy. Bounded caches reuse motif phases,
read profiles and flank results. Synthetic candidate expansion reuses previously
calculated read/candidate alignment scores. Recruitment seed lengths and evidence
thresholds are unchanged. Ambiguous-molecule audit rows are streamed to their
existing TSV instead of accumulating in RAM. Uniquely recruited evidence remains
in memory for inference; output volume and fitting cost still scale with evidence.
No new dependencies or command-line switches are required.

Progress messages now separate recruitment, repeat fitting, likelihood output,
sequence reconstruction and genotype output. `reconstruction_metadata.json`
adds `performance.stage_seconds` and `performance.recruitment` (actual worker
count, chunk size, examined pairs and locus tests). These timings exclude the
preceding FASTQ QC and subsequent reference classification. Each process has
its own bounded caches, so memory still increases with worker count.

A local synthetic benchmark (12 loci, 150-base reads, 12/60-base motifs) measured:

| Measurement | Before | After |
| --- | ---: | ---: |
| Same serial classification workload, 500 pairs | 9.52 s | 1.42 s |
| Optimized recruitment, 3,000 pairs | 8.46 s with 1 worker | 2.55 s with 4 workers |

The first comparison was against commit `40f3d96`, without profiler overhead;
retained classification signatures were identical. Serial and parallel recruitment
also returned identical evidence, anchor metadata and ordering. These are local
recruitment benchmarks, not whole-sample or SLURM runtime guarantees. Database
classification and audit-file I/O can remain substantial costs.

Use [the recruitment benchmark](../../scripts/benchmark_short_read_recruitment.py)
on a representative subset in your allocation to measure actual scaling:

```bash
python scripts/benchmark_short_read_recruitment.py \
  --reads1 R1.fastq.gz --reads2 R2.fastq.gz \
  --primers panel.tsv --database reference_build \
  --pairs 10000 --threads 1 4 8 --output recruitment_timing.json
```

Choose worker counts within the allocated CPUs. The benchmark includes worker
startup and evidence hashing and verifies an identical evidence fingerprint for
all requested worker counts. It uses the default SR QC thresholds; it does not
run repeat inference or reference classification.
