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
The core optimizations require no additional dependencies or switches. Optional
native I/O libraries described below further accelerate compressed FASTQs.

Progress messages now separate recruitment, repeat fitting, likelihood output,
sequence reconstruction and genotype output. `reconstruction_metadata.json`
adds `performance.stage_seconds` and `performance.recruitment` (actual worker
count, chunk size, examined pairs and locus tests). Recruitment now includes
streamed input/QC time. Subsequent reference classification is timed separately
in `short_read_run_metadata.json`. Each process has
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

### Streaming QC and native gzip libraries

QC now feeds recruitment directly, so workers start before the complete sample
has been filtered. Filtered paired FASTQs are only materialized when needed for
database classification, `--keep-intermediates`, or the legacy competitive engine.
With `--database`, writes overlap recruitment and the recruitment stage no longer
decompresses those files again. Orphans are still processed after paired survivors
to retain the existing evidence order. QC thresholds and read support are unchanged.

Optional accelerators can be installed in the active environment:

```bash
conda install -c conda-forge rapidgzip python-isal
# Or install the Python extra from this checkout:
python -m pip install '.[fastq]'
```

[Rapidgzip](https://github.com/mxmlnkn/rapidgzip) decompresses ordinary gzip input
in native code using its Python API. The adapter feeds a bounded OS pipe into
the existing pysam/HTSlib parser; it keeps wrapped FASTQ support and pair validation.
[ISA-L](https://python-isal.readthedocs.io/en/stable/) accelerates gzip output
through `isal.igzip`. Small files and installations without rapidgzip retain
the direct HTSlib reader; without ISA-L, writes retain standard Python gzip.
No native extension is built by mlvamaps itself.

Automatic rapidgzip selection requires a regular `.gz` file of at least 8 MiB,
an installed Python module, `/dev/fd`, and enough allocated CPUs to leave at
least half for recruitment. For two eligible inputs and `--threads 32`, the budget
is 25 recruitment workers, two decoder threads per input, two pipe feeders, and
the QC parent. Rapidgzip may also create internal coordinating threads. Decoder
parallelism is always explicit; it never defaults to the entire SLURM node.
Later fitting and mapping can use the full allocation after input workers exit.

Conda-forge's [rapidgzip package](https://github.com/conda-forge/rapidgzip-feedstock)
currently lists Linux x86-64 and Intel macOS builds, plus Windows. Apple Silicon
can use its PyPI wheel via the `fastq` extra. This adapter uses POSIX `/dev/fd`;
other platforms keep the original reader. `python-isal` is also included in
`environment.yml` and supports the project's Linux/macOS platforms.

Flank alignment now uses Parasail's SIMD byte lanes with an explicit wider-lane
fallback on saturation. Impossible anchor scores skip traceback decoding, using
a bound derived from the existing scoring and edit cutoff. Worker processes
format ambiguous audit TSV chunks before sending them to the parent. These
changes preserve evidence rather than subsampling it or changing anchor cutoffs.

`short_read_run_metadata.json` records `performance.io`, the selected gzip writer,
QC/input/write time, recovery time, reference classification, report time, and
total time. QC time is nested within streamed recovery, so those entries must
not be added together. The competitive engine still performs QC before recovery.

On a local synthetic dataset of 300,000 paired 150-base reads (102 MB compressed),
input/QC/materialization plus feeding the next stage took **7.76 s at `4b90032`
versus 4.07 s with streaming and both accelerators**. Both processed the same
600,000 reads with identical QC counts. This approximately 1.9× result is for
that part of the workflow, not complete-sample or cluster throughput. Native
libraries alone provided a smaller gain; eliminating the redundant pass and
per-read NumPy reduction accounted for additional savings.

Measure I/O independently on your allocation and storage:

```bash
python scripts/benchmark_fastq_io.py \
  --reads1 R1.fastq.gz --reads2 R2.fastq.gz \
  --pairs 300000 --threads 8 \
  --scratch-dir "$SLURM_TMPDIR" --output fastq_io_timing.json
```

The benchmark compares HTSlib/standard gzip, HTSlib/ISA-L, and rapidgzip/ISA-L,
including current QC, output writes, and retained-read hashing. It checks equal
QC counts and retained-read hashes. Temporary outputs are removed automatically.
Use a sufficiently large subset: decoder startup can outweigh gains on tiny inputs.

### Repeat fitting and reconstruction

The recovery stage also reduces native alignment work. Repeat likelihoods use
Parasail's 16-bit SIMD scores, with a 32-bit fallback if the narrower calculation
saturates. An exact read substring already attains the maximum possible score
and requires no dynamic-programming alignment.

For long repeat candidates, the scorer reuses mathematically equivalent scores.
Under the existing match reward (2) and gap extension cost (1), a positive local
alignment of an L-base read spans fewer than 3L reference bases. A temporary
scoring target can therefore remove whole motif periods beyond 3L plus one
period, preserving the right-boundary phase and all possible scoring windows.
The full candidate grid, length likelihoods, output coordinates, repeat counts,
and molecule weights are retained. This shortening is never used for traceback
coordinates or reconstructed locus sequences. Range expansion reuses already
computed scores, including equivalent long targets.

In an unprofiled synthetic long-repeat case (500 flanking molecules, 60-base
motif, 201 candidate states), repeat fitting decreased from **38.30 s to 0.42 s**.
All 100,500 molecule/candidate scores, joint scores and posterior values were
exactly equal. This roughly 90× improvement applies to that fitting workload;
it is not a whole-sample runtime prediction. Short motifs, recruitment-heavy
samples, and reference classification have different cost profiles.

Reconstruction reuses raw alignments across identical reads and flank projections
across identical sequence/quality pairs, with bounded caches. Every original
molecule still contributes its own vote and support count. Full-product
haplotypes are grouped in one pass instead of rescanning all molecules for each
variant. Likelihood TSVs use the native CSV row writer directly, preserving all
records and the existing columns.

Progress reports individual locus fits and likelihood-output row counts.
`reconstruction_metadata.json` also records `performance.loci`, including each
locus's fitting time, reconstruction time, evidence count and candidate count.

Reproduce the exhaustive-versus-optimized fitting comparison:

```bash
python scripts/benchmark_repeat_fitting.py \
  --molecules 500 --unit 60 --maximum 100 --output repeat_fitting_timing.json
```

The benchmark verifies exact equality of every inference field, including
mixture fractions and confidence, and excludes QC, recruitment and file writing.
