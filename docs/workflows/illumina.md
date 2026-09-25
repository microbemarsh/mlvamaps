# Illumina short-read workflow

The default `--sr-engine competitive` maps reads against candidate MLVA allele
contexts with minimap2 and applies shared allele inference. See the
[competitive workflow](competitive-fastq.md) for evidence and output details.

```bash
mlvamaps call --loci panel.tsv -i sr --fq1 reads_1.fastq.gz --fq2 reads_2.fastq.gz \
  --sample-mode metagenome -o results
```

## Optional repeat-likelihood engine

Selecting `--sr-engine repeat-likelihood` recruits molecules by non-repeat
flanks, then combines E (enclosing), S (spanning pair), F (flanking), and anchored
FRR evidence. Repeat states are synthetic and can be novel or outside the panel's
expected range. Haploid likelihood and the existing microbial EM model retain
distinct components, including supported SNP haplotypes with the same count.

```bash
mlvamaps call --loci panel.tsv -i sr --fq1 reads_1.fastq.gz --fq2 reads_2.fastq.gz \
  --sr-engine repeat-likelihood --insert-mean 400 --insert-sd 35 \
  --sample-mode metagenome -o results
```

Insert overrides must describe independently measured library statistics. They
are optional when repeat-independent flank pairs can estimate them or enclosing
reads determine the repeat directly. Single-end data can use E and F evidence;
FRRs without a locus anchor are excluded. Unidentifiable lengths remain blank.

Read-supported non-repeat sequence is reconstructed after repeat inference;
uncovered bases become N. Products enter the same genotype and SNP-mask path as
assembly and LR products. See [formulas, diagnostics, CLI semantics and
limitations](../concepts/locus-reconstruction.md).

## Performance and CPU allocation

The recruitment and repeat-fitting details below apply to
`--sr-engine repeat-likelihood`. Native FASTQ I/O options also apply to competitive runs.

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
thresholds are unchanged. Ambiguous-molecule diagnostics are streamed instead of
accumulating in RAM (see compact/full audit modes below). Uniquely recruited evidence remains
in memory for inference; output volume and fitting cost still scale with evidence.
The core optimizations require no additional dependencies or switches. Optional
native I/O libraries described below further accelerate compressed FASTQs.

Progress messages now separate recruitment, repeat fitting, likelihood output,
sequence reconstruction and genotype output. `reconstruction_metadata.json`
adds `performance.stage_seconds` and `performance.recruitment` (actual worker
count, chunk size, examined pairs, locus tests, skipped locus tests, unique,
ambiguous and unmatched pair counts, audit mode and pairs per second). Recruitment now includes
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
startup and evidence hashing and compares both full and compact audit modes by
default. It verifies identical retained evidence and insert calibration for all
requested worker counts and modes, and identical diagnostics across worker
counts within each mode. Use `--audit-modes compact` to time only this engine's
default audit mode. It uses the default SR QC thresholds; it does not
run repeat inference or reference classification.

### Resolving competing flank matches

Recruitment accounts for all candidate loci before assigning a pair. A short
incidental flank match no longer causes a much stronger locus match to be
discarded. Assignment uses the sum of the strongest flank alignment per mate;
the winner must exceed the runner-up by both 12 score points and 20% of its own
score. Close matches remain ambiguous. These guards are not calibrated
probabilities. The seed scan stops once all loci are already candidates.

For panels with more than three candidate loci, fast score-only alignments
bound the maximum possible anchor score on both strands. Candidates are checked
in descending bound order. Full anchor checks stop only when the remaining
bounds cannot change the winning assignment, or (in compact mode) the two
strongest ambiguous matches. Ignoring identity and repeat-only exclusions in
these bounds makes them conservative; those exclusions still apply to every
accepted anchor. Small candidate sets and panels with at most three distinct
flank contexts avoid the extra scoring pass.

Full repeat evidence is built only for the retained locus and requested full
ambiguity rows. Bounded caches reuse mate scores and quality-independent anchor
alignments; each retained molecule still supplies its own qualities, evidence
and insert contribution. `locus_tests` counts exact anchor checks and
`skipped_locus_tests` counts candidates ruled out by the score bounds.

In compact mode, `molecule_candidate_evidence.tsv` contains uniquely recruited
evidence. `ambiguous_molecules.tsv` records one row per excluded pair with its
sample ID, molecule ID, two witness locus IDs and `candidate_search_complete`
(`yes` when the candidate ranking is resolved, including by score bounds).
The witnesses are the two
strongest matching loci, not an exhaustive list of matches.
`--sr-recruitment-audit full` writes the
exhaustive per-locus `ambiguous_locus` rows and alignments in
`molecule_candidate_evidence.tsv`. This option applies to the repeat-likelihood
engine; it does not change the legacy competitive engine.

Progress now reports `Read pairs scanned`, unique and ambiguous counts, and
pairs/s. Compact mode reduces ambiguity output and avoids unnecessary full
evidence construction. It never rejects a pair merely because two loci match:
a later candidate may be the correct locus. Both audit modes use the same
assignment rule and insert calibration. Runtime depends on panel ambiguity;
use `scripts/benchmark_short_read_recruitment.py` on representative FASTQs to
measure throughput and verify matching evidence across worker counts.
An existing running process must be restarted with the updated installation to
use these changes.

Regression tests cover perfect 150-base paired reads across 25 loci with
incidental competing matches, near-identical loci that must remain ambiguous,
and agreement between compact/full audits and serial/parallel recruitment.

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
format full-mode ambiguous audit TSV chunks before sending them to the parent. These
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

### Long database-derived motifs during recruitment

When a panel lacks a concrete nucleotide motif, its database-derived motif can
be a whole observed repeat tract. Comparing every rotation of a long tract with
every background read is expensive, even before repeat fitting begins.

For primitive motifs of at least 128 bases, recruitment first searches for exact
read pieces using native string search. A read allowing `floor(length/10)`
substitutions is divided into one more piece than that error budget. Every
accepted cyclic match must contain at least one exact piece, so this filter
cannot remove a phase accepted by the existing 90% identity rule. Proposed phases
are verified with NumPy; low-complexity cases with many hits use the previous
dense comparison. No anchor cutoff, evidence class, or read count is changed.

A synthetic recruitment workload with 12 loci, 1,800-base database-style motifs,
and 500 paired 150-base reads decreased from 23.83 s to 3.08 s in a single worker.
All 6,000 locus tests and the complete evidence fingerprint were identical.
This measures recruitment with long motifs, not the runtime of a real sample.
The existing recruitment benchmark can measure your actual panel and inputs.

Recovery updates now go to **stdout**, alongside `Recovering loci using
repeat_likelihood`, and are flushed immediately. The log distinguishes template
loading, recruitment, fitting, likelihood output, and reconstruction even when
SLURM sends stderr elsewhere. `reconstruction_metadata.json` includes template
loading time and `performance.recruitment.template_motif_lengths` so long
database-derived motifs are visible in completed runs.

### Canonical genotype and evidence output

The output stage separates repeat/SNP masking, canonical sequence alignment,
molecule memberships, allele predictions, and locus summary tables. Each stage
reports progress to stdout with immediate flushing. The subsequent reference
classification message also flushes immediately, so buffered SLURM output does
not leave the earlier canonical-output message visible during classification.
`reconstruction_metadata.json` records `performance.output` timings for
genotyping, product tables, canonical alignments, and compatibility tables.
These are breakdowns of `stage_seconds.genotype_and_evidence_output`, not
additional time to sum with it.

The shared repeat masker checks each candidate unit once using bounded NumPy
blocks, then joins adjacent accepted units by phase. This replaces repeated
rescanning of long tracts while preserving the mismatch allowance, IUPAC
semantics, and earliest-start tie rule. Primer/flank searches use compiled
native regular expressions. Repeat-count-only tables use the existing length
calibration directly, without calculating unused SNP haplotypes. Identical
canonical SNP sequences use their exact ungapped alignment instead of running
dynamic programming; differing sequences still use Parasail.

Membership and prediction tables now stream in 8,192-row batches through the
native CSV writer. Output columns, row order, molecule IDs and evidence are
unchanged. On a local synthetic workload with two repeat variants and 100,000
molecule memberships, canonical output took **0.49 s before and 0.18 s after**.
Peak Python allocation during membership/prediction output fell from **45.7 MiB
to 0.77 MiB**, excluding the already-resident products and molecule IDs. All
output files matched after gzip decompression. A separate 4,000-base perfect
repeat masking example fell from 0.392 s to 0.00056 s with identical boundaries.
These timings measure these specific output/masking workloads; they do not
predict recruitment time or whole-sample runtime on shared cluster storage.
