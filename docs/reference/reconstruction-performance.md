# Short-read reconstruction performance

The high-depth pileup implementation at baseline `5bdcbc6` spent most of its
reconstruction time in Python overlap comparisons. Preparation also assembled
some loci a second time during final calling, and CPU-bound locus recovery used
threads. The revised path uses native NumPy mismatch counting and molecule
voting, reuses unchanged preparation results, and runs preparation and recovery
in spawned process workers.

One pool per sample is reused across these locus stages. Its size is bounded by
the sample's `--threads` allocation and number of loci with recruited evidence.
Tasks and results remain ordered with bounded submissions. Recruitment and
reconstruction computation run in successive stages; the reconstruction pool
is idle during additional read recruitment. Single-locus and empty inputs avoid
starting this pool. Worker-side molecule classifications are copied back into
the audit, preserving serial/process equivalence.

The native vote matrix is bounded to 8 MiB per worker. Temporary arrays,
retained reads and process runtimes consume additional memory; this is not a
total worker memory limit. Conflicting mates still contribute one unknown-base
vote per molecule, and uncovered positions receive no votes.

## Local measurements

Identical generated FASTQs: 12,000 pairs, six primer-only loci, varied read
starts, substitution errors and 25% off-target background. Input generation and
initial parent imports are excluded; calling and output writing are included.
Each configuration was measured once on the same host, sequentially. These are
synthetic measurements, not a prediction for server libraries.

| Configuration | CPU budget | Wall time |
| --- | ---: | ---: |
| Baseline `5bdcbc6` | 4 | 16.09 s |
| Native operations and process workers | 4 | 7.09 s |
| Native operations, serial | 1 | 13.15 s |

The same-budget run was 2.27 times faster. All three runs returned the same six
repeat counts and PASS statuses. Process concurrency reduced runtime relative
to the new serial path as well. Parent peak RSS was 146.3 MiB before and
117.4 MiB after; maximum child peak RSS was 67.5 and 92.9 MiB respectively.
These are separate process measurements, not aggregate concurrent memory.
Multiprocessing adds worker memory alongside its speed benefit.

[Raw measurements and environment](reconstruction-performance.json) retain
timings, allocations, counts and per-process resource measurements.

## Reproduce

```bash
python scripts/benchmark_reconstruction.py --output bench/pileup-1 --molecules 12000 --scenario pileup --threads 1
python scripts/benchmark_reconstruction.py --output bench/pileup-4 --molecules 12000 --scenario pileup --threads 4
python -m pytest -q
```

To compare an older checkout, run its caller on the same generated `panel.tsv`,
`reads1.fq` and `reads2.fq` with a separate output directory. The new benchmark
includes reconstruction stage timings and parent/child peak RSS. Regression
checks cover reuse of prepared products, process/serial call and audit equality,
native vote chunk boundaries, and the high-depth FASTQ/assembly count matches.

## Competitive caller audit against `b7f6a60`

The competitive caller already used spawned processes for read recruitment and
locus recovery, Parasail SIMD for alignment, HTSlib for parsing/pileups, and
NumPy for likelihoods. Three expensive stages were still serial: sample-anchor
learning, depth graph construction, and whole-input depth counting. These now
reuse the existing locus process pool with bounded submissions. Depth counting
encodes exact 15/21-mers as NumPy uint64 values and counts only graph targets;
quality failures, ambiguous bases and read boundaries break k-mers. No minimum
depth or molecule-support requirement was added.

A serial CPU profile also found overlap searches dominating the depth fixture:
1.33 million `_overlaps` calls and 26.2 million native substring searches.
Overlap verification now searches only enough disjoint seed blocks to exceed
the maximum allowed substitution count. This preserves accepted overlaps while
reducing repeated work in periodic tails. Inferred flank projection also reuses
its existing bounded alignment cache across molecules, with quality and anchor
checks still evaluated for each molecule.

The same deterministic 12,000-pair, six-locus `depth` fixture was run before and
after with four CPUs. These single-run synthetic timings exclude input
generation and initial imports, and include calling and output writing.

| Stage | `b7f6a60` | Revised |
| --- | ---: | ---: |
| Entire call | 14.85 s | 10.35 s |
| Sample learning and recruitment | 4.24 s | 2.94 s |
| Depth graphs and counting | 1.64 s | 0.57 s |
| Locus recovery, including depth and PCR | 8.60 s | 5.58 s |

Whole-call runtime decreased by 30%. The depth row is part of locus recovery;
these rows must not be summed. Biological TSV/FASTA outputs, learned graphs and
depth estimates matched; the summary differed only in the input file paths.
The provisional estimates remain 69, 24.5, 5, 12, 14 and 9.5, against true counts
70, 26, 5, 11, 14 and 9. This benchmark measures performance equivalence, not
exact recovery or superiority over Shovill. Raw measurements are under
`competitive_audit` in [the measurement JSON](reconstruction-performance.json).

### Remaining limits and the 32-CPU SLURM invocation

One sample discovered by `call -i ... --short-reads --threads 32` receives the
full 32-CPU budget. Locus workers are capped by the number of populated loci;
one difficult locus cannot use every CPU during its graph reconstruction.
The overlap graph still has a quadratic worst case within its node limit.
Sample-anchor index construction remains serial and bounded to 256 sequences
per locus. Small reconstructed-contig PCR and output preparation also remain
serial; PCR took under 0.01 seconds in this depth fixture.

FASTQ parsing uses HTSlib, but the parent feeds/QCs reads. Parallel gzip decoding
requires the optional `fastq` dependencies and sufficiently large gzip inputs.
With rapidgzip available, two large gzip mates and 32 CPUs reserve seven CPUs
for I/O and feeding, leaving 25 recruitment workers. Without it, the input log
reports `htslib` and decoding is serial. Gzip intermediate writing uses ISA-L
when installed, otherwise Python's gzip wrapper around native zlib.

`--database` additionally writes filtered FASTQs and runs downstream reference
classification. Minimap2 receives the thread budget; BAM-to-evidence conversion,
per-molecule reference bookkeeping and some matrix preparation remain serial.
The numerical classification calculation uses NumPy. Database size and the
number of competing alignments can therefore add substantial work not measured
by the database-free fixture. Inspect `reference_classification` in
`short_read_run_metadata.json` separately from reconstruction stage timings in
`reconstruction_metadata.json`.

Symlinking `/blue` FASTQs into `$SLURM_TMPDIR` still reads shared storage on every
scan. To measure local-disk performance, copy the two inputs into task-local
scratch before calling. Running into task-local output and copying the completed
sample directory back also avoids shared-storage intermediate writes. Only do
this when the allocation provides enough scratch space. This repository audit
does not measure the user's server filesystem or reference database.

### Further reuse without changing calling decisions

A second pass reduced the same four-CPU depth benchmark from 10.35 to 8.24
seconds: another 20% reduction, or 44% below `b7f6a60`. Overlap assembly now
computes repeat-only prefix/suffix masks once per node with NumPy cumulative
counts, and reuses exact seed positions across node pairs. Both use bounded
caches. The seed-position cache holds only eight nodes because graph traversal
processes one left node at a time. The 90% repeat exclusion rule, mismatch
tolerance, accepted offsets and retained molecules are unchanged.

Classification reuses alignment statistics only when the edit strings, mapped
coordinates, repeat coordinates, unit length and repeat-count availability
match. Each molecule still contributes independently. Equivalent reference
columns use NumPy byte signatures instead of tuples of Python floats; signed
zero is normalized to preserve the previous equality rule.

On a separate synthetic fixture of 80,000 alignments, 1,000 molecules and 80
competing references, likelihood construction and reference grouping decreased
from 0.430 to 0.234 seconds. All likelihoods, error estimates, groups and
diagnostics matched exactly. This excludes minimap2, BAM decoding and filesystem
I/O; it is not a whole-classification speedup estimate for the server database.
The fixture cycles perfect, substitution and repeat-deletion alignments to
exercise repeated alignment patterns.

The depth benchmark's biological outputs also matched the preceding pass.
Regression checks compare every prefix/suffix length against scalar repeat
testing, overlap offsets against exhaustive comparison, and classification
reuse across different repeat contexts. The full suite passed 473 tests, with
eight optional-backend skips. Raw results are in `reuse_audit` in the measurement
JSON. These checks demonstrate preserved behavior on the tested inputs, not
universal biological accuracy.
