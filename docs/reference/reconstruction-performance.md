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
