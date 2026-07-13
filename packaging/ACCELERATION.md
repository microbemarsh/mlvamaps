# acceleration backends

MLVA Seer should use biological sequence tooling implemented in Rust or C where
possible, and avoid generic fuzzy-string packages for core matching.

Current backend policy:

- `amplirust` is the required in-silico PCR backend for assembly, contig, and
  reference product extraction.
- `sassy` is the preferred Rust/SIMD approximate DNA matcher for primer-style
  searches. The assembly extraction path already gets sassy-backed matching via
  `amplirust`; the direct FASTQ read-assignment path uses the `sassy-rs>=0.2.4`
  batched API and gives Sassy ownership of all requested worker threads.
- `edlib` is only a C-backed fallback for approximate primer/flank alignment
  when `sassy-rs` is not available.
- `minimap2` plus `pysam` are the preferred low-level alignment stack for
  future read-to-reference-amplicon and assembly evidence.

Default threading policy:

- CLI options use 32 threads by default. Users can explicitly pass `--threads 0`
  to use all available CPUs.
- `0` means auto-detect available CPUs for MLVA Seer workers or delegate
  auto-detection to the external backend, such as `amplirust`.
- Native backends receive the resolved thread count directly. MLVA Seer does
  not place a Python thread pool around Sassy's internally threaded batch
  search, which avoids nested parallelism and CPU oversubscription.
