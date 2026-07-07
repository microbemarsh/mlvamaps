# acceleration backends

MLVA Seer should use biological sequence tooling implemented in Rust or C where
possible, and avoid generic fuzzy-string packages for core matching.

Current backend policy:

- `amplirust` is the required in-silico PCR backend for assembly, contig, and
  reference product extraction.
- `sassy` is the preferred Rust/SIMD approximate DNA matcher for primer-style
  searches. The assembly extraction path already gets sassy-backed matching via
  `amplirust`; the direct FASTQ read-assignment path prefers the `sassy-rs`
  Python binding when it is installed.
- `edlib` is only a C-backed fallback for approximate primer/flank alignment
  when `sassy-rs` is not available.
- `minimap2` plus `pysam` are the preferred low-level alignment stack for
  future read-to-reference-amplicon and assembly evidence.

Default threading policy:

- CLI options use `--threads 0` by default.
- `0` means auto-detect available CPUs for MLVA Seer workers or delegate
  auto-detection to the external backend, such as `amplirust`.
