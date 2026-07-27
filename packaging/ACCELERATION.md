# acceleration backends

MLVAMaps should use biological sequence tooling implemented in Rust or C where
possible, and avoid generic fuzzy-string packages for core matching.

Current backend policy:

- `sassy-rs>=0.2.6` is the Rust/SIMD engine for in-silico PCR, paired-primer
  FASTQ assignment, and bounded flank localization. MLVAMaps owns deterministic
  IUPAC expansion, MLVA_finder-compatible strand fallback and product pairing,
  product-length constraints, and result selection.
- `parasail` computes exact Needleman-Wunsch global tracebacks between each
  mapped read repeat and its diagnostic product-group representative.
- `spoars` performs SIMD-accelerated partial-order assembly of complete reads
  in the dominant mapping-derived product group.
- NumPy performs quality-score reductions, batched repeat-motif comparisons,
  and per-read repeat-count likelihood vectors in compiled loops.
  `pysam.FastxFile` delegates FASTA/FASTQ parsing to htslib.
- `minimap2` performs competitive FASTQ locus/product recruitment and maps
  locus reads back to assembly-PCR-resolved SPOARS products; `pysam` parses
  SAM evidence.
- `minimap2` also maps accurate reads to extracted assembly products for depth
  support. `pysam` handles existing SAM/BAM support supplied by the user.
- MUMmer4 `dnadiff` performs exact whole-genome alignments only when an assembly
  query has tied exact marker matches. Independent tied-reference comparisons
  share the configured worker pool.

Default threading policy:

- CLI options use 32 threads by default. Users can explicitly pass `--threads 0`
  to use all available CPUs.
- `0` means auto-detect available CPUs for MLVAMaps workers.
- Native backends receive the resolved thread count directly. MLVAMaps does
  not place a Python thread pool around Sassy's internally threaded batch
  search, which avoids nested parallelism and CPU oversubscription.
- Native mapping and phylogenetic tools receive the resolved thread count.
