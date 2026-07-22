# acceleration backends

MLVAMaps should use biological sequence tooling implemented in Rust or C where
possible, and avoid generic fuzzy-string packages for core matching.

Current backend policy:

- `sassy-rs>=0.2.6` is the Rust/SIMD engine for in-silico PCR, paired-primer
  FASTQ assignment, and bounded flank localization. MLVAMaps owns deterministic
  IUPAC expansion, MLVA_finder-compatible strand fallback and product pairing,
  product-length constraints, and result selection.
- `vsearch>=2.30` performs exact dereplication and abundance-sorted global
  clustering independently for each locus. Its SIMD alignment and identity
  calculation include gaps. MLVAMaps consumes UC memberships and uses the
  observed cluster seed as the representative; it never uses a consensus.
- `parasail` computes exact Needleman-Wunsch global tracebacks between each
  unique repeat sequence and its observed VSEARCH centroid. Independent
  alignments are distributed across the configured Python thread pool.
- NumPy performs quality-score reductions, batched repeat-motif comparisons,
  and per-read repeat-count likelihood vectors in compiled loops.
  `pysam.FastxFile` delegates FASTA/FASTQ parsing to htslib.
- `minimap2` maps FASTQ locus reads back to dominant observed VSEARCH
  representative amplicons; `pysam` parses the resulting SAM for base depth
  and reference-relative SNP evidence.
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
- Loci are submitted to VSEARCH sequentially, and each `cluster_size` process
  receives the full resolved thread count. This avoids competing native thread
  pools while keeping locus membership isolated.
- After VSEARCH clustering, unique Parasail alignments share the resolved
  thread pool; identical repeat sequences reuse the same traceback metrics.
