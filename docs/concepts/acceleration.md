# Acceleration backends and threading

mlvamaps should use biological sequence tooling implemented in Rust or C where
possible, and avoid generic fuzzy-string packages for core matching.

Current backend policy:

- `minimap2` performs competitive alignment of both Illumina and long reads
  against versioned candidate allele contexts; `pysam` parses all retained
  primary and secondary alignments.
- The Bioconda `sassy>=0.2.2` command-line tool is the Rust/SIMD search engine
  used by in silico PCR, paired-primer FASTQ assignment, and bounded flank
  localization.
  mlvamaps owns deterministic IUPAC expansion, MLVA_finder-compatible strand
  fallback and product pairing, product-length constraints, and result selection.
- `parasail` computes exact Needleman-Wunsch global tracebacks between each
  mapped read repeat and its diagnostic product-group representative.
- `spoars` performs SIMD-accelerated partial-order assembly of complete reads
  in the dominant mapping-derived product group.
- NumPy performs quality-score reductions, batched repeat-motif comparisons,
  and per-read repeat-count likelihood vectors in compiled loops.
  `pysam.FastxFile` delegates FASTA/FASTQ parsing to htslib.
- SPOARS corrects and represents sequence variants after direct molecule
  evidence has already entered shared allele inference.
- `minimap2` also maps accurate reads to extracted assembly products for depth
  support. `pysam` handles existing SAM/BAM support supplied by the user.
- MUMmer4 `dnadiff` performs exact whole-genome alignments only when an assembly
  query has tied exact marker matches. Independent tied-reference comparisons
  share the configured worker pool.

Default threading policy:

- CLI options use 32 threads by default. Users can explicitly pass `--threads 0`
  to use all available CPUs.
- `0` means auto-detect available CPUs for mlvamaps workers.
- Native backends receive an appropriate share of the resolved thread count.
  The Sassy CLI adapter uses temporary FASTA inputs and runs each search with
  one Sassy thread; assembly records can be distributed across mlvamaps worker
  processes.
- Native mapping and phylogenetic tools receive the resolved thread count.
