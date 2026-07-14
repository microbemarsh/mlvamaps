# acceleration backends

MLVA Seer should use biological sequence tooling implemented in Rust or C where
possible, and avoid generic fuzzy-string packages for core matching.

Current backend policy:

- `amplirust` is the required in-silico PCR backend for assembly, contig, and
  reference product extraction and for paired-primer assignment of FASTQ reads
  after a lossless FASTA projection. It owns IUPAC primer interpretation,
  product-length constraints, strand handling, and primer alignment CIGARs.
- `sassy` is the preferred Rust/SIMD approximate DNA matcher for primer-style
  searches within already assigned amplicons, such as optional locus-flank
  localization. Degenerate primer pairing is delegated to `amplirust`.
- `vsearch>=2.30` performs exact dereplication and abundance-sorted global
  clustering independently for each locus. Its SIMD alignment and identity
  calculation include gaps. MLVA Seer consumes UC memberships and uses the
  observed cluster seed as the representative; it never uses a consensus.
- Parasail's SIMD C implementation provides exact end-to-end Needleman-Wunsch
  tracebacks used to annotate per-read indels against the selected observed
  repeat representative. Its native calls release the GIL, so independent
  unique sequences are aligned with threads rather than serialized processes.
  Identical repeat sequences are aligned once per representative, and independent
  unique sequences use the resolved worker allocation.
- NumPy performs quality-score reductions, batched repeat-motif comparisons,
  and per-read repeat-count likelihood vectors in compiled loops.
  `pysam.FastxFile` delegates FASTA/FASTQ parsing to htslib.
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
- Loci are submitted to VSEARCH sequentially, and each `cluster_size` process
  receives the full resolved thread count. This avoids competing native thread
  pools while keeping locus membership isolated.
- After VSEARCH clustering completes, Parasail alignments for independent unique
  repeat sequences use up to the resolved number of worker threads. This phase
  does not overlap VSEARCH.
