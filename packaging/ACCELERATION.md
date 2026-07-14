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
- `savont>=0.6.1` is the required multithreaded ASV clustering, variant
  detection, consensus, and EM abundance backend for FASTQ calls. All locus
  FASTQs are submitted in one pooled invocation.
- WFA2, through `pywfa`, provides exact end-to-end Levenshtein tracebacks used
  to annotate per-read indels against final Savont repeat consensuses.
  Identical repeat sequences are aligned once per consensus, and independent
  unique sequences use the resolved worker allocation.
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
- Savont likewise receives the full resolved thread count in one process; MLVA
  Seer does not launch one competing Savont process per locus.
- After Savont completes, WFA2 alignments for independent unique repeat
  sequences use up to the resolved number of worker processes. This phase does
  not overlap Savont.
