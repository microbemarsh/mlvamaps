# bioconda packaging notes

This directory is a staging area for the future bioconda recipe. The actual
recipe should be submitted to the external `bioconda/bioconda-recipes`
repository after MLVAMaps has a tagged GitHub release.

Current packaging assumptions:

- MLVAMaps is pure Python, so the recipe should use `noarch: python`.
- `mlvamaps` is the console script.
- `amplirust` is a runtime dependency because assembly extraction and FASTQ
  paired-primer assignment delegate IUPAC matching to its executable.
- `minimap2` should be a runtime dependency because assembly calls can use
  `--reads` to add read-depth support.
- `minibwa` should be a runtime dependency because FASTQ calls map all
  primer-oriented locus reads to the dominant observed representative amplicon
  and report reference-relative SNP evidence by default.
- `sassy-rs` must be packaged separately for Bioconda because bounded flank
  localization within assigned amplicons uses its Python binding.
- `parasail` supplies SIMD C global-alignment tracebacks for per-read
  substitutions and indels against each observed cluster representative.
- `pysam`/htslib handles FASTA, FASTQ, SAM, and BAM parsing, while NumPy moves
  quality, motif, and repeat-likelihood reductions out of Python loops.
- Release source should use a stable GitHub archive URL for a tag, with a
  `sha256` checksum filled in before submission.
- bioconda tests should avoid requiring local data. Import checks and CLI help
  checks are enough for package installation validation.

Before submitting to bioconda:

1. Create a GitHub release tag, for example `v0.1.0`.
2. Download the release tarball and calculate its `sha256`.
3. Copy `meta.yaml` into `bioconda-recipes/recipes/mlvamaps/meta.yaml`.
4. Replace the placeholder checksum.
5. Run the local bioconda recipe tests from the `bioconda-recipes` checkout.

## Release checklist moved from the main README

The recipe in this directory is not meant to be consumed directly from this
repository. bioconda recipes are submitted to the external
`bioconda/bioconda-recipes` repository after a tagged release is available.

Before submission we will need to:

1. Publish a stable GitHub release tag.
2. Replace the recipe URL/checksum placeholders with the release tarball and
   `sha256`.
3. Confirm the package name is not already present in upstream conda channels.
4. Run the bioconda recipe tests.
