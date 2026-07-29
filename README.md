# mlvamaps

mlvamaps is a panel-driven toolkit for analyzing microbial MLVA/VNTR typing
data from accurate amplicon or long reads containing the target products,
assemblies, and assemblies with supporting accurate reads.

It is designed to be organism-agnostic: the same pipeline can be used for any
microbe with a defined MLVA primer panel and enough repeat metadata to interpret
the products. mlvamaps does not hard-code one species, typing scheme, or profile
database.

mlvamaps returns conventional repeat-copy-number fingerprints while preserving
the sequence evidence behind them:

- Primer-supported locus detection in compressed or uncompressed FASTQ.
- Built-in, SIMD-accelerated in-silico PCR in either orientation, with matching
  behavior designed for MLVA_finder compatibility.
- Competitive mapping of long reads to locus and repeat-product classes.
- Per-read substitution and indel evidence within mapped product groups.
- SPOARS partial-order assembly of dominant mapped locus products.
- Emu-inspired expectation-maximization estimates of meaningful variant
  fractions from mapped-read counts.
- minimap2 mapping to sample-derived POA amplicons and
  reference-relative SNP evidence.
- Optional read-depth support for assembly products.
- Dedicated individual-locus repeat-count tables and report graphics, separate
  from the generated gel and SNP evidence.
- Optional per-locus phylogenetic placement against a sequence database using
  MAFFT, RAxML-NG, and EPA-ng.
- Probabilistic locus calls, profile matching, and a
  self-contained HTML report.

## Install

The recommended installation uses conda through
[Miniforge](https://github.com/conda-forge/miniforge):

```bash
git clone https://github.com/microbemarsh/mlvamaps.git
cd mlvamaps
conda env create -f environment.yml
conda activate mlvamaps
python setup.py install
```

The environment includes
[Sassy's Rust/Python bindings](https://github.com/RagnarGrootKoerkamp/sassy),
[Python regex](https://github.com/mrabarnett/mrab-regex),
[Parasail's Python bindings](https://github.com/jeffdaily/parasail-python),
[SPOARS' Python bindings](https://github.com/fg-labs/spoars),
[Deacon](https://github.com/bede/deacon),
[minimap2](https://github.com/lh3/minimap2),
[MAFFT](https://mafft.cbrc.jp/alignment/software/),
[RAxML-NG](https://github.com/amkozlov/raxml-ng),
[EPA-ng](https://github.com/pierrebarbera/epa-ng),
[MUMmer4](https://github.com/mummer4/mummer), and the Python/native
sequence libraries declared in `environment.yml`.

## Quick start

Analyze amplicon or other primer-spanning reads:

```bash
mlvamaps call primers.tsv sample.fastq.gz
```

Analyze an assembly:

```bash
mlvamaps call primers.tsv assembly.fasta
```

Analyze every supported FASTA or FASTQ file in a directory:

```bash
mlvamaps call primers.tsv sequence_files/ -o results
```

Directory input is non-recursive and may contain a mixture of FASTA and FASTQ
files. Each file is treated as one sample and written to
`results/<filename-stem>/`. Unrelated files are ignored.

Analyze an assembly with FASTQ read-depth support:

```bash
mlvamaps call primers.tsv assembly.fasta --reads sample.fastq.gz
```

Use an existing assembly-aligned SAM/BAM for depth support:

```bash
mlvamaps call primers.tsv assembly.fasta --bam assembly_reads.bam
```

Compare the fingerprint with known profiles:

```bash
mlvamaps call primers.tsv sample.fastq.gz --profiles profiles.tsv
```

Pre-screen a metagenomic FASTQ against a target-taxon Deacon pangenome index:

```bash
mlvamaps call primers.tsv metagenome.fastq.gz \
  --taxon-screen-index target_taxon.idx \
  --profiles profiles.tsv
```

Supply an existing target index to `--taxon-screen-index`. See
[bede/deacon-indexes](https://github.com/bede/deacon-indexes) for information
on building indexes.

The screen runs before MLVAMaps loads or quality-filters reads. Deacon receives
the shared `--threads` CPU budget and writes the retained reads plus its JSON
summary under `results/taxon_screen/`. The default Deacon retention thresholds
are two shared minimizers and a 1% relative match.

Place each callable query locus in fixed reference trees:

```bash
mlvamaps call primers.tsv sample.fastq.gz \
  --database reference_build \
  --profiles profiles.tsv
```

The recommended database is the top-level output from `mlvamaps
build-reference`. MLVAMaps reuses its fixed reference alignments, trees, and
selected models, then runs only MAFFT query addition and EPA-ng placement.
Independent locus placements run concurrently, and the available CPU budget is
shared between them. This is more efficient than assigning all OpenMP threads
to each short single-query placement. The
`reference_build/database` subdirectory is also accepted.

For a sequence-only database, use a directory containing one FASTA per locus.
Filename stems must match panel locus IDs, and FASTA headers must use the same
reference IDs across loci:

```text
reference_sequences/VNTR_01.fasta
reference_sequences/VNTR_02.fasta
```

Long-form sequence TSV and combined FASTA databases are also supported; see
[input formats](docs/reference/input-formats.md#phylogenetic-sequence-database).

Results are written to `results/` by default. Start with:

- `calls.tsv` for compact per-locus calls.
- `locus_repeat_counts.tsv` for explicit individual-locus repeat counts.
- `mlva_fingerprint.tsv` for the conventional wide fingerprint.
- `profile_matches.tsv` for ranked, metadata-rich profile comparisons.
- `profile_match_loci.tsv` for one machine-readable row per profile and locus.
- `report.html` for the visual summary.

With `--database`, `phylogeny/phylogenetic_matches.tsv` ranks complete
references by summed distance. Every database locus retains its MAFFT alignment
and fixed RAxML-NG reference tree/model; callable query loci also retain an EPA-ng
`.jplace` result. Rich panels additionally mask the tandem-repeat tract from
the SNP tree, preserve repeat count and repeat-unit haplotype separately, and
rank the combined normalized evidence in
`phylogeny/combined_marker_matches.tsv`. For assembly queries, exact marker
ties are resolved by canonical whole-genome identity followed by MUMmer4
`dnadiff` SNPs, indel bases, and one-to-one aligned fraction. Optional dated/geocoded reference
metadata can be joined to `phylogeny/combined_markers.tree` in MYOGA.

FASTQ runs additionally provide native primer-pair evidence under
`in_silico_pcr/`, mapping-derived variant groups, read memberships,
EM-estimated mixture abundance, minimap2 mapping coverage, and SNP evidence.
Assembly runs provide the same native primer-match evidence, extracted
products, and optional read support. No Amplirust executable is required or
invoked.

## Supported data

| Input | What mlvamaps assesses |
| --- | --- |
| Directory of FASTA/FASTQ files | Runs each supported top-level file as a separate sample under its own output subdirectory. |
| High-accuracy long-read FASTQ/FASTQ.GZ | Competitively recruited full and partial locus reads, presence evidence, local products, assembly-equivalent repeat counts, variants, and SNP evidence. |
| Accurate long-read WGS/metagenomic reads | Complete products are genotyped directly; repeat-spanning partial reads can provide provisional alleles and locus-specific partial reads establish untyped presence. |
| Assembly FASTA | In-silico primer products, product coordinates, sizes, and repeat counts. |
| Assembly plus accurate FASTQ | Assembly calls plus minimap2 read count and mean coverage for extracted products. |
| Assembly plus SAM/BAM | Assembly calls plus overlap-based read support from existing alignments. |
| Known profile TSV | Closest MLVA profiles, mismatched loci, distance, and comparison confidence. |
| Per-locus sequence database | Fixed-tree phylogenetic placement and a closest-reference ranking across callable loci. |

FASTQ mode competitively recruits reads to complete locus products before
primer pairing. Database products are preferred; rich panels can synthesize
auditable fallback templates. Presence-only mappings are never promoted to
repeat calls, while reads spanning both repeat boundaries can provide
provisional genotypes even when the complete product or whole genome cannot be
assembled. The mapping paths are scoped to accurate reads; noisy long-read
mapping is not supported.

FASTQ calls default to mean Q17 or better (approximately 98% per-base
accuracy) and retain singleton locus evidence. For every locus with complete
products, SPOARS builds a partial-order-alignment consensus from the dominant
read cluster. The normal assembly in-silico PCR and legacy product caller are
then run on that local contig. Individual reads determine confidence and
variant evidence; they do not redefine the assembled primary repeat count.
Low-depth calls remain in the fingerprint with an explicit `LOW_DEPTH` status.
Metagenome interpretation is the default and conservatively flags meaningful
secondary alleles; use `--sample-mode isolate` for cultured material.

Primary allele confidence increases when multiple reads in the dominant
sequence cluster agree. Secondary variants are evaluated separately:
single-read candidates remain visible for rapid detection, confirmed
secondaries trigger mixture interpretation, and neither is averaged into the
primary signature.

The FASTQ `report.html` includes a prominent **FASTQ Local Assembly
Concordance** table. It shows the raw read-product length range and mode,
SPOARS consensus length, assembly-PCR product length, raw/final repeat counts,
support, measurement source, and fallback status for every locus.

## MLVA_finder-compatible in-silico PCR

mlvamaps includes its own paired-primer engine for assembly extraction and
FASTQ locus assignment. Sassy performs SIMD-accelerated approximate matching
through its Rust Python binding. A small compatibility layer then resolves
fuzzy-alignment ties with the historical Python `regex` behavior used by
[i2bc/MLVA_finder](https://github.com/i2bc/MLVA_finder). This keeps the expensive
sequence scan in native code while retaining legacy match selection.

Compatibility behavior includes:

- deterministic expansion of IUPAC-degenerate primer bases;
- treating `N` in an input assembly as an error, rather than as a wildcard;
- successive per-primer error rounds from zero through
  `--max-primer-mismatches`;
- forward-strand-first matching with reverse-complement fallback at each error
  round;
- preference for equal-length fuzzy matches when both equal-length and indel
  matches are available at the same threshold;
- the MLVA_finder product-size formula, which uses configured primer lengths
  even when an observed primer match contains an insertion or deletion; and
- legacy assembly result selection rules, including FASTA record order and the
  smallest eligible unrounded allele on the final matching record.

The engine writes normalized primers, extracted products, coordinates, edit
costs, identities, CIGAR strings, strand, and product sequence to
`in_silico_pcr/`. These native files replace the former `amplirust/` evidence
directory.

## How it works

For FASTQ data, mlvamaps:

1. Filters reads by length and quality.
2. Uses the built-in Sassy-backed engine to pair degenerate primers and orient
   each MLVA_finder-compatible product.
3. Locates the repeat region and measures repeat/motif evidence.
4. Groups reads by their competitive locus/product mapping.
5. Uses mapped-read counts to distinguish dominant and secondary repeat-product
   groups.
6. Builds a dominant per-locus SPOARS POA contig and sends it through the same
   in-silico PCR, product-size calculation, and repeat caller as an assembly.
7. Maps locus reads back to that POA product for support and SNP evidence.
8. Uses supporting reads to determine confidence without redefining the
   assembly-derived allele.
10. Builds the fingerprint, compares profiles, and writes a
    plot-first HTML report.
11. When `--database` is supplied, separates the tandem-repeat tract from the
    SNP-bearing sequence, aligns repeat-masked references with MAFFT, infers a
    maximum-likelihood tree with RAxML-NG, and places the masked query with
    EPA-ng. It then combines normalized SNP-tree distance with the separately
    retained repeat-count distance. Only references present at every placed
    locus are ranked, so missing loci cannot produce an artificially small
    total.

For assemblies, mlvamaps:

1. Finds paired-primer products with the built-in Sassy-backed,
   MLVA_finder-compatible engine.
2. Selects valid products and converts size into repeat count where the panel
   provides enough metadata.
3. Optionally adds minimap2 or existing SAM/BAM read support.
4. Builds the same fingerprint, profile comparison, and report formats.
5. Optionally performs the same per-locus MAFFT, RAxML-NG, and EPA-ng fixed-tree
   placement from extracted assembly products.

FASTQ and assembly profile matches use repeat-count distance and matched-locus
count in the same order. FASTQ allele probabilities are used only to break an
otherwise equal profile match, so uncertainty cannot displace a closer
assembly-equivalent signature.

The minimap2 mapping coordinates are positions within the sample-derived
representative amplicon, not chromosome coordinates. The SNP table is
transparent within-sample evidence rather than a whole-genome or clinical VCF.

## Bring any microbial MLVA scheme

A minimal panel needs:

```text
locus_id
forward_primer
reverse_primer
```

Repeat-unit length, nominal repeat count, expected product size, repeat motif,
flanks, and valid size/count ranges make the resulting calls more informative.
No species name is required by the software.

See [adapting a panel for a new organism](docs/guides/new-organism-panel.md) for
recommended metadata, validation, and profile-table setup.

## Documentation

- [Documentation index](docs/README.md)
- [FASTQ and amplicon workflow](docs/workflows/fastq.md)
- [Assembly workflow](docs/workflows/assembly.md)
- [Input and panel formats](docs/reference/input-formats.md)
- [Output file reference](docs/reference/outputs.md)
- [CLI and thresholds](docs/reference/cli.md)
- [Representative mapping and SNP evidence](docs/concepts/representative-mapping.md)
- [Variant mixture abundance](docs/concepts/variant-mixtures.md)
- [Repeat-aware SNP placement and phylogeography](docs/concepts/repeat-snp-phylogeography.md)
- [Allele calling and profiles](docs/concepts/calling-and-profiles.md)
- [Adding a new organism or MLVA scheme](docs/guides/new-organism-panel.md)

## Additional commands

Build a reference sequence database and one maximum-likelihood phylogeny per
locus from a directory of assemblies. Assembly basenames must match the
metadata identifier unless the metadata has an `assembly_file`, `filename`, or
`path` column:

```bash
mlvamaps build-reference \
  --assemblies reference_assemblies/ \
  --primers primers.csv \
  --metadata metadata.csv \
  -o reference_build \
  -t 16
```

The metadata identifier may be named `reference_id`, `genome_id`, `sample_id`,
`strain`, `accession`, or `id`. The command writes:

- `reference_build/database/LOCUS.fasta`: raw amplicons, one FASTA record per
  reference, suitable for `mlvamaps call --database`;
- `reference_build/database/reference_sequence_index.tsv`: canonical sequence
  hashes for the default exact-reference fast path;
- `reference_build/database/reference_assemblies.tsv`: stable reference-ID,
  whole-genome assembly path, and canonical assembly SHA-256;
- `reference_build/phylogeny/LOCUS.tree`: a portable Newick tree for each locus;
- `reference_build/reference_build_manifest.tsv`: extraction and ambiguity QC;
- `reference_build/database/reference_metadata.tsv`: normalized placement metadata;
- `reference_build/myoga_metadata.csv`: metadata whose `genome_id` matches tree tips.

Multiple products at the same locus are excluded by default because an
unresolved paralog is unsafe as a phylogenetic reference. Review the manifest,
or use `--multiple-products best` only when choosing the best primer match is
appropriate. A rich `--loci` TSV containing the repeat motif or bounding flanks
is preferable to a primer-only CSV: it lets the tree builder mask the tandem
repeat for the SNP tree while retaining the unmasked amplicon in the database.

Simulate amplicon reads for pipeline testing:

```bash
mlvamaps simulate \
  --loci examples/mlva_loci.example.tsv \
  --sample-id SIM1 \
  --depth 500 \
  -o simulated
```

Extract MLVA_finder-compatible primer products from an assembly:

```bash
mlvamaps extract-amplicons \
  --input assembly.fasta \
  --primers examples/seer_lab_Ba/mlvamaps_primers.example.tsv
```

mlvamaps uses 32 threads by default. Pass `-t N` or `--threads N`;
`--threads 0` uses all available CPUs. Use `--quiet` to suppress progress.
External executables can be overridden with `--minimap2-bin`,
`--mafft-bin`, `--raxml-ng-bin`, `--epa-ng-bin`, and `--dnadiff-bin`.
RAxML-NG uses its `DNA` model-selection set by default to choose a nucleotide
model independently for each locus; override it with `--raxml-model`.

## Motivation and recognition

mlvamaps was created to make microbial MLVA data faster to analyze and easier
to inspect across laboratories, organisms, and sequencing approaches. Its
design was influenced by
[MLVA_finder](https://github.com/i2bc/MLVA_finder) and
[Sassy](https://github.com/RagnarGrootKoerkamp/sassy). Earlier mlvamaps
versions used Amplirust as an external in-silico PCR backend; the built-in
compatibility engine replaces that dependency.
