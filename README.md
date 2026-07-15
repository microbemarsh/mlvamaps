# mlvamaps

mlvamaps is a panel-driven toolkit for analyzing microbial MLVA/VNTR typing
data from accurate amplicon or long reads containing the target products,
assemblies, and assemblies with supporting accurate reads.

It is designed to be organism-agnostic: the same pipeline can be used for any
microbe with a defined MLVA primer panel and enough repeat metadata to interpret
the products. MLVAMaps does not hard-code one species, typing scheme, or profile
database.

mlvamaps returns conventional repeat-copy-number fingerprints while preserving
the sequence evidence behind them:

- Primer-supported locus detection in compressed or uncompressed FASTQ.
- Assembly in-silico PCR in either orientation.
- Fast, gap-aware VSEARCH clustering of observed VNTR sequences.
- Per-read substitution and indel evidence without replacing observed reads
  with a consensus.
- minibwa mapping to sample-derived representative amplicons and
  reference-relative SNP evidence.
- Optional read-depth support for assembly products.
- Probabilistic locus calls, profile matching, novelty summaries, generated gel
  comparisons, and a self-contained HTML report.

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
[Amplirust](https://github.com/erdikilic/amplirust),
[VSEARCH](https://github.com/torognes/vsearch),
[minibwa](https://github.com/lh3/minibwa), and the Python/native sequence
libraries declared in `environment.yml`.

## Quick start

Analyze amplicon or other primer-spanning reads:

```bash
mlvamaps call primers.tsv sample.fastq.gz
```

Analyze an assembly:

```bash
mlvamaps call primers.tsv assembly.fasta
```

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

Results are written to `results/` by default. Start with:

- `calls.tsv` for compact per-locus calls.
- `mlva_fingerprint.tsv` for the conventional wide fingerprint.
- `report.html` for the visual summary.

FASTQ runs additionally provide VSEARCH variants, read memberships, minibwa
alignments, mapping coverage, and SNP evidence. Assembly runs provide extracted
primer products and optional read support.

## Supported data

| Input | What mlvamaps assesses |
| --- | --- |
| Amplicon FASTQ/FASTQ.GZ | Primer-supported VNTR reads, sequence variants, repeat counts, representative mapping, and SNP evidence. |
| Accurate primer-spanning WGS reads | Reads containing both required primers and a valid product are analyzed through the same FASTQ path. |
| Assembly FASTA | In-silico primer products, product coordinates, sizes, and repeat counts. |
| Assembly plus accurate FASTQ | Assembly calls plus minibwa read count and mean coverage for extracted products. |
| Assembly plus SAM/BAM | Assembly calls plus overlap-based read support from existing alignments. |
| Known profile TSV | Closest MLVA profiles, mismatched loci, distance, confidence, and novelty context. |

FASTQ assignment requires both primers to occur in a valid product. Ordinary
shotgun reads that do not span the complete target amplicon are therefore not
treated as locus calls. Assembly mode is the appropriate route when the target
is represented across multiple non-spanning reads. The mapping paths are scoped
to accurate reads; noisy long-read mapping is not supported.

## How it works

For FASTQ data, mlvamaps:

1. Filters reads by length and quality.
2. Uses Amplirust to pair degenerate primers and orient each product.
3. Locates the repeat region and measures repeat/motif evidence.
4. Dereplicates and clusters reads by locus with VSEARCH.
5. Preserves observed cluster representatives and annotates per-read edits.
6. Maps locus reads to the dominant observed amplicon with minibwa.
7. Reports quality-filtered, representative-relative SNP evidence.
8. Combines read probabilities into per-locus repeat-count calls.
9. Builds the fingerprint, compares profiles, scores novelty, and writes HTML.

For assemblies, mlvamaps:

1. Finds paired-primer products with Amplirust.
2. Selects valid products and converts size into repeat count where the panel
   provides enough metadata.
3. Optionally adds minibwa or existing SAM/BAM read support.
4. Builds the same fingerprint, profile comparison, and report formats.

The minibwa mapping coordinates are positions within the sample-derived
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
- [Allele calling, profiles, and novelty](docs/concepts/calling-and-profiles.md)
- [Adding a new organism or MLVA scheme](docs/guides/new-organism-panel.md)

## Additional commands

Simulate amplicon reads for pipeline testing:

```bash
mlvamaps simulate \
  --loci examples/mlva_loci.example.tsv \
  --sample-id SIM1 \
  --depth 500 \
  -o simulated
```

Export Amplirust products from an assembly:

```bash
mlvamaps extract-amplicons \
  --input assembly.fasta \
  --primers examples/seer_lab_Ba/mlvamaps_primers.example.tsv
```

MLVAMaps uses 32 threads by default. Pass `-t N` or `--threads N`;
`--threads 0` uses all available CPUs. Use `--quiet` to suppress progress.

## Motivation and recognition

mlvamaps was created to make microbial MLVA data faster to analyze and easier
to inspect across laboratories, organisms, and sequencing approaches. Its
design was influenced by
[MLVA_finder](https://github.com/i2bc/MLVA_finder) and
[Amplirust](https://github.com/erdikilic/amplirust).
