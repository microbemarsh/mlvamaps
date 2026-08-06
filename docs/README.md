# mlvamaps documentation

mlvamaps separates the typing result from the evidence used to make it. The
root README provides the quick start; this documentation explains each input
path, intermediate stage, output, and interpretation boundary.

## Workflows

- [FASTQ and amplicon sequencing](workflows/fastq.md)
- [Paired/single-end Illumina sequencing](workflows/illumina.md)
- [Assemblies and assembly read support](workflows/assembly.md)
- [Reference building](workflows/reference-building.md)

## Concepts

- [Variant mixture abundance](concepts/variant-mixtures.md)
- [Repeat-aware SNP placement and phylogeography](concepts/repeat-snp-phylogeography.md)
- [Representative mapping and SNP evidence](concepts/representative-mapping.md)
- [Allele calling and profiles](concepts/calling-and-profiles.md)

## Reference

- [Input and panel formats](reference/input-formats.md)
- [Output files](reference/outputs.md)
- [CLI options and thresholds](reference/cli.md)

## Guides

- [Adapting a panel for a new organism or MLVA scheme](guides/new-organism-panel.md)

## Scope

mlvamaps is organism-agnostic and panel-driven. It can assess microbial MLVA
data when the user supplies the primers and interpretation metadata for that
scheme. It does not infer an unknown MLVA panel from a genome, and it does not
claim that profiles from different laboratories or schemes are interchangeable.

The FASTQ path is best suited to accurate amplicon or long reads that span both
primers. The assembly path is best suited to assembled products and can add
accurate-read support without requiring individual reads to span the complete
locus. Noisy long-read mapping is outside the supported scope.
