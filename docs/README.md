# mlvamaps documentation

mlvamaps separates the typing result from the evidence used to make it. The
root README provides the quick start; this documentation explains each input
path, intermediate stage, output, and interpretation boundary.

## Workflows

- [Architecture and resource flow](workflows/architecture.md)
- [Accurate long-read and amplicon FASTQ](workflows/fastq.md)
- [Paired/single-end Illumina sequencing](workflows/illumina.md)
- [Genome assemblies and optional read support](workflows/assembly.md)
- [Reference building](workflows/reference-building.md)
- [Dataset aggregation and MYOGA export](workflows/myoga-export.md)

## Concepts

- [Variant mixture abundance](concepts/variant-mixtures.md)
- [Repeat-aware SNP placement and phylogeography](concepts/repeat-snp-phylogeography.md)
- [Representative mapping and SNP evidence](concepts/representative-mapping.md)
- [Allele calling and profiles](concepts/calling-and-profiles.md)
- [Acceleration backends and threading](concepts/acceleration.md)

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

The accurate-long/amplicon pathway operates directly on FASTQ reads and can
retain locus-presence evidence when a read does not span the complete product.
The assembly pathway operates on FASTA contigs and can add read-depth support
without changing the selected assembly product or its repeat count. Noisy
long-read mapping is outside the supported scope.
