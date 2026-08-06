# Representative mapping and SNP evidence

mlvamaps adds mapping without introducing an external species reference.
Instead, each sample supplies its own SPOARS-assembled locus references.

## Reference selection

Competitive recruitment maps reads to locus and repeat-product classes.
mlvamaps selects the dominant mapped product group and builds a SPOARS POA
consensus. The standard assembly PCR caller must resolve that consensus before
it becomes the primary mapping reference.

Two representative files serve different purposes:

- `mapped_variant_representatives.fasta.gz` contains diagnostic repeat
  representatives for mapped product groups.
- `locus_mapping_references.fasta.gz` contains the complete primer-oriented
  assembly-PCR-resolved POA product at each locus.

## Mapped-group read alignment

Before locus-wide SNP mapping, Parasail globally aligns each unique repeat
sequence in a mapped product group to its diagnostic repeat representative.
These exact end-to-end tracebacks populate `mapped_read_memberships.tsv` with
substitutions, insertions, deletions, and edit distance.

Global Needleman-Wunsch alignment consumes both complete repeat sequences, so
this stage has no local-alignment or soft-clipping path.

## Dominant-locus mapping

All usable reads assigned to a locus are written to the internal mapping FASTQ,
including low-depth observations. minimap2 maps those reads against the
collection of assembly-PCR-resolved POA products. mlvamaps accepts a primary
alignment as evidence only when it maps to that read's assigned locus
representative and passes the MAPQ threshold.

The SAM is retained as `locus_read_alignments.sam`.

## Coverage

For accepted alignments, mlvamaps walks aligned query/reference base pairs.
Insertions and deletions do not contribute a base at a reference position.
Bases below `--min-base-quality` and non-ACGT observations are excluded from
the depth used by the SNP evidence stage.

`locus_mapping_summary.tsv` reports:

- Total usable reads.
- Accepted mapped reads and mapping rate.
- Mean mapping quality.
- Mean quality-filtered depth across the reference.
- Covered bases and coverage percentage.
- Number of retained SNP rows.

## SNP evidence

At every A/C/G/T reference position, an alternate A/C/G/T observation is
reported when it passes all of:

- Minimum accepted mapping quality.
- Minimum base quality.
- Minimum total depth.
- Minimum alternate-read count.
- Minimum alternate frequency.

`locus_snps.tsv` includes the locus, reference variant, 1-based position,
reference and alternate bases, depth, alternate depth/frequency, and mean
alternate-base quality.

Multiple alternate bases can produce separate rows at one position when each
passes the thresholds.

## Coordinate and interpretation limits

Positions are relative to the sample-derived POA amplicon. They are
not chromosome coordinates, even when optional coordinates exist in the panel.

This lightweight table makes read evidence inspectable and useful for
within-sample diversity or troubleshooting. It is not a whole-genome variant
caller, clinical VCF, or substitute for mapping all reads to a validated genome
reference.

Indels are reported separately in `mapped_read_memberships.tsv`, where global
repeat-region alignments preserve insertions and deletions within each mapped
product group.
