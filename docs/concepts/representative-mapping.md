# Representative mapping and SNP evidence

MLVAMaps adds mapping without introducing an external species reference.
Instead, each sample supplies its own observed locus references.

## Reference selection

VSEARCH returns one or more retained sequence variants per locus. MLVAMaps
selects the variant with the greatest read support. Its representative is the
observed read chosen by VSEARCH, not a generated consensus.

Two representative files serve different purposes:

- `vntr_asv_representatives.fasta` contains the repeat region for every
  retained variant.
- `locus_mapping_references.fasta` contains the complete primer-oriented
  amplicon for the dominant variant at each locus.

## Cluster-member alignment

Before locus-wide SNP mapping, minibwa separately indexes every retained
VSEARCH representative amplicon and maps that cluster's reads back to it.
MLVAMaps walks each SAM CIGAR across the recorded query and representative
repeat boundaries to reconstruct the gapped repeat alignment. These alignments
populate `vntr_asv_memberships.tsv` with substitutions, insertions, deletions,
and edit distance relative to the read's own selected cluster representative.

This stage requires both repeat regions to be fully spanned. MLVAMaps stops
with a clear error instead of silently treating a soft-clipped repeat as an
exact alignment.

## Dominant-locus mapping

All usable reads assigned to a locus are written to the internal mapping FASTQ,
including reads that did not enter a retained cluster. minibwa maps those reads
against the collection of dominant amplicons. MLVAMaps accepts a primary
alignment as evidence only when it maps to that read's assigned locus
representative and passes the MAPQ threshold.

The SAM is retained as `locus_read_alignments.sam`.

## Coverage

For accepted alignments, MLVAMaps walks aligned query/reference base pairs.
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

Positions are relative to the sample-derived representative amplicon. They are
not chromosome coordinates, even when optional coordinates exist in the panel.

This lightweight table makes read evidence inspectable and useful for
within-sample diversity or troubleshooting. It is not a whole-genome variant
caller, clinical VCF, or substitute for mapping all reads to a validated genome
reference.

Indels are reported separately in `vntr_asv_memberships.tsv`, where global
repeat-region alignments preserve insertions and deletions relative to each
retained representative.
