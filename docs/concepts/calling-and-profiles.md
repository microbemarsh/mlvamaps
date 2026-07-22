# Allele calling, profiles, and novelty

## Read-level predictions

Every read in a retained VSEARCH cluster receives:

- Best integer or half-unit repeat count and probability.
- Alternate repeat count and probability.
- Its unrounded repeat-count measurement and quality-derived uncertainty.
- Variant membership.
- Representative-relative substitution and indel counts.
- Evidence weight.

Reads with more edits receive less weight. These records are written to
`read_level_allele_predictions.tsv`.

## Locus posterior

The novel caller evaluates each unrounded read measurement against an explicit
half-unit allele grid with a Gaussian error model. Read quality, primer/flank
agreement, and representative-alignment identity control the contribution of
each read. The normalized sum is an abundance-like allele probability
distribution, so depth supports the modal allele without averaging a genuine
two-allele mixture into a false intermediate call.

The full ranked distribution is written to
`allele_probability_distribution.tsv`; the compact representation is retained
in `allele_calls.tsv`. A midpoint or a genuinely split distribution is marked
`AMBIGUOUS` rather than being resolved by a rounding tie-break.

`allele_calls.tsv` includes raw read depth and effective read depth. Effective
depth is the sum of edit-aware evidence weights and can be lower than the
number of contributing reads.

## Status logic

A locus is `LOW_DEPTH` when it has fewer reads than `--min-depth`. With
enough reads, it is `AMBIGUOUS` when the top posterior is below
`--min-posterior` or leads the second-best count by less than 0.2.

An otherwise decisive call is `OUT_OF_RANGE` when it lies outside the panel's
expected count range. It is `MULTIPLE_VARIANTS` when at least two variants pass
the EM mixture threshold and the dominant estimated fraction is below 0.8.
Raw low-count clusters that the mixture model classifies as trace evidence do
not force this status.

`allele_calls.tsv` records both the raw retained-ASV count and the number of
EM-meaningful variants, plus the dominant estimated fraction. See
[variant mixture abundance](variant-mixtures.md) for the model and thresholds.

## Fingerprints

`mlva_fingerprint.tsv` contains one row with one column per panel locus.
`mlva_fingerprint_probabilistic.tsv` contains long-form repeat counts and
posterior values.

Assembly calls use the same fingerprint shape and half-unit grid. With one
product, its exact length is evaluated using a one-base-resolution measurement
model. With multiple products and FASTQ/BAM support, their mapped-read counts
weight the allele distribution and select the supported product. Assembly
confidence therefore reflects the evidence instead of always being 1.0.

## Profile matching

For every profile and every locus with both a called and expected value,
MLVAMaps calculates the absolute repeat-count difference. Exact values count as
matches; nonzero values contribute to distance and appear in
`mismatched_loci`.

`profile_matches.tsv` contains the closest 20 rows sorted by total distance
and then matching-locus count. Confidence is the fraction of compared loci that
match exactly.

Profiles from different panels or repeat-number conventions should not be
compared without a documented conversion.

## Novelty summary

The novelty score combines:

- Distance from the nearest supplied profile.
- Low posterior calls.
- Non-`PASS`/`LOW_DEPTH` statuses.

The interpretation is:

- Below 0.25: `known-like`
- 0.25 to below 0.6: `uncertain`
- 0.6 or above: `potentially novel profile`

Without a supplied profile database, the score uses a neutral distance
component. Novelty is a prioritization summary, not proof that a strain,
lineage, or allele is biologically novel.
