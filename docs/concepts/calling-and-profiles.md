# Allele calling, profiles, and novelty

## Read-level predictions

Every read in a retained VSEARCH cluster receives:

- Best repeat count and probability.
- Alternate repeat count and probability.
- Variant membership.
- Representative-relative substitution and indel counts.
- Evidence weight.

Reads with more edits receive less weight. These records are written to
`read_level_allele_predictions.tsv`.

## Locus posterior

The Bayesian caller accumulates best and alternate read probabilities across
the configured repeat-count range. It normalizes the evidence to produce the
best and second-best posterior.

`allele_calls.tsv` includes raw read depth and effective read depth. Effective
depth is the sum of edit-aware evidence weights and can be lower than the
number of contributing reads.

## Status logic

A locus is `LOW_DEPTH` when it has fewer reads than `--min-depth`. With
enough reads, it is `AMBIGUOUS` when the top posterior is below
`--min-posterior` or leads the second-best count by less than 0.2.

An otherwise decisive call is `OUT_OF_RANGE` when it lies outside the panel's
expected count range. It is `MULTIPLE_VARIANTS` when multiple retained
variants exist and the dominant variant frequency is below 0.8.

## Fingerprints

`mlva_fingerprint.tsv` contains one row with one column per panel locus.
`mlva_fingerprint_probabilistic.tsv` contains long-form repeat counts and
posterior values.

Assembly calls use the same fingerprint shape. Assembly `PASS` calls receive
confidence 1.0 because the count is derived from an observed product, while
unknown or absent calls receive 0.0.

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
