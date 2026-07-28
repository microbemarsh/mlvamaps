# Allele calling and profiles

## Read-level predictions

Every complete read in a mapping-derived product group receives:

- Best integer or half-unit repeat count and probability.
- Alternate repeat count and probability.
- Its unrounded repeat-count measurement and quality-derived uncertainty.
- Variant membership.
- Representative-relative substitution and indel counts.
- Evidence weight.

For the default assembly-equivalent FASTQ convention, the raw observation is
calculated from the complete primer-bounded product using the same panel
calibration as an assembly product. Historical assembly rounding supplies the
measurement center used by the probabilistic caller. Both the raw value and
the rounded measurement center are retained for auditability.

Reads with more edits receive less weight. These records are written to
`read_level_allele_predictions.tsv`.

## Locus posterior

The caller first identifies the dominant mapped product group and builds a
SPOARS consensus. The standard assembly PCR caller defines the primary allele
from that product. Read quality, primer/flank agreement, and group-alignment
identity control confidence contributions. Effective confidence depth is
capped at 25 by default to limit overconfidence from correlated amplification
or systematic errors.

Secondary clusters never average the primary call toward a false intermediate
allele. They are retained and interpreted separately by the mixture model.

The full ranked distribution is written to
`allele_probability_distribution.tsv`; the compact representation is retained
in `allele_calls.tsv`. A midpoint or a genuinely split distribution is marked
`AMBIGUOUS` rather than being resolved by a rounding tie-break.

`allele_calls.tsv` includes total and primary-cluster read depth, their
edit-aware effective depths, capped confidence depth, confirmed and candidate
variant counts, and a compact secondary-allele representation.

## Status logic

A locus is `LOW_DEPTH` when its dominant cluster has fewer reads than
`--min-depth`, which defaults to one. A single informative molecule is
retained as `SINGLE_MOLECULE_PROVISIONAL`. With enough primary reads, a locus
is `AMBIGUOUS` when the top posterior is below `--min-posterior` or leads the
second-best count by less than 0.2.

An otherwise decisive call is `OUT_OF_RANGE` when it lies outside the panel's
expected count range. It is `MULTIPLE_VARIANTS` when at least two variants pass
the EM mixture threshold. Metagenome interpretation is the default and flags
every meaningful secondary variant; explicit isolate mode requires the
dominant estimated fraction to be below 0.8. Raw low-count clusters that the
mixture model classifies as trace evidence do not force this status.

`allele_calls.tsv` records both the raw mapped-group count and the number of
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
`mismatched_loci`. Both workflows rank by this distance and matched-locus count;
FASTQ probabilities only break otherwise equal matches. Unobserved loci are
not imputed.

Let \(q_\ell\) be the query repeat count at locus \(\ell\), \(p_\ell\) the
profile repeat count, and
\(C(q,p)=\{\ell:q_\ell\text{ and }p_\ell\text{ are both observed}\}\).
The reported MLVA distance is the Manhattan (L1) distance on repeat counts:

```latex
D_{\mathrm{MLVA}}(q,p)
  = \sum_{\ell \in C(q,p)} \left|q_\ell-p_\ell\right|.
```

A one-repeat disagreement therefore contributes 1, while a half-repeat
disagreement contributes 0.5. Missing loci do not contribute a zero or a
penalty; they are listed in `uncompared_loci`. For that reason, distance should
be interpreted together with `compared_loci` and `recovered_loci`.

The exact-match confidence reported alongside distance is:

```latex
\operatorname{confidence}(q,p)
  = \frac{\sum_{\ell \in C(q,p)}
      \mathbf{1}\!\left[q_\ell=p_\ell\right]}
    {|C(q,p)|}.
```

When per-locus FASTQ probability distributions are available, MLVAMaps also
calculates the mean negative log likelihood of the profile:

```latex
\overline{\operatorname{NLL}}(q,p)
  = -\frac{1}{|C(q,p)|}
      \sum_{\ell \in C(q,p)}
      \log\!\left(\max(P_\ell(p_\ell),10^{-9})\right).
```

Profiles are sorted by distance ascending, exact matching loci descending,
mean NLL ascending, and profile ID. The NLL is a tie-breaker; it does not
change the repeat-count distance.

`profile_matches.tsv` contains every ranked reference match. The `match_type`
column distinguishes direct `mlva_profile` comparisons from
`sequence_reference` matches produced by `--database`; the latter retain their
combined SNP/repeat marker fields from the HTML report. Direct profile rows
include aggregate comparison fields and compact per-locus allele vectors.
`profile_match_loci.tsv` expands those comparisons to one machine-readable row
per profile and locus.

Profiles from different panels or repeat-number conventions should not be
compared without a documented conversion.
