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

The primary FASTQ caller maps molecules to observed database loci and measures
primer-bounded products, resolved repeat flanks, or supported local assemblies.
Only measured molecule counts contribute the primary allele likelihoods. Reference mapping
without a measurement establishes presence. Measurements can introduce alleles
outside the reference catalog. A separate fallback reports an explicitly labelled
`ESTIMATED` count when direct sizing is unresolved: first use competing measured
observations, then flank-supported repeat-spanning CIGAR evidence, then a
reference-conditioned approximation from partial mappings. Estimates do not
become measured molecule evidence for reference classification.

Independent molecules contribute to the locus posterior. Duplicate reference
contexts and two mates from one fragment do not multiply support. Conflicting
measurements of one molecule remain unresolved. The shared primary results are
recorded in `common_locus_calls.tsv` and projected into `allele_calls.tsv`,
`calls.tsv`, and fingerprints. SPOARS and read-level predictions remain
confirmation outputs; they do not replace the shared primary call.

## Status logic

Primary statuses distinguish called, low-coverage, ambiguous, mixed,
detected-but-unresolved, and not-found loci. `--min-depth` defaults to one for
long reads and three for Illumina. Low-depth measured counts remain visible.
A call is ambiguous when the leading probability fails its configured threshold
or leads the next allele by less than 0.2.

Mixture detection uses independently measured molecule winners. A secondary
allele needs at least 20% of informative molecules and at least two molecules
for long reads or three for Illumina. Reference classification separately uses
isolate support or mixture EM according to `--sample-mode`.

`ESTIMATED` retains the best available approximation at low confidence, including
ambiguous measured calls. `estimation_method`, `estimation_reason`, and
`estimate_support` explain its source, uncertainty, and molecule support.
Partial-mapping estimates are limited to supported catalog alleles. Duplicate
reference matches do not multiply support, and repeat-only cross-mappings do
not qualify for an alignment estimate.

`PRESENT_COUNT_UNKNOWN` retains presence with an empty count when neither direct
sizing nor estimation is possible. Missing data are
never zero-filled. `presence_support` records detected molecules and
`molecule_support` records informative molecules in the common table. The
probabilities describe model support, not a calibrated probability of correctness.

## Fingerprints

`mlva_fingerprint.tsv` contains one row with one column per panel locus.
`mlva_fingerprint_probabilistic.tsv` contains long-form repeat counts and
posterior values.
Estimated counts are included; consult the call status and estimation fields
before treating a fingerprint value as measured. The HTML report presents every
panel locus in an always-visible count table, marking estimates with `≈`.

Assembly calls use the same fingerprint shape and historical integer/half-unit
rounding convention. Product selection follows the MLVA_finder-compatible rule
described in the assembly workflow. Optional FASTQ/BAM mapping reports depth
for the selected product but does not choose a different product or alter its
size-derived repeat count; the assembly call confidence is therefore `1.0`.

## Profile matching

For every profile and every locus with both a called and expected value,
mlvamaps calculates the absolute repeat-count difference. Exact values count as
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

When per-locus FASTQ probability distributions are available, mlvamaps also
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
