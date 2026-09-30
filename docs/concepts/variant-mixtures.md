# Variant mixture abundance

Competitive mapping can assign locus reads to multiple repeat-product classes.
At deep coverage, these may represent a dominant allele, genuine secondary
alleles, or small error-compatible groups.

mlvamaps adds a count-based expectation-maximization layer inspired by
[Emu](https://github.com/treangenlab/emu), which uses alignment likelihoods and
an abundance-dependent prior to estimate microbial composition. mlvamaps does
not run Emu's taxonomic workflow or require a taxonomy database; it adapts the
same iterative abundance idea to the observed VNTR representatives at one
locus.

## Model inputs

For each locus, the model uses:

- The mapped-read support count for every product group.
- A diagnostic repeat representative for every group.
- Substitution and indel totals from Parasail group-member alignments.

The within-cluster edits provide a smoothed locus error-rate estimate. Exact
Parasail alignments between representative sequences convert their edit
distances into relative observation likelihoods. Very different
representatives consequently have little assignment ambiguity; close
representatives can borrow support according to their current abundance.

## EM iterations

The model starts from smoothed mapped-group count fractions. Each iteration:

1. Calculates the probability that each observed mapped-group count arose from each
   candidate variant using sequence likelihood multiplied by current abundance.
2. Distributes the cluster count across candidates using those probabilities.
3. Normalizes the distributed counts to update variant fractions.

Iterations stop when both abundance and likelihood changes converge or after
the iteration limit. This operates on aggregated counts rather than expanding
millions of identical read assignments, so runtime depends mainly on the number
of retained variants per locus.

## Floors and meaningful variants

No depth-adaptive component floor or minimum read-support cutoff is applied.
All observed components participate in the abundance fit.

The final `--min-mixture-fraction` threshold, 0.01 by default, determines which
variants have enough relative abundance for interpretation. A secondary
variant can alter mixture status with one observation; the legacy
`--min-secondary-reads` option is ignored. The most abundant component is
always retained. Evidence tiers are:

- `DOMINANT`: highest estimated fraction.
- `CONFIRMED_SECONDARY`: an observed secondary passing the abundance threshold;
  this legacy label does not imply replicated read support.
- `CANDIDATE`: retained for compatibility with older output files.
- `TRACE`: below the meaningful threshold.

In default metagenome mode, `MULTIPLE_VARIANTS` is assigned only when a
confirmed secondary remains. Candidate and trace variants never change the
primary allele or force mixture status. Isolate mode additionally requires the
dominant estimated fraction to be below 0.8. No tier is discarded from the
abundance TSV; the HTML report shows candidates separately and combines trace
components visually.

## Output

`vntr_mixture_abundance.tsv` reports observed and EM-estimated fractions,
estimated read counts, abundance class, evidence tier, thresholds, inferred
error rate, iteration count, and convergence for every retained variant. Fractions describe retained
variant-supporting reads at that locus, not organism abundance or absolute cell
counts.
