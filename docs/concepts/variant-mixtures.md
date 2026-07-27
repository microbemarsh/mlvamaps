# Variant mixture abundance

VSEARCH deliberately retains observed clusters. At deep coverage, a locus can
therefore contain a dominant variant, genuine secondary variants, and small
clusters compatible with sequencing or amplification error. Treating every
retained cluster equally makes `MULTIPLE_VARIANTS` difficult to interpret.

MLVAMaps adds a count-based expectation-maximization layer inspired by
[Emu](https://github.com/treangenlab/emu), which uses alignment likelihoods and
an abundance-dependent prior to estimate microbial composition. MLVAMaps does
not run Emu's taxonomic workflow or require a taxonomy database; it adapts the
same iterative abundance idea to the observed VNTR representatives at one
locus.

## Model inputs

For each locus, the model uses:

- The retained VSEARCH support count for every variant.
- The observed representative repeat sequence for every variant.
- The substitution and indel totals from Parasail cluster-member alignments.

The within-cluster edits provide a smoothed locus error-rate estimate. Exact
Parasail alignments between representative sequences convert their edit
distances into relative observation likelihoods. Very different
representatives consequently have little assignment ambiguity; close
representatives can borrow support according to their current abundance.

## EM iterations

The model starts from smoothed VSEARCH count fractions. Each iteration:

1. Calculates the probability that each observed cluster count arose from each
   candidate variant using sequence likelihood multiplied by current abundance.
2. Distributes the cluster count across candidates using those probabilities.
3. Normalizes the distributed counts to update variant fractions.

Iterations stop when both abundance and likelihood changes converge or after
the iteration limit. This operates on aggregated counts rather than expanding
millions of identical read assignments, so runtime depends mainly on the number
of retained variants per locus.

## Floors and meaningful variants

After initial convergence, MLVAMaps applies Emu's depth-adaptive component
floor: `1 / (reads + 1)` through 1,000 retained reads and `10 / reads` above
1,000 reads. Components below that floor are removed and the model is refit.

The final `--min-mixture-fraction` threshold, 0.01 by default, determines which
remaining variants have enough abundance for interpretation. A secondary
variant must also have at least `--min-secondary-reads` observations, two by
default, before it can alter mixture status. The most abundant component is
always retained. Evidence tiers are:

- `DOMINANT`: highest estimated fraction.
- `CONFIRMED_SECONDARY`: passes both abundance and read-support thresholds.
- `CANDIDATE`: passes the abundance threshold but lacks independent read
  support; retained for rapid low-coverage detection without changing the
  primary signature.
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
