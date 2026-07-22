# Repeat-aware SNP placement and phylogeography

MLVAMaps treats tandem-repeat evolution and nucleotide substitution as linked
but distinct signals. Standard nucleotide models do not model VNTR copy-number
changes, and alignment gaps generally do not contribute like substitutions.
Combining both signals inside one unlabelled tree distance would therefore lose
biological meaning.

Reference builds persist canonical SHA-256 keys for the oriented unmasked
amplicon, repeat-masked SNP sequence, and combined SNP/repeat marker. Query
placement checks these keys first. A reference or tied reference group matching
at every callable locus receives distance zero and reuses its existing tree
location without launching MAFFT or EPA-ng. This is the default fast path;
novel and non-exact queries continue through alignment and phylogenetic
placement.

Reference-tree construction groups identical repeat-masked SNP sequences and
uses one deterministic representative ID per haplotype whenever at least two
distinct haplotypes remain. Placement distances are expanded back to every
member reference afterward. This removes redundant zero-length tree tips while
preserving all reference IDs in rankings, metadata joins, and reports.

For loci with defined flanks, MLVAMaps identifies and removes the repeat tract
before MAFFT reference alignment, RAxML-NG tree inference, and EPA-ng query
placement. It retains the removed tract as explicit evidence:

- raw and reported repeat count;
- repeat sequence and ordered repeat-unit haplotype;
- masking coordinates and boundary method;
- likelihood-weighted SNP-tree distance and placement entropy.

Each per-locus SNP comparison retains two complementary measurements. EPA-ng's
likelihood-weighted tip distance is divided by the median pairwise distance of
the reference tree. Direct query/reference divergence is measured in the
shared MAFFT alignment and divided by the median positive direct divergence
from the query at that locus. An exact repeat-masked SNP-sequence match is zero
by definition; otherwise the two normalized SNP measurements are averaged.
This prevents placement ambiguity from ranking a sequence below a reference
with which it is exactly identical while retaining tree context for non-exact
comparisons.

Repeat-count differences are divided by the reference repeat dispersion, with
a minimum half-unit scale. This prevents one fast-evolving locus or
intrinsically variable VNTR from dominating only because of its numeric scale.

The combined score is:

```text
(--phylogeny-snp-weight × summed normalized SNP distance)
+ (--phylogeny-repeat-weight × summed normalized repeat distance)
```

Both components remain in the output. Weights should be calibrated against
known epidemiological links when those data are available; equal weights are a
transparent starting point, not a universal biological constant.

If the EPA/tree-only score would rank a complete exact marker match below
another reference, `combined_marker_matches.tsv` records
`EXACT_MATCH_OVERRIDES_PLACEMENT` and the HTML report displays a warning. Tied
exact references retain the same rank because sequence evidence cannot
distinguish them.

Combined ranking requires every candidate reference to be present at all placed
loci and to have repeat information wherever the query has a repeat call. This
prevents missing marker components from producing an artificially small score.

`reference_metadata.tsv` can add collection date, coordinates, location, and
source to the combined ranking. `combined_markers.tree` is a plain Newick
neighbor-joining tree whose pairwise distances use the same normalized SNP and
repeat components. Its tip labels are the exact reference IDs plus the query
sample ID, allowing MYOGA to join the tree to external metadata. The long-form
locus table supports downstream analyses that distinguish repeat changes on
stable SNP backgrounds from SNP divergence within a shared MLVA profile.
