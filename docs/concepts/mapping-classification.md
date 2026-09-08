# Mapping-based VNTR classification and profile trees

`mlvamaps call --database DB` uses original alignment evidence for reference
classification. It does not average the nearest reference distances within a
taxon. The mixture calculation is inspired by
[Emu's alignment-likelihood and EM approach](https://github.com/treangenlab/emu),
with a separate VNTR repeat-length term and molecule-level accounting.
This is an independent implementation using existing NumPy, minimap2, pysam,
and parasail dependencies; installing Emu is unnecessary.

## Evidence and scoring

Read inputs are mapped competitively to **observed reference amplicons**. This is
a second mapping pass: allele-calling contexts contain synthetic repeat states
and their evidence projection discards competing sequence backgrounds. Those
selected or synthetic sequences are not used as query observations here.
Identical observed locus sequences share one mapping target and retain all
reference IDs. The classifier requests secondary competitors up to the number
of targets, with minimap2's secondary score-ratio filter disabled; candidates
that the mapper does not report are not assumed to have a measured alignment.

For assemblies, Sassy-backed in silico PCR recovers primer-oriented products.
The selected observed products are aligned end-to-end to corresponding reference
loci with Parasail. Each product contributes one locus observation. The reported
repeat-count tables remain unchanged; reference-relative mismatch counts, CIGARs
and repeat differences are retained in `assembly_reference_evidence.tsv`, with
query sequences in `query_amplicons.fasta`, under `classification/`.

For each retained alignment:

- Mismatches, including repeat substitutions, and indels outside the repeat
  interval contribute sequence-error terms. Error rates are estimated from the
  best alignment per read using pseudocounts and bounds, not from MAPQ.
- At least six aligned bases outside the repeat are required. Alignments that
  leave more nonrepeat sequence unexplained than the competing alignments are
  penalized for that deficit.
- The net indel length inside the repeat interval is excluded from the sequence-gap term; opposing indels with no net length change still contribute sequence errors.
  When the alignment spans the repeat with at least three flanking bases on
  each side, the net repeat-length difference contributes
  `-abs(repeat_delta) / classification_repeat_scale` to the log score.
  Partial coverage and paired geometry alone do not establish repeat length.
- Both mates contribute sequence evidence, but their repeat-length disagreement
  is scored once per molecule. Duplicate alignments do not add observations.
  Molecules tied between different loci are excluded rather than counted twice.

For an isolate, scores are summed across observations **for the same reference**.
Each molecule receives weight `min(locus_depth, 20) / locus_depth`, limiting a
single deep locus to twenty effective observations. References with identical
likelihood columns share one component and one prior. A softmax of the summed
scores gives the relative model weights; indistinguishable references spanning
different taxa remain taxonomically ambiguous.

`--sample-mode metagenome` runs weighted EM over those components. One component
fraction is shared across loci. EM alternates normalized molecule assignments
and fraction updates until fractions change by less than `1e-8`, with a limit
of 500 iterations. The objective history and convergence status are recorded.
These depth-capped fractions describe the retained marker evidence; they are
**not estimates of organism abundance** without a locus sampling/copy-number
model. The historical default sample mode remains `metagenome`; use
`--sample-mode isolate` when one source is expected. Assemblies use isolate mode.

An explicit unclassified component has log score `-max(4, 0.1 * aligned_flank_bases)`
per molecule. An unreported reference alignment receives the molecule's lowest
reported score minus eight. These are provisional background/floor settings,
not calibrated novel-taxon probabilities. The error-rate pseudocounts are 100
matches and one of each error class, bounded to `[1e-4, 0.25]`.

## Coverage and interpretation

Low depth does not impose a new recovery-fraction requirement. Informative
sequence evidence can support a closest taxon without a resolved repeat count.
The optional missing-locus heuristic retains its defaults: three informative
molecules at 80% of panel loci, then one penalty unit for an undetected query
locus recorded in that reference. A locus detected by the new mapping pass is
not penalized even if the allele caller had no call. The penalty is subtracted
once from each reference's joint log score; in mixture mode it is distributed
across the weighted observations. It is not evidence of confirmed absence.

Identification follows the highest-ranked **reference group**, using exactly the
same order and support as `mapping_reference_matches.tsv`. Support is never
summed across a species to choose the identification. The `assignment` names the taxon of the
highest-ranked reference group, with its reference IDs retained separately for
traceability. Groups spanning multiple taxa are explicitly ambiguous; missing
metadata do not produce an invented taxon call. Reference identification works without
taxon metadata. The evidence table has one row per reference group.

A model-supported unique reference requires support of at least 0.9, convergence,
at least two catalog references for comparison, and two observed loci by default
(`--taxon-min-loci` overrides the locus count). A strongly supported group of
indistinguishable references is `AMBIGUOUS_REFERENCES`, even when all its members
have the same taxon label. Multiple reference groups with at least 0.05 fraction
in metagenome mode are `MIXED_REFERENCES`, including mixtures within one species;
the headline names the taxon of the highest-ranked group, with other components listed
in the evidence table. Below the support requirement the closest group remains
visible as `CLOSEST_REFERENCE_LOW_CONFIDENCE`. No observations produce
`INSUFFICIENT_EVIDENCE` and an unresolved assignment.

These are model summaries, not calibrated identification probabilities. Mixtures
of distant VNTR loci can also be indistinguishable from an unrepresented mosaic
profile without molecule linkage. The unclassified fraction measures the explicit
background component; ambiguity between named references is represented by the
listed equivalent IDs. Unmapped input reads are outside the fitted matrix.

Distance-based ranking, phylogenetic placement, and calibrated target testing
are not part of the typing workflow.

## Newick trees for MYOGA

Classification and profile similarity answer different questions. EM mixture
fractions and query likelihoods are not symmetric sample-to-sample distances.
They are therefore **not used to construct the tree**.

Each classified sample with at least two numeric, nonmixed repeat calls can
produce `classification/mlva_profiles.tree`. It is a neighbor-joining Newick
tree, using the existing implementation, for the query and up to twenty ranked
reference profiles complete at every selected query locus. All pairs use that
same locus set. Branch-length units are mean absolute repeat-count differences.
Missing calls are not imputed to zero; a status table explains when no tree can
be written. Matching metadata and the distance matrix are also exported.

This is a **phenetic MLVA profile-similarity tree**, not an inferred evolutionary
phylogeny, transmission tree, or time-scaled genealogy. Use it within comparable
panels and taxonomic groups. Repeated profile changes need not imply common
ancestry, and classification rank is not a tree branch length.

For comparing many samples, reuse the existing export rather than rebuilding
reference phylogenies:

```bash
mlvamaps export-myoga --results results/ --metadata metadata.tsv \
  --metadata-id sample_id --distance repeat \
  --min-pairwise-loci 3 --min-pairwise-fraction 0.8 -o myoga_mlva/
```

This writes `mlva_nj.tree` and `mlva_metadata.tsv`, plus matrices, overlap and
exclusion tables. The three-locus/80% example is adjustable to the panel. Low
coverage classification can remain available even when the evidence is too
sparse for a meaningful profile tree. Unlike the per-sample tree, the multi-sample
export uses shared calls per pair and explicitly removes unsupported overlaps.
The `.tree` files are Newick; MYOGA versions may additionally require a GGR
file or a tree-only import mode. No synthetic Parsnp/GGR data is invented.

MAFFT, RAxML-NG and EPA-ng are not required for mapping classification or these
repeat-profile trees. Reference builds and typing do not invoke sequence-tree tools.

## Validation and limits

Regression tests cover repeat-length-only differences, sequence errors,
partially observed repeats, paired molecules, duplicate references, cross-locus
ambiguity, a monotone EM objective, no-match evidence and Newick output. An
optional real-minimap2 test checks synthetic isolates at one and three molecules
per locus and a known 70:30 mixture. These checks verify implementation behavior;
they do not replace held-out, taxon-specific validation across sequencing
technologies, missing-reference scenarios and uneven locus recovery.

The likelihood matrix is currently dense, like the existing alignment evidence
path. Large catalogs/high-depth inputs can require substantial memory. Held-out validation on representative data is still required.

Manifest batches with only legacy classification outputs are rerun with the new classifier. Use `--force` when changing settings on already completed mapping-classification results.
