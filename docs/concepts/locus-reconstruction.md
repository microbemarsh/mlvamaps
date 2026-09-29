# Targeted paired-end locus reconstruction

Illumina FASTQ has one standard pathway: competitive locus recruitment, direct
measurement, targeted reconstruction, and haploid likelihood inference when a
complete locus cannot be recovered. Assembly and long-read recovery feed the
same `LocusProduct` and `genotype_product` interpreter.

## Recruitment and mate rescue

One streaming pass reads paired FASTQ, applies QC, and looks up both strands in
a combined hash index of panel anchors. Each locus contributes one panel
template; reference database sequences never supply or expand these templates.
Native SIMD flank alignments verify seed hits. A pair retains both mates even
when only one has an anchor; the second mate need not map at all. Full and
compact ambiguity audits are available with `--sr-recruitment-audit`.

Legacy primer-only panels use full primer matches with up to two edits (capped
for short primers), including IUPAC primer symbols. Mates without a primer are
oriented using their anchored partner. Short or degenerate primers bypass the
four-base seed filter when a valid match need not contain a concrete seed.
Only complete observed products provide length calls for these panels;
unknown internal flanks and motifs are not borrowed from reference sequences.

Assignment uses only non-repeat anchor scores. Multiple contexts of one locus
count as one competitor. A winning locus must exceed the runner-up by both 12
alignment-score points and 20% of its own score. Tied loci remain unassigned.
Repeat copy number and the frequency of a reference in the database do not
increase recruitment support. This is an in-process combined target index;
allele recovery does not launch minimap2 or create intermediate SAM files.

Only candidate molecules undergo motif testing. Cyclic templates and native
gapped alignment handle rotations, reverse complements, substitutions and
small indels. `--short-repeat-fraction` (default 0.7) sets the matched-base
fraction required for repeat-rich tagging. Unanchored repeat-only molecules
cannot identify a locus and are excluded rather than assigned to every locus
with that motif.

## Evidence definitions

Tags are nonexclusive, counted once per molecule:

| Tag | Deterministic condition |
| --- | --- |
| `FULL_SPAN` | A read has accepted anchors at both repeat boundaries, or a uniquely merged molecule contains both primers. |
| `LEFT_BOUNDARY` | A terminal left-flank anchor leads into motif-compatible sequence. |
| `RIGHT_BOUNDARY` | Motif-compatible sequence leads into an initial right-flank anchor. |
| `FLANK_PAIR` | Opposite mates anchor to opposite flanks with inward orientation. |
| `REPEAT_RICH` | At least one mate meets the motif-compatible fraction threshold. |
| `ANCHORED_REPEAT` | A repeat-rich mate belongs to a flank-anchored pair. |
| `SOFTCLIP_LEFT` | Sequence continues beyond a left boundary anchor without a right boundary anchor. |
| `SOFTCLIP_RIGHT` | Sequence precedes a right boundary anchor without a left boundary anchor. |
| `UNINFORMATIVE` | A recruited molecule meets none of the preceding conditions. |

The soft-clip categories describe clipping relative to a flank alignment, not
literal mapper CIGAR operations. Same-strand discordant pairs and inconsistent
spans are retained as diagnostics but excluded from recovery. An accepted
anchor has at least 20 aligned bases (or the entire shorter flank) and at most
15% edits. Anchors of 10–19 bases are rescued when adjacent repeat sequence
independently supports a boundary. Boundary anchors can extrapolate at most three bases.

## Direct and reconstructed sequence

Complete high-quality primer-bounded reads are measured first. Mate merging
requires inward orientation and exactly one exact overlap of at least 20 bases
that is not repeat-only. A periodic overlap cannot identify its own offset.

A read spanning both repeat boundaries can establish length even without both
primers. Unobserved flank tails are `N`, not reference bases. When other reads
can fill those tails, bounded reconstruction attempts a complete sequence.
Complete products retain SNPs and indels, and supported distinct sequences are
retained as possible mixtures. Singleton errors remain trace evidence; if only
singleton same-length sequences exist, a 70%-agreement consensus is used.
All bases used for reconstruction must be Q20 or higher when qualities exist.

Unresolved small pools enter an exact overlap graph. Contained sequences are
collapsed; each edge must have a unique non-repeat-only overlap. Paths start
at the forward primer and must reach the reverse primer. Cycles, multiple
valid products, inconsistent paired placements and nonunique overlap offsets
are unresolved. Reliable library statistics additionally constrain uniquely
placed mate distances. The graph is capped at 256 distinct sequences and 64
path expansions. Hitting either cap defers to inference, never an arbitrary
repeat traversal. No whole-genome assembly or external assembler is used.

## Haploid candidate fallback

Only unresolved loci with rich panel templates receive candidate scoring.
Synthetic products resize the panel repeat interval across a bounded range.
Primer-only loci skip this fallback and remain unresolved if direct observation
and local reconstruction cannot establish a complete product.
The range starts from panel bounds plus padding and observed repeat lower
bounds, expands when probability reaches its edge, and stops at
`--short-max-candidate-repeat-count` (default 100). Half-repeat states are
included where the unit is at least two bases. Direct observations are not
restricted by this inference ceiling.

Each molecule contributes native alignment scores from its reads, a penalty
for violating its observed repeat lower bound, and, for inward flank pairs,
a Student-t fragment-length likelihood. Scores are summed in log space and
normalized across candidates. No diploid model or component EM is involved.
Boundary-only and repeat-rich evidence supply bounds but cannot establish an
exact count. A ceiling-truncated likelihood is not an exact call. Widely
separated, substantially supported fragment clusters are marked `MIXED`
without inventing a single sequence.

Fragment statistics come from uniquely positioned, ordinary inward pairs on
the same non-repeat flank, independently of unknown VNTR length. At least ten
observations are needed. Median/MAD outlier filtering yields mean, SD and 5th/
95th percentiles. A library without enough such pairs contributes no fragment
term. Independently measured `--insert-mean BP --insert-sd BP` may be supplied.

Inferred sequence retains only flank bases with Q20, accepted flank placement,
at least three molecules and 90% base agreement. Uncovered positions and the
unobserved repeat remain `N`. An inferred repeat count is not a measured repeat
haplotype, and should not be used as evidence for intra-repeat SNP absence.

## Shared interpretation and outputs

Assembly, spanning long reads and short-read recovered products pass through
`locus_products.genotype_product`, using `assembly_equivalent_product_allele`
for the same product-length calibration and repeat rounding. Full recovered
sequence preserves intra-locus changes. Canonical product tables also record
`repeat_sequence`, `repeat_motif_cigar` and `repeat_motif_edits`: substitutions
and indels relative to an equal-length motif template. Unknown repeat bases
supply no such calls. Motif coordinates are not unique genomic placements;
nonidentical intervals above 10 kb retain sequence without global alignment.
The existing repeat-masked SNP marker
and reference comparison remain separate from physical repeat measurement.
Reference classification occurs afterward and cannot change the measured
allele. Assembly and long-read output semantics are unchanged.

`calls.tsv` adds `call_method` and `confidence`. Detailed evidence is in
`short_read_repeat_evidence.tsv`, with `n_full_span`, `n_left_boundary`,
`n_right_boundary`, `n_flank_pairs`, `n_repeat_rich`, `n_anchored_repeat`,
`n_softclip_left`, `n_softclip_right`, and `n_uninformative`. It also records
read length, motif length, amplicon/VNTR lengths, depth, the amplicon-to-insert
ratio, likelihood intervals and failure reasons.

Methods are `DIRECT`, `RECONSTRUCTED`, `INFERRED`, `AMBIGUOUS`, `NO_CALL`, or
`MIXED`. Existing PASS/LOW_DEPTH/UNRESOLVED-style status columns remain.
Aggregate likelihood TSVs contain rows only for loci that needed inference.
Per-molecule diagnostics retain evidence and placements rather than a large
molecule-by-candidate likelihood grid. Confidence
is an evidence score, not an externally calibrated error probability.

## Performance and limits

FASTQs are never rescanned per locus. Recruitment uses bounded process batches;
local reconstruction runs after recruitment within the same thread allocation.
Only recruited locus pools are retained, not the complete input FASTQ dataset.
Identical reads share cached anchor/motif work. Graph complexity is deliberately
bounded. Temporary QC files are closed and removed on both success and failure,
unless `--keep-intermediates` is requested. A supplied reference database may
also require a separate combined mapping pass for downstream classification.

Pure repeats longer than the reads often cannot be reconstructed uniquely.
Short flanks may prevent empirical insert calibration. The exact-overlap
assembler is conservative with sequencing errors and very diverse pools.
Recruitment memory grows with on-target depth; the graph cap does not cap the
retained evidence pool. The small 4-mer recruitment filter can lose selectivity
for large panels. Partial SNP calls are conservative and unphased. Shared masked SNP markers
keep their established semantics; the additive repeat-motif edit fields
describe observed repeat variation without asserting unique genome coordinates. Mixture sensitivity and numeric
confidence require validation on real microbial libraries.

See [Illumina commands and validation](../workflows/illumina.md).
