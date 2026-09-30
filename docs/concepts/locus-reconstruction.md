# Targeted paired-end locus reconstruction

The single short-read pathway uses competitive locus recruitment, direct
measurement, targeted reconstruction, haploid likelihood inference, and
read-depth length estimation when a complete locus cannot be recovered.
Assembly and long-read recovery feed the
same `LocusProduct` and `genotype_product` interpreter.

## Reference recruitment and mate retention

SR FASTQ commands require a reference database. Minimap2 maps the QC-filtered
reads against candidate contexts for all panel loci. Per-mate scores are
collapsed across duplicate reference/allele contexts. Locus assignment needs a
margin of both 12 score points and 20%; cross-locus ties remain unassigned.

Supported reference backgrounds supply flank coordinates. Native flank
alignment classifies the original read pairs, including mates without mapping
records. Repeat-only molecules without a usable flank do not identify a locus.
Up to eight equally supported backgrounds are fitted; conflicting lengths remain
unresolved. Observed spans, local reconstruction and fragment likelihoods use
the measurement rules below. Database counts are never copied into missing
sample calls, and unobserved sequence remains `N`.

The lower-level sample-only recruitment and graph-learning helpers remain for
synthetic validation. They are not the user-facing SR FASTQ calling pathway.
Assembly and default accurate long-read recovery need no database. Taxonomic
classification remains a separate database-dependent step in every mode.

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
Contiguous Q20 segments of at least 20 bases remain usable when another part
of their read fails this threshold. A quality-validated complete product is
not rejected because of low-quality bases outside its primers.

Unresolved pools enter a seeded overlap graph. High-depth reads first form
HTSlib coordinate pileups through `pysam`: full-read placements with unique
offsets consolidate substitution errors, with each molecule voting once per
position. Different-length quality-trimmed reads can be placed inside longer
reads with at most 2% substitutions and a unique offset in the chosen container.
If several containers match equally, one holds the observation; the contained
read does not join those containers. Original read coordinates and base votes
are retained, so short reads cannot vote across unobserved extensions.
Indels and ambiguous repeat offsets remain separate. Contained sequences are
collapsed while retaining their molecule support; each edge must have a unique
non-repeat-only overlap, allowing at most 2% substitutions. Exact 15-base seeds
bound the candidate-offset search. Consistent overlap components form a shared
coordinate pileup without enumerating redundant paths through every read start.

Consensus uses the validated overlap coordinates and preserves their observed
length; disagreement below 70% support becomes `N`. Paths start at the forward
primer and must reach the reverse primer. Nonunique overlaps are excluded
without vetoing independently anchored complete paths. Conflicting product
lengths, cycles and inconsistent paired placements remain unresolved. Reliable
library statistics additionally constrain uniquely placed mate distances.
Separate components and alternative paths that agree on product length and
primer ends retain that length, masking disputed interior bases as `N`.
Resource limits apply to 2048 compacted nodes and 64 alternative path expansions,
not the original number of recruited read sequences. Hitting a limit sets
`limit_reached=true` and defers to likelihood inference and read-depth estimation.
No whole-genome assembly or external assembler is used.

## Read-depth length estimates

Unresolved short-read loci also enter a compact weighted de Bruijn graph. Its
edges are observed Q20 21-mers; if no primer-connected graph exists, 15-mers are
tried. Sequence multiplicity increases edge weights instead of creating more
read-fragment nodes. The overlap assembler's 2,048-node and 64-path limits do
not limit this graph, and no paths through repeat cycles are enumerated.
Coordinate-based pileup consolidates substitution errors while retaining
singleton observations. Edges that cannot connect the observed primer endpoints
are removed. There is no minimum edge depth. Motif discovery and a complete
contig are not prerequisites.

One combined replay of the complete, QC-filtered input counts graph k-mers for
all eligible loci, including repeat-only reads omitted during recruitment.
Both strands contribute to the same sequence count. The two unbranched primer
arms define single-copy k-mer depth. Product length is estimated as
`k - 1 + sum(distinct canonical k-mer counts) / single-copy depth`, bounded
below by the shortest primer-connected graph path. This uses read-derived
sequence multiplicity, following the coverage-normalization principle described
by [Kuśmirek and Nowak (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6052550/).
It is a local length estimator, not an implementation of their assembler.

Measured products, identifiable likelihood calls and mixtures retain priority.
Successful calls incur no additional input replay. Depth estimation imposes no
minimum depth, molecule count or read-start diversity. Observed primer-start
position counts are recorded in the diagnostics. Identical primer-start amplicon
reads violate the uniform-coverage assumption and can bias these estimates.
No connected graph or no usable single-copy depth still leaves length
unresolved; abundance alone cannot identify a length without a coverage model
or spanning evidence.

`KMER_DEPTH` values are diagnostic only. Matching sequence outside the locus
can inflate the ratio even when both arms have equal coverage and no other
eligible locus shares a graph edge. Primer connectivity establishes a path,
but cannot assign internal repeat-only reads to that locus. A numeric depth
estimate therefore cannot establish an exact repeat call or replace existing
likelihood evidence. When no length-based point estimate exists, the reporting
layer retains it as `ESTIMATED`, with `confidence_kind=unvalidated_depth` and
confidence 0. Zero means unvalidated, not a probability of being wrong.

Diagnostic lengths are converted to repeat units using shared panel calibration
or sample-learned boundaries. Their intervals describe coverage sensitivity,
with at least 10% length variation, and are not genotype confidence intervals.
They do not modify the caller's likelihood interval. Shared k-mers across
eligible loci are flagged; a zero shared-edge count does not prove uniqueness.

`repeat_length_estimates.json` records the k-mer size, graph size, flank depths,
estimated amplicon length, repeat count, sensitivity intervals and shared-edge
count, with `status=diagnostic_only`, or the reason estimation could not run.
Depth estimates also populate call tables and report bars as provisional
estimates. Their intervals are labeled `coverage_sensitivity`; original
likelihoods remain unchanged in `repeat_likelihoods.tsv`. A depth estimate below
an observed repeat lower bound is rejected.

If no read-derived point estimate exists, `PANEL_PRIOR` uses the panel nominal
count, or the midpoint of its expected range when no nominal count is supplied.
Observed lower bounds constrain it; a computational search ceiling does not
become an observed upper bound. These estimates have `confidence_kind=prior_only`
and confidence 0. They do not create synthetic sequence, raw measurements,
exact-concordance successes, or profile/reference matches. Every detected locus
with a repeat definition thus has an estimate and confidence; missing repeat
metadata remains `PRESENT_COUNT_UNKNOWN`.

Observed complete products and reconstructed contigs then pass through
`run_in_silico_pcr_loci`, `pcr_rows_to_products` and
`legacy_assembly_call_rows`, exactly as in assembly calling. This includes
primer-indel product-size correction, mismatch-round selection and the legacy
allele eligibility rules. A PCR hit for a different locus cannot inherit the
contig's recruited molecules. Inferred sequence is not used as independent
Sassy evidence. The report includes the contig length, PCR-calibrated length,
repeat count and measurement source; raw PCR matches remain available for audit.

## Haploid candidate fallback

Unresolved loci with rich panel templates or successfully learned sample
graphs receive candidate scoring. Graph paths resize the repeat interval
between observed arms across a bounded range. Primer-only loci without both
learned boundaries remain unresolved if reconstruction cannot establish length.
The range starts from panel bounds plus padding and observed repeat lower
bounds, expands when probability reaches its edge, and stops at
`--short-max-candidate-repeat-count` (default 100). Half-repeat states are
included where the unit is at least two bases. Direct observations are not
restricted by this inference ceiling.

Native alignments establish locus and boundary coordinates once. A boundary
observation excludes shorter lengths and assigns equal likelihood to every
compatible longer length. Repeated observations of the same bound cannot
sharpen its unobserved end. Inward flank pairs contribute Student-t fragment
likelihoods, summed in log space and normalized over half-repeat length states.
Different learned flank phases are converted to the same physical product
coordinates before comparing evidence. Repeated fragment geometries share
vector computation while retaining their molecule multiplicities. No candidate
sequences, per-candidate alignments, diploid model or component EM are required.
Graph loop counts and uncertainty intervals are converted to the same
product-length calibration used for the reported MLVA allele; motif phase
choices do not change the reported allele coordinates.
Boundary-only and repeat-rich evidence supply bounds but cannot establish an
exact count. A ceiling-truncated likelihood is not an exact call. Widely
separated, substantially supported fragment clusters are marked `MIXED`
without inventing a single sequence.

Fragment statistics come from uniquely positioned, ordinary inward pairs on
the same non-repeat flank, independently of unknown VNTR length. Any observed
pair can contribute. Median/MAD outlier filtering yields mean, SD and 5th/
95th percentiles. The automatic predictive SD includes uncertainty in the mean
and shrinks variance toward a broad prior (10% coefficient of variation, at
least 10 bp, with weight equivalent to two observations). This is a modeling
assumption, not empirically calibrated confidence or a depth cutoff. A singleton
still contributes but cannot imply a 1-bp library SD. Explicit library overrides
retain their supplied mean and SD.
A library without such pairs contributes no fragment
term. Independently measured `--insert-mean BP --insert-sd BP` may be supplied.

Inferred sequence retains only flank bases with Q20, accepted flank placement,
90% base agreement, without a minimum molecule count. Uncovered positions and the
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
`MIXED`. A single observation can be called; no depth cutoff emits `LOW_DEPTH`.
Legacy status values remain readable in existing outputs.
Aggregate likelihood TSVs contain rows only for loci that needed inference.
Per-molecule diagnostics retain evidence and placements rather than a large
molecule-by-candidate likelihood grid. Confidence
is an evidence score, not an externally calibrated error probability.

## Performance and limits

FASTQs are never rescanned per locus. Reference recruitment uses a combined
mapping pass and replay; unresolved depth estimation can add a combined replay.
Recruitment uses bounded process batches;
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

The synthetic internal-read rescue benchmark can be run with
`python scripts/benchmark_reconstruction.py --scenario rescue --molecules 10000 --threads 1 --output /tmp/mlvamaps-rescue`.
It tests primer-only panels with missing internal reads and off-target
background; it does not establish accuracy or throughput on real libraries.
