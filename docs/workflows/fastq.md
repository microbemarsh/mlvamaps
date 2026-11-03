# FASTQ and amplicon sequencing workflow

FASTQ calling requires a completed `--database`. Short and long reads map
directly to its observed reference loci and share molecule-level repeat
inference and reference classification. Every valid configured locus receives a
numerical prediction when an eligible prior exists. Direct measurements take
precedence; uncertainty changes confidence rather than suppressing the count.

```bash
mlvamaps call --database reference_build -i sample.fastq.gz -o results --sample-mode isolate
```

Use `--sample-mode isolate` for the isolate-versus-assembly comparison. The
default `metagenome` mode fits mixtures and should be used for mixed samples.

For paired reads, see the [Illumina workflow](illumina.md).
For current uncertainty rules, limits, and same-input old/new comparisons, see
the [mandatory-count validation report](fastq-complete-count-validation.md). The historical
runtime tables below describe earlier implementations.

## 1. Load the database panel and filter reads

The saved `reference_panel.tsv` supplies primers, locus names, repeat
calibration, and fingerprint order unless `-p` supplies an explicit panel.
The database's per-locus FASTAs supply the observed mapping targets.
Current-schema databases remain usable; their older `competitive_mapping/`
assets are retained for compatibility but do not supply the primary FASTQ
mapping targets or numeric allele estimates.

An optional Deacon target screen can precede QC. In metagenome mode, an available
database recruitment index is used automatically. Reads then pass the configured
length and mean-quality filters; original qualities remain attached.
The outputs include `qc_summary.tsv` and `filtered_reads.fastq.gz`.

## 2. Recover mapped locus segments and rescue unresolved reads

The boundary-refinement approach is informed by
[mahi-mahi](https://github.com/nanoporetech/mahi-mahi), implemented with the existing
mapping, Sassy and SPOARS components. Mahi-mahi is not a runtime dependency;
MLVA product selection and repeat calibration remain assembly-compatible.

Mapping runs first. Both read workflows use CIGAR coordinates to constrain local
flank refinement and measure the original sequence in its mapped orientation.
Available reference flanks are searched in bounded windows, trying up to 75, 35,
and 21 bases with 24-base padding. Reference sequence locates the interval but
never supplies missing sample bases or determines its repeat count.

Long reads without a measured complete product enter the Sassy whole-read PCR
rescue step. Reads long enough to contain another product also enter rescue.
Disjoint locus intervals on one genomic read contribute once to each locus,
including intervals on opposite strands. Overlapping competing loci remain
ambiguous. Primer orientation, pairing, IUPAC handling, and product-length rules
remain unchanged. Tied loci or product sizes do not become accepted assignments.
Rescue outputs include `primer_rescue.fasta.gz` and `in_silico_pcr/` when needed.

## 3. Map to observed reference loci

Minimap2 maps retained reads competitively against observed database amplicons.
Identical sequences within a locus share a target while retaining all source
reference IDs. Alternative reference mappings are preserved. Synthetic repeat
expansions do not define the sample's allele hypotheses in this path.

For short reads, flank recruitment and mate rescue run during QC and feed the
same measurement stage. Surviving orphan mates are mapped as unpaired reads;
they retain the original molecule identity. Conflicting mate measurements remain
alternatives from one molecule, and a complete product takes precedence over a
partial measurement. The appropriate technology preset is used for mapping.
The observed target FASTA, metadata, and provenance are recorded under
`candidate_mapping/`; minimap2 builds its target index in memory.

SAM is streamed through htslib to temporary BAM. The applicable intermediate or
debug flag retains `candidate_mapping/candidate_alignments.bam` for short reads.
Long reads retain this BAM for recruitment, replacing the separate recruitment
mapping. Decoded alignments remain available to reference classification.

## 4. Measure molecules and call repeats

Complete primer products use the same calibrated product-length functions as
assemblies. Raw lengths and measurement uncertainty are retained through
inference; assembly rounding is applied to the final allele. Isolate mode fits
a single-allele likelihood with a broad error component for occasional outliers.
Metagenome mode uses the existing weighted EM solver to estimate per-locus
allele fractions. A secondary allele requires sufficient fraction and effective
molecule support, convergence, and a BIC-penalized likelihood improvement over
the single-allele model. EM convergence
and fractions are reported separately from single-allele model support.

Identical observations share likelihood calculations while retaining molecule
multiplicity. Long-read anchor measurements are cached by locus, sequence and
quality, including failed measurements, rather than repeated for each reference
allele. These caches are bounded and local to one inference run.

Independently resolved repeat flanks can also measure
a repeat interval. For historical and minimal panels, a read missing a primer
uses the reference products' nonrepeat length to convert that observed interval
into assembly allele units. Configured flanks may be short boundary anchors;
their lengths are not assumed to cover the whole nonrepeat interior. The
conversion requires unambiguous reference anchors and agreement across reference
nonrepeat lengths. Otherwise the flank span does not supply a direct calibrated
allele, and the existing reconstruction/estimate paths remain available. Complete
sample products retain their own measured calibration. No reference bases are
exported as sample sequence, and this conversion cannot establish an unobserved
sample-specific indel outside the repeat.

Observed reference labels use this same product calibration, including correction
for primer indels. Minimal panels retain the assembly convention for nonrepeat
spacer sequence outside the configured flank anchors.

When a read endpoint truncates a configured flank longer than 21 bases, the
boundary can be recovered from its terminal 21 bases if reference nonrepeat
calibration is unique. Both boundary anchors must be ordered, unambiguous,
contain at least three base types, and be at least 90% identical; ambiguous or weak partial anchors
remain unresolved. Strong anchors can support novel repeat lengths outside the
panel's expected range. This does not resolve arrays longer than the available
fragment span.

The regression `test_calibrated_flank_reads_preserve_assembly_distance_matrices`
checks an offset panel on matching assembly and short-read samples, including
fractional alleles, and requires identical repeat-only and combined-marker
distance matrices. This is a regression check, not a real-cohort accuracy claim.

Paired reads can contribute a uniquely merged fragment or a
supported local assembly, as described in the Illumina workflow.
For concrete DNA anchors in windows of at most 512 bases, an in-process edit-distance
existence check rejects impossible matches before launching Sassy. Accepted
approximate matches still use Sassy's coordinates and traceback.

Measured evidence supplies direct alleles. If direct sizing remains unresolved,
repeat-spanning CIGAR projections with meaningful aligned support at both flanks
can provide `ESTIMATED` counts. Conflicting measurements and partial reference
matches retain best estimates with candidate distributions. Reference copies
and mate records do not multiply molecule support. Estimates are kept separate
from measured counts used in reference classification.
Novel and fractional measured alleles remain callable outside the catalog range.
User-supplied or conservatively flank-calibrated paired-fragment lengths provide another
labelled estimate when no direct measurement is available; see the
[Illumina guide](illumina.md).

`common_locus_calls.tsv` distinguishes `presence_support` from informative
`molecule_support`. `primary_molecule_support` and `secondary_molecule_support`
count independent observations matching the corresponding rounded allele;
neither is inferred by multiplying a fitted fraction by total depth. These can
be smaller than total informative support, especially for noisy sizing or a
consensus-derived estimate. `molecule_candidate_evidence.tsv` retains reference IDs,
alignment diagnostics, evidence types, and inferred molecule counts.

`genotype_evidence.tsv` distinguishes exact, estimated, bounded, unresolved,
ambiguous-locus, and possible-mixture evidence. Its interval is an observation
envelope, not a calibrated confidence interval. `sequence_alleles.tsv` reports
identical complete products observed in at least three independent molecules;
it does not manufacture a polished sequence from uncertain bases. Isolate
bimodality retains the best primary count with possible secondary
evidence. Neither that label nor legacy mixture support confirms strain mixtures.

Low-depth calls, mixtures, ambiguous measurements, detected-but-unresolved loci,
and missing loci retain distinct statuses. The finalizer uses a locus-specific
prior when sample repeat evidence is absent, including when every read fails QC.
These predictions carry `prior_only=True` and zero data-evidence confidence.
Missing evidence is never converted to an arbitrary zero. A locus with neither
informative evidence nor an eligible prior fails with a configuration error.
See [prior and confidence definitions](../reference/outputs.md#complete-fastq-predictions-and-confidence).

## 5. Confirm representative sequences

A dominant observed product takes the fast path through assembly PCR. Otherwise,
supported mapped segments are clustered by measured allele, or by observed length
when unmeasured, and corrected with the installed SPOARS POA. Correction uses at
most 100 quality-ranked independent reads. It can recover primer sites before
requiring successful PCR on individual reads. Short reads also retain their
bounded fragment graph for products not spanned by an individual read or pair.

A supported consensus supplies an assembly-PCR measurement without adding votes
or replacing original measurements used for reference classification. In isolate
mode, nearby uncertain raw measurements can contribute to the same reconstruction
even when they straddle a rounding boundary. A consensus evaluated against molecule
likelihoods retains that likelihood and its confidence; a rejected consensus cannot
overwrite the numerical prediction. A conflict with a confident molecule call
produces `AMBIGUOUS` and `consensus_status=disagrees`. Mixtures remain separate.
Products recovered without a measurable molecule likelihood use
`probability_kind=reconstruction_support_fraction`; this is observed support,
not calibrated accuracy. Only products agreeing with final resolved calls enter
`query_loci.fasta`; `query_locus_reconstruction.tsv` records the decision.

Recruitment reuses the observed-reference BAM. Optional representative SNP
mapping and the historical evidence tables remain available.

Relevant outputs include `local_locus_products.fasta.gz`,
`local_assembly_concordance.tsv`, mapped-variant tables, and locus SNP evidence.

## 6. Classify against references using the shared evidence

Reference classification reuses the observed-reference alignments from step 3;
it does not map the FASTQ a second time. Sequence mismatches and measured repeat
length differences both inform reference likelihoods. A measured paired-fragment
count contributes once per molecule, even when both mates align. Where no measured
count exists, the existing alignment-derived repeat evidence remains available
to classification. A reported estimate does not become an independent measured
count in that classification model.

For long genomic reads, classification also retains independently assigned,
non-overlapping locus intervals. Repeated mappings do not add molecule support,
and overlapping locus ties remain unresolved.

Equivalent references remain grouped. Isolate support and metagenome mixture EM
retain an unclassified component, and weak or ambiguous evidence remains explicit.
Model support is conditional on the available reference catalog.

Outputs include `classification/molecule_reference_likelihoods.tsv.gz`, reference
matches, and taxonomic identification tables when applicable. Primer-only recovery
without an observed-reference alignment can support a locus count but does not
by itself establish a reference identity.

## 7. Export results

Start with `report.html`, `calls.tsv`, `mlva_fingerprint.tsv`, and
`locus_repeat_counts.tsv`. Fingerprints contain all predictions; biological
comparisons exclude prior-only values. Use the same panel calibration and `--assembly-round-tolerance` when
comparing FASTQ and assembly fingerprints.

Low-confidence predictions supported by sample evidence remain available to
profile and combined-marker distances. Their status, uncertainty, and support
remain explicit; including them does not make the estimate exact. For FASTQ
combined-marker exports, the current reconstruction audit governs sample sequence
eligibility. Rejected or mismatched products cannot re-enter through stale
classification or local-assembly files. Empty runs clear the current query
sequences and evidence table. A disabled distance component contributes neither
distance nor callable overlap.

Synthetic regression coverage includes single-end, paired, overlapping,
shotgun, and long reads, novel and fractional alleles, zero-repeat products,
primer indels, and primer rescue when mapping misses reads. Performance and
accuracy on a particular sample still require comparison with its known truth.

## FASTQ concordance sweep (2026-11-03)

The complete native-enabled suite passed **552 tests**. Regressions cover orphan
mates, conflicting measurements, independent loci on both strands of a genomic
read, primer indels, partial flanks, consensus disagreements, empty-run artifacts,
and low-confidence count/depth retention through exports. Matching assembly and
FASTQ samples produce identical repeat-only and combined-marker distance matrices
in both historical and minimal-panel calibration fixtures.

A fresh-process comparison used identical retained inputs and the working source
saved immediately before this sweep (including the preceding flank-calibration
fix). One trial per layout used two copies and isolate mode:

| Fixture | Before exact | After exact |
| --- | ---: | ---: |
| Clean paired short reads | 11/11 | 11/11 |
| Clean shotgun short reads | 11/11 | 11/11 |
| Clean long-read products | 11/11 | 11/11 |
| Difficult paired short reads, all known-truth loci | 6/14 | 9/14 |

Partial-flank recovery corrected the homopolymer call from 15 to 16, increasing
informative support from 6 to 18 molecules without changing total presence (24).
The noisy locus changed from a lower-bound estimate of 12.5 to a measured estimate
of 16 with three informative molecules out of 40 detected. It remains
`AMBIGUOUS`; its low confidence does not erase the count. The long-motif lower-bound
estimate now uses assembly rounding (raw 7.8889, reported 8).

All 15 challenge loci retain numerical predictions, including the locus with no
truth allele. Five known-truth discrepancies remain: expansions of 40 and 400
repeats reported as lower-bound estimates of 17.5, two indistinguishable loci
with truth 5 and prior predictions 4, and a zero-coverage locus with truth 6 and
prior prediction 4. These remain in the accuracy denominator. The fragment span
does not resolve either expansion; no calibrated insert distribution was supplied.
The noisy fixture's Q40 qualities do not describe its simulated 2–4% error rate,
so these results do not calibrate real-data confidence.

The [validation record](../reference/fastq-concordance-sweep.json) retains source
and input hashes, environment versions, timings, support metrics, and every
disagreement. Single-trial timings are descriptive, not speedup evidence. No noisy
simulated ONT challenge was run. The inaccessible server cohort and Parsnp
matrices were not evaluated, so these checks do not establish 100% real-cohort
agreement.

## Synthetic runtime check

On 2026-10-19, three alternating runs of the pre-optimization and updated workflows
used two threads and fresh processes on the same machine. Median calling times,
including output generation but excluding Python imports, were:

| Input | Before | After | Time reduction |
| --- | ---: | ---: | ---: |
| 11,100 pairs of 100-base reads | 7.60 s | 5.05 s | 34% |
| 132 long-read amplicons | 2.69 s | 1.16 s | 57% |

Both inputs used an 11-locus synthetic panel. The paired input contained 1,100
locus-associated pairs and 10,000 random background pairs; the long-read input
contained 12 copies of each product. The baseline already used observed-reference
mapping. All 28 short-read and 37 long-read output tables matched across versions
in each trial after normalizing output-directory paths and decompressing gzip
tables. These small synthetic timings measure reduced Sassy process overhead;
they do not establish throughput or accuracy on real server datasets.

## Raw-length inference runtime check

On 2026-10-20, the previous caller and the raw-length/EM implementation were
compared on the same machine with two threads and five runs per workflow.
The benchmark times FASTQ calling and output generation, excluding database
construction and assembly truth generation.

| Synthetic input | Previous median | Updated median | Time reduction |
| --- | ---: | ---: | ---: |
| 4,400 full-product pairs | 0.723 s | 0.722 s | Essentially unchanged |
| 29,000 shotgun pairs | 3.751 s | 3.507 s | 6.5% |
| 4,400 long-read products | 14.643 s | 10.196 s | 30.4% |

All resulting fingerprints matched the assembly truth. Inputs use an 11-locus
panel covering novel, fractional, zero-repeat and primer-indel cases, with
replicated high-quality reads. These timings establish performance for this
synthetic isolate workload, not real noisy reads or mixtures.

Reproduce with dev dependencies and minimap2/Sassy on PATH:

```bash
python scripts/benchmark_fastq_inference.py --repetitions 5 --copies 100
```

The [benchmark record](../reference/fastq-likelihood-benchmark.json) contains
individual timings, the baseline commit, environment and input conditions.
Regression tests additionally cover noisy raw-length aggregation, allele
mixtures, duplicate-reference invariance, uncertain consensus and calibrated
paired-fragment estimates.

## Mapping-guided recovery runtime check

The shared boundary-refinement update was compared with the preceding commit
using three alternating serial trials, each in a fresh interpreter, with 100
copies of the 11-locus concordance fixture, two threads and isolate mode.
Database construction and assembly truth generation were excluded from timing.

| Workflow | Previous median | Updated median | Reduction |
| --- | ---: | ---: | ---: |
| Paired short reads | 0.934 s | 0.836 s | 10.5% |
| Shotgun short reads | 4.503 s | 4.409 s | 2.1% |
| Long reads | 10.275 s | 8.925 s | 13.1% |

All fingerprints matched assembly truth. Short-read differences are small and
vary between runs; these replicated high-quality fixtures do not establish
real-data accuracy or production speedups. End-to-end regression tests separately
cover primer-error consensus recovery in both read workflows, novel repeat
lengths, reverse orientation, and conflicting evidence.

The [benchmark record](../reference/mapping-guided-benchmark.json) includes the
baseline commit and individual timings. Use the existing benchmark command above
on each checkout in fresh interpreters to repeat the comparison.
