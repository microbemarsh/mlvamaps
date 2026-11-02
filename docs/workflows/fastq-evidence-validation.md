# FASTQ evidence changes and validation, 2026-11-02

The historical audit below describes commit `6851a68` and its predecessor.
Its numerical-suppression behavior is superseded by the
[mandatory-count implementation and benchmark](fastq-complete-count-validation.md).
Current FASTQ output always retains a best numerical estimate for valid configured
loci with eligible priors. Confidence still reflects uncertainty; prior-only
counts have zero sample-evidence confidence and are excluded from biological
comparisons. The earlier measurements below remain unchanged as a historical record.

This change hardens the existing caller against duplicate evidence, ambiguous
recruitment, unsupported reference-count estimates, and excessive computation.
It adds observed sequence alleles and explicit uncertainty outputs. It does
**not establish publication-level accuracy**. The synthetic comparison finds
fewer incorrect short-read calls at a lower call rate, one remaining incorrect
short-read estimate, and lower long-read call rate without an accuracy gain.

The [benchmark record](../reference/fastq-evidence-benchmark.json) contains
source/input hashes, individual trials, accuracy/call-rate curves, coverage,
motif-length, array-length and error-rate strata, ablations, and profiles.
The [disagreement record](../reference/fastq-evidence-disagreements.json) retains
both versions' rows, known synthetic truth, wrong calls, and missing calls.
Neither agreement with the baseline nor agreement with an assembly alone is
treated as truth.

## Original architecture and implementation scope

The audited baseline is commit `d468d380c4adae3f63d4243b0f6f63ff4495327a`.
The repository is a Python package built with setuptools, with pytest tests and
the `mlvamaps` CLI. Existing dependencies are NumPy, Parasail, pysam/htslib,
regex, and SPOARS; minimap2 and Sassy provide native mapping and anchor search.
No language, mandatory dependency, external service, or separate genotyper was
added. The removed GangSTR-inspired path was not reintroduced.

Assembly input already used Sassy PCR and calibrated product sizing. FASTQ paths
already shared `CandidateAlignment`, `CandidateEvidence`, competitive mapping to
observed database amplicons, product/flank measurement, raw-length likelihoods,
and reference classification. Short reads already had flank seeds, mate rescue,
unique overlap merging and bounded local assembly. Long reads already had local
flank searches, primer rescue and SPOARS consensus. The implementation extends
these components rather than replacing them.

Panels represent primers, oriented flanks, motif/unit length, expected range,
and product-size calibration. Observed reference contexts preserve reference
and taxon membership when identical targets are collapsed. Reference typing
reuses the competitive alignments and original measured molecule counts.
Existing repeat/SNP combined-marker export remains in place; it masks repeats
in the SNP component and does not ingest the new sequence table as extra SNPs.

## Staged gates

Each stage's relevant tests passed before the next stage started. Later review
findings returned to implementation and testing before accepting a benchmark.

| Stage | Gate recorded during development |
| --- | --- |
| 1. Audit and baseline | Full native-tool baseline: 389 passed; comparison-harness subset: 26 passed. |
| 2. Evidence correctness | 120 relevant tests passed. |
| 3. Long-read safeguards | 90 relevant tests passed, followed by a 62-test refinement gate. |
| 4. Short-read safeguards | 123 relevant tests passed. |
| 5. Reporting and compatibility | 219 relevant tests passed. |
| 6. Profiling, differential review and fixes | Final full suite: 422 passed, no skips; all 75 final comparison workers succeeded. |
| 7. Documentation and package review | Wheel build succeeded; reproduction commands and output semantics documented. |

The final full suite includes assembly, short/long FASTQ, metagenome inference,
classification, SNP/combined-marker export, reporting and malformed-input tests.
New checks cover omitted secondary sequence, reverse/hard-clipped supplementary
alignments, cross-locus ties, duplicate/conflicting fragments, enormous CIGAR
deletions, likelihood underflow, resource refusal, multimodal isolates, insert
calibration, partial bounds, sequence identity, atomic TSV failure, and benchmark
missingness. Repeated native benchmark runs check deterministic artifact hashes.

## Algorithms and source map

| Files and key functions | Change and purpose |
| --- | --- |
| `alignment_evidence.py`: `CandidateAlignment`, `CandidateEvidence`, `evidence_row`; `models.py`: `AnchorMeasurement` | Compact slotted evidence with molecule identity, strand, supplementary status, local-query offset, candidate-locus ambiguity, boundary flags and refinement provenance. Unknown repeat counts remain missing. |
| `minimap_mapping.py`: `parse_candidate_alignments`; `mapped_reconstruction.py`: `recruit_mapping_fragments`, `recruitment_identity` | Preserve secondary/supplementary evidence, restore missing sequence/quality and hard-clip coordinates, retain cross-locus ambiguity, and tolerate repeat indels left-aligned just outside the nominal boundary. |
| `locus_measurement.py`: `reference_interval_to_query_interval`, `find_anchor`, `measure_locus_product`; `calling.py`: `gaussian_allele_probabilities` | Project CIGAR intervals in O(operations), validate slices/qualities, retain tied boundaries, use local allele grids and stable likelihood normalization, and handle a true zero-repeat allele without assigning it empty-slice quality. |
| `long_read_evidence.py`: `flank_error_model`, `refine_repeat_sequence`, `long_read_measurement_metadata`, `extract_long_read_evidence` | Bounded native motif scoring, shared regularized ONT error rates, explicit boundary/limit uncertainty, and the same uncertainty rules for mapped and primer-rescued products. |
| `short_read_evidence.py`: `estimate_insert_distribution`, `partial_repeat_lower_bound`, `extract_short_read_evidence` | Independent fragment geometry, conservative automatic insert calibration, exact anchored lower bounds, meaningful two-flank support for CIGAR spans, and preserved ambiguity flags. |
| `allele_inference.py`: `InferenceThresholds`, `infer_alleles`, `estimate_unresolved_calls` | Bounded likelihood allocation/EM, one vote per molecule, one derived-product likelihood, haploid bimodality detection, observation envelopes, raw score gaps, and refusal to fill unresolved calls from partial reference matches. |
| `evidence_reporting.py`: `write_evidence_reports` | Explicit uncertainty statuses and independently repeated observed sequence alleles. |
| `io.py`: `ReadSpool`, `read_fastq`, `read_fastq_pairs`, `write_tsv`; `qc.py`: `filter_reads` | Stream input to disk, reject malformed/duplicate input identities, and atomically replace individual TSV files. |
| `pipeline.py`, `unified_fastq.py`, `recruitment.py`, `locus_assignment.py` | Connect the disk spool, local mapped sequences, rescue metadata, supplemental outputs and limits through the existing pipeline. |
| `local_assembly.py`: `sequence_consensus`; `mapped_variants.py`; `mixture.py` | Apply resource bounds to downstream consensus and diagnostic variant/mixture alignment as well as primary inference. |
| `alignment.py`; `sassy_cli.py`; `locus_measurement.py` | Reuse native implementations, bypass identical-sequence alignment, bound the shared search cache, and batch plausible small-product anchor searches. |
| `cli.py`, `short_reads.py`, `short_read_mapping.py` | Pass four advanced limit options through both FASTQ entry points. |
| `scripts/benchmark_fastq_inference.py`, `scripts/fastq_challenge_fixture.py` | Same-input, fresh-process comparison; independent synthetic truth; retained failures and ablations. |

### Fragment and long-read interpretation

Molecule IDs join mates, alternative references and supplementary records.
Multiple rows never become independent biological observations. Unresolved
cross-locus assignments remain evidence but do not vote for an allele.
Reconstructed products carry their supporter IDs, but one reconstructed product
does not become many independent likelihood observations or classification votes.

Long-read CIGAR coordinates locate bounded flank-search windows; the observed
sequence supplies product/repeat length. Full-span alignment estimates require
aligned support from both flanks, normally 21 bases where available. An isolated
base of the opposite flank no longer establishes a complete span. Ambiguous
flank boundaries increase measurement uncertainty and cannot qualify as exact.

Canonical and imperfect products retain the existing repeat-count convention.
Poor motif fits can invoke Parasail banded scoring around the observed length:
a radius of two repeat units expands to four if needed, with explicit limits.
These scores are diagnostic and widen uncertainty; they do not silently rewrite
the observed count into the best motif fit. Compound-repeat polishing and a
validated indel-correction model are not implemented by this change.

For ONT/ONT-HQ, primary high-identity flank alignments provide pooled substitution,
insertion and deletion rates. Repeat indels do not train those rates. Fitting
requires at least 20 molecules and 2,000 flank bases and uses a 2,000-base prior;
otherwise conservative prior rates are retained. This is a shared uncertainty
proxy, not a basecaller- or sequence-context-calibrated model. Rescue products
receive the same adjustment as mapped products.

Insert calibration requires at least 100 independent, high-quality inward pairs
wholly within nonrepeat flanks, robust outlier trimming, nonzero dispersion and
distance from target edges. It commonly remains unavailable for short reference
amplicons. User-supplied mean/SD remain supported. Opposite-flank pair estimates
remain `ESTIMATED`: mapping/recruitment selection bias is **not calibrated**.
Repeat-internal depth is presence evidence, not a standalone length estimator.
The partial-read fallback currently supplies exact left-anchor terminal motif
lower bounds, including reverse-oriented reads; it is not a general partial-read
joint likelihood model.

Isolate mode keeps a haploid primary model. Two sufficiently separated modes
with independent support yield possible secondary-allele evidence and can
withhold a primary count. Existing bounded metagenome EM remains available.
Neither path establishes that a secondary mode is a distinct strain.

## Confidence, sequences and compatibility

`genotype_evidence.tsv` and the appended `genotype_status` field in
`common_locus_calls.tsv` are the evidence-aware interpretation. The full
[status definitions](../reference/outputs.md#fastq-call-statuses) distinguish
`EXACT`, `ESTIMATED`, `LOWER_BOUND`, `UNRESOLVED`, `AMBIGUOUS_LOCUS` and
`POSSIBLE_MIXTURE`. A count alone, even legacy `PASS`, is not exact evidence.
Intervals are observation envelopes; likelihood differences and legacy fields
named probabilities are uncalibrated model support.

`sequence_alleles.tsv` requires at least three independent identical complete
products, resolved flanks and concrete DNA bases. It reports the observed
product hash, repeat sequence, raw/conventional count, motif pattern,
interruptions and terminal partial bases. Same-count sequence differences remain
distinct. This conservative output cannot recover a sequence when independent
reads disagree; it does not label an uncertain polished consensus as observed.

Conventional fingerprints, per-locus counts, QC, SNP and combined-marker formats
remain available. `common_locus_calls.tsv` and detailed evidence append fields.
The intentional behavioral change is missing numeric counts where conflicting
or partial reference-conditioned evidence previously supplied an estimate.
Consumers that require exact evidence should join by sample/locus to the new
status table, rather than treating legacy `PASS` or a nonempty cell as exact.
Missing alleles remain empty, never zero. Genuine zero-repeat measurements remain
supported.

The effective insert calibration is in `library_insert_distribution.tsv`.
Legacy short-read QC/run-metadata insert fields still describe user-supplied
settings only. Per-file TSV writes are atomic; the full output directory is not
transactional. A failed run must not be treated as a completed result directory.

## Computation and memory limits

| Resource | Default | Behavior at the limit |
| --- | ---: | --- |
| Candidate alleles | 2,048 | Preserve evidence/envelope; no forced point estimate. |
| Unique observations × candidates | 262,144 cells | Refuse dense likelihood allocation and report fallback. |
| Expensive refinement/consensus reads | 100 per locus | Deterministic bounded selection; other observations remain evidence. |
| Expensive alignment sequence | 10,000 bases | Skip refinement with a diagnostic; cheap direct measurement can continue. |
| Estimated alignment work | 2,000,000 cells | Bound banded scoring, differing-sequence alignment and POA work before entry. |
| Refinement iterations | 200 | Retain convergence/fallback diagnostics. |

The four CLI options expose candidate count, likelihood cells, refinement reads
and alignment bases. These are computational defaults, not validated biological
maximum repeat lengths. They do not exclude longer directly measurable arrays.
The existing short-read graph retains its separate node/path bounds.

Input FASTQ is streamed; long-read sequence/quality live in a temporary disk spool
with ID/offset indexing. Long-read alignment records retain local locus windows.
Failed primer rescue releases original sequence. Retained mapped evidence,
accepted rescue records and read-ID indexes still grow with input/target depth;
the complete pipeline is not constant-memory. Entry-count-bounded caches are
not byte-budgeted. Large real-world panels and depth need additional scaling tests.

## Measured benchmark results

Environment: macOS arm64, Python 3.13.15, minimap2 2.31-r1302, Sassy 0.2.1;
dependency versions are in the JSON record. Each version runs three times in
alternating serial fresh processes with the same environment and exact input
bytes. Reported wall/CPU exclude imports and fixture preparation; CPU includes
native children. Peak RSS is the larger of parent and largest child, **not their
concurrent sum**. No tests ran concurrently with final timings.

The clean fixture has 11 explicitly specified alleles, including zero, novel and
fractional counts. Each layout retained all 11 correct calls in both versions.
The challenge fixture has 15 loci: 14 known present loci (including one with zero
coverage) and one absent locus. It includes a 400-unit array, imperfect and
partial repeats, homopolymer errors, similar loci and a 5/9 mixture. A mixture's
numeric count is scored against either true allele, which is lenient toward the
old single-count result. No difficult loci were excluded.

| Workflow | Input fragments/reads | Wall seconds, old → new | CPU seconds, old → new | Peak MiB, old → new |
| --- | ---: | ---: | ---: | ---: |
| Clean paired | 88 | 0.136 → 0.149 | 0.086 → 0.094 | 77.2 → 76.6 |
| Clean shotgun paired | 580 | 0.837 → 1.005 | 0.547 → 0.629 | 77.3 → 77.1 |
| Clean long | 88 | 0.916 → 0.697 | 0.720 → 0.577 | 80.6 → 80.3 |
| Challenge paired | 314 | 0.687 → 0.628 | 0.513 → 0.494 | 83.0 → 83.1 |
| Challenge long | 314 | 5.756 → 4.978 | 4.921 → 4.609 | 273.3 → 183.5 |

All entries are medians of three trials. Throughput, stage timings and tier
fractions for each trial are retained in the benchmark record. On updated
challenge long reads, median mapping was 1.574 s, recruitment/extraction 2.459 s,
and inference/reconstruction 0.018 s. These internal stages exclude later
classification and report generation. The clean paired and long cases used
Tier 1 for every locus; the clean shotgun case used Tier 3 for 8/11 loci.
Tier fractions report the highest tier, so Tier 2 work can be included in Tier 3.

| Challenge result | Old | New |
| --- | ---: | ---: |
| Paired: numeric calls / present loci | 11/14 (78.6%) | 6/14 (42.9%) |
| Paired: correct / numeric calls | 6/11 (54.5%) | 5/6 (83.3%) |
| Paired: mean absolute repeat error among calls | 40.818 | 0.167 |
| Long: numeric calls / present loci | 11/14 (78.6%) | 9/14 (64.3%) |
| Long: correct / numeric calls | 11/11 (100%) | 9/9 (100%) |
| Long: mean absolute repeat error among calls | 0 | 0 |

Both versions had zero numeric calls for the one absent locus. Unique assignment
accuracy was 100% among assigned molecules in both versions, with assignment
rates 79.6% paired and 86.9% long. These small artificial panels do not establish
real-species specificity. All five reported paired sequence alleles and all
eight long-read sequence alleles (seven loci) matched synthetic products.
Missing sequence alleles are not scored as correct.

The updated short-read homopolymer count is **15 for truth 16**. It is
`ESTIMATED` in the evidence table despite legacy `PASS`; this remains an error.
The old 400-unit short-read array estimate of 4 is now a lower bound, not a
successful 400-unit call. In long-read mode, the new uncertainty model withholds
the correctly called 400-unit array and the isolate mixture's primary count.
Thus the measured short-read gain is mainly refusal of unsupported counts;
there is no measured long-read accuracy gain here.

At model support ≥0.95 the updated paired challenge has 5/14 calls, all correct.
This is a descriptive accuracy/call-rate point on the same synthetic fixture,
not an independently calibrated threshold recommendation.

### Ablations and performance investigation

Disabling the ONT error model restores the correct long-read 400-unit numeric
call (10/14 calls), without reducing accuracy on this fixture. Disabling motif
refinement changes no numeric challenge calls. Disabling partial bounds removes
lower-bound reporting for four difficult paired loci without changing numeric
accuracy or call rate. The three ablations do not isolate every coupled
recruitment/correctness change and do not demonstrate a causal accuracy benefit
for each new mechanism.

All checked TSV artifacts were byte-identical across repetitions, including
calls, common calls, evidence, sequence alleles, fingerprints and SNP tables
where produced. Timing fields and gzip headers are not part of this check.

The clean shotgun case regressed by about 20% (0.168 s). Separate profiles on
the same input found additional FASTQ validation (0.025 vs 0.002 s), sequence
evidence reporting (0.042 s), TSV writing (0.048 vs 0.034 s), and native search
time (0.647 vs 0.589 s). Profile timings include instrumentation overhead and are
not the performance table. The regression is retained and unresolved as a speed
tradeoff; no general short-read speedup is claimed. Long-read profiling still
identifies alignment parsing and native anchor searches as the main costs;
motif refinement consumed about 0.001 s in that profiled challenge run.

## Reproduce and inspect

Use a clean baseline checkout/archive and the same environment for both sources:

```bash
git archive d468d380c4adae3f63d4243b0f6f63ff4495327a -o /tmp/mlvamaps-baseline.tar
mkdir /tmp/mlvamaps-baseline
tar -xf /tmp/mlvamaps-baseline.tar -C /tmp/mlvamaps-baseline
python -m pytest -q
python scripts/benchmark_fastq_inference.py \
  --baseline-source /tmp/mlvamaps-baseline \
  --workdir /tmp/mlvamaps-comparison \
  --repetitions 3 --copies 2 --challenge \
  --ablations no_error_model no_motif_refinement no_partial_bounds
python -m pip wheel . --no-deps --no-build-isolation -w /tmp/mlvamaps-wheels
```

The comparison directory must be new/empty. Native tools must be on `PATH`;
otherwise native integration tests are skipped. The retained local run is
`/private/tmp/mlvamaps-stage6-release`, including every worker's arguments,
FASTQs, truth, BAMs, evidence and logs. Temporary directories are not a durable
publication archive; retain that directory externally before system cleanup.
The checked-in JSON summaries preserve disagreements and hashes independently.

To profile one exact retained command, copy its `arguments.json`, change the
output argument to a fresh directory, and run:

```bash
python -m cProfile -o /tmp/mlvamaps.prof scripts/benchmark_fastq_inference.py \
  --worker /path/to/checkout /tmp/profile-arguments.json /tmp/profile-measurement.json
python -m pstats /tmp/mlvamaps.prof
```

Normal CLI examples:

```bash
mlvamaps call --database reference_build -i sr \
  --fq1 R1.fastq.gz --fq2 R2.fastq.gz --sample-mode isolate -o paired_results
mlvamaps call --database reference_build -i reads.fastq.gz \
  --read-technology ont --sample-mode isolate -o long_results
mlvamaps call --database reference_build -i reads.fastq.gz \
  --read-technology ont --fastq-max-candidates 2048 \
  --fastq-max-likelihood-cells 262144 --fastq-max-refinement-reads 100 \
  --fastq-max-alignment-bases 10000 --debug-disagreements -o audited_results
```

## Remaining publication-validation work

The blocker to validated probability, context-aware ONT correction and a
selection-bias-corrected partial-fragment model is the absence of independent
truth/calibration data in this task. Conservative fallbacks are implemented;
those statistical features are not represented as complete.

Use held-out isolates with independently established repeat and sequence alleles,
multiple ONT chemistries/basecallers, matched Illumina libraries, and measured
library geometry. Evaluate mixtures, related taxa, difficult homopolymers,
compound repeats, long expansions, reference omissions and downsampled coverage.
Keep calibration and evaluation samples separate. Report errors jointly with
call rate, exact/estimated/bounded statuses, sequence yield and assignment
ambiguity. Investigate every old/new disagreement against truth and read-level
evidence. Test whole-genome/metagenome inputs, deep target coverage and large
panels for RAM, disk and runtime scaling. Repeat the short-read performance
comparison on representative server workloads before claiming practical speed
improvements. Archive inputs, tool versions, parameters and all failures.
