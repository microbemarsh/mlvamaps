# Mandatory FASTQ repeat counts: implementation and validation

This is a historical validation record. Current FASTQ runs retain these
predictions diagnostically and leave rejected primary counts blank; see the
[current FASTQ acceptance rules](fastq.md). The results below describe the
earlier mandatory-count behavior and have not been reinterpreted as accepted calls.

Every successfully completed FASTQ run now contains one finite, nonnegative
repeat-count prediction per configured locus when an eligible locus prior is
available. Confidence thresholds change the evidence classification rather than
erase the number. No-data predictions remain explicitly prior-only, with zero
sample-evidence confidence. Assembly inference and assembly output schemas are
unchanged.

Numerical completion is not genotyping accuracy. In particular, short reads
that do not span a large repeat cannot determine its full length merely because
the software is required to report a number. The difficult synthetic results
below retain these wrong predictions and include all truth-known loci.

## Baseline and validation method

The baseline is commit `6851a68` (`Beef up fastq inference`). Before editing,
the native-tool suite passed **422 tests** and both source copies were run on
the same five synthetic layouts. The baseline source archive and initial
comparison remain at `/private/tmp/mlvamaps-complete-count-baseline` and
`/private/tmp/mlvamaps-complete-count-initial` in the development environment.

Final validation uses the existing organism-independent concordance and challenge
fixtures. Clean layouts contain 11 known alleles. The challenge has 15 configured
loci: 14 with known repeat truth, including a zero-coverage locus, and one absent
locus with no truth allele. It includes a 400-unit repeat, partial and interrupted
units, sequencing errors, low coverage, indistinguishable loci and a 5/9 mixture.
Either true mixture allele counts as correct. No difficult loci are excluded.

Each comparison runs the baseline and updated source in alternating serial fresh
processes on identical input bytes, with two native threads. Default timings are
medians of three trials at two fixture copies; additional coverage, metagenome
and forced-limit checks use one trial each. No tests ran concurrently with the
reported timings. Wall and CPU exclude imports and fixture preparation. Peak
RSS is `max(parent, largest child)`, **not the sum of concurrent processes**.
The [machine-readable record](../reference/fastq-complete-count-benchmark.json)
retains source/input hashes, all trials, distributions of errors, calibration
bins, coverage strata and baseline/new disagreements.

## Changed functions and preserved behavior

| Area | Functions and behavior |
| --- | --- |
| Shared inference | `infer_alleles`, `_infer_measurements`: retain the selected count for ambiguous, mixed and low-support calls; enforce bounded candidate work. Limits use retained observed modes and a precision-weighted candidate, scored without allocating the rejected dense matrix. |
| Incomplete evidence | `estimate_unresolved_calls`: joint likelihood across independent molecules, marginalizing distinct conflicting measurements or tied calibrated pair lengths within a molecule; maximum alignment score per candidate avoids duplicate-reference votes. Paired lengths precede spanning alignments, partial bounds and reference-conditioned evidence. |
| Completion | `locus_repeat_prior`, `finalize_repeat_count_calls`: deterministic eligible priors, partial-bound conditioning, one finite count per unique locus, explicit configuration errors, confidence and provenance fields. |
| Partial calibration | `partial_repeat_lower_bound`: uses the complete configured flank length for product calibration, even when only its terminal 21 bases anchor the search. |
| Orchestration | `run_unified_fastq_inference`: finalizes after fragment and consensus reconstruction. `run_call`: also finalizes when no reads survive QC. Mapping and extraction systems are reused. |
| Reconstruction | `apply_reconstructed_calls`: consensus disagreement lowers confidence while retaining the selected molecule estimate. |
| Reporting | `write_evidence_reports`, `common_calls_to_compatibility`, `simple_call_rows_from_alleles`, `_allele_rows`, `build_fingerprint`: propagate counts and uncertainty; depth zero cannot erase a prediction; presence remains independent. HTML labels prior-only predictions. |
| Biological comparisons | `match_profiles`, `profile_match_locus_rows`, `build_call_matrix`: exclude prior-only values. The filtered matrix also supplies SNP-plus-repeat distances. Classification and per-sample profile trees never receive prior-only predictions as measurements. Exported numerical profiles retain primary counts. |
| Benchmark | `challenge_metrics`, `compare_versions`: all-locus exact accuracy, completion across configured loci, evidence strata, error and support diagnostics, explicit missingness, metagenome and resource-limit arguments. |

No dependency, mapping pass, new external service or organism-specific model was
added. Existing recovery, assembly PCR calibration, molecule identities and
reference classification are retained.

## Prior and likelihood assumptions

The prior uses equal mass on distinct observed reference alleles. Reference
copies and multiple backgrounds with the same count do not increase that mass.
Without observed alleles, a positive calibrated panel nominal count is eligible;
otherwise an explicit non-default repeat range defines a discrete uniform prior.
The default `0..100` search range and default nominal zero are not informative.
Observed zero alleles and an explicit zero range remain valid. Ties select the
smaller count. Very large priors use a bounded representative grid and preserve
the full prior envelope in the reported interval.

Exact anchored terminal repeat runs impose one lower-bound constraint per
molecule. The greatest bound conditions the prior. An exponential continuation
with scale `max(1, prior maximum - prior minimum)` permits expansions outside the
catalog; at most 32 tail states are retained, within the candidate budget. The
mode can be a fractional bound in calibrated repeat units. Without a prior, an
informative bound can supply the minimum compatible estimate. This continuation
is an explicit heuristic; it does not identify an unobserved upper boundary or
create a validated biological maximum. More depth of the same partial sequence
does not make its unsequenced repeat length measurable.

Paired-fragment inference retains the existing calibrated Gaussian geometry
model and its sufficient-statistic reduction. Candidate or cell limits use a
bounded likelihood over retained measurements rather than filling from a prior.
Long-read error sigma, boundary flags and refinement limits lower evidence
confidence and preserve/widen measurement envelopes; the maximum-likelihood
count remains visible. Isolate and metagenome inference retain their existing
primary models. Incomplete-evidence fallbacks do not establish strain fractions.

`model_support` is relative to represented candidates. `data_evidence_confidence`
is an uncalibrated evidence-strength score, discounted by independent depth,
source reliability, `1/(1 + sigma)`, quality flags, ambiguity and resource limits.
It is exactly zero for prior-only predictions, even with model support one.
Intervals are observation or bounded candidate envelopes, not calibrated
confidence intervals. Source reliability factors are explicit heuristics:
measurement 1, reconstruction .5, spanning alignment .4, paired length .35,
conflicting measurement .25, and partial/reference-conditioned evidence .1.

## Regression coverage

The full suite passes **439 tests**, including existing assembly tests. The new
checks cover zero-depth export propagation, default-prior rejection, genuine
zero priors, duplicate-reference invariance, ambiguous-locus exclusion, joint
conflicting measurement and paired-length likelihoods, partial-bound calibration and conditioning, large noisy
measurements under candidate/cell limits, measurement uncertainty, and prior
exclusion from profile/distance comparisons.

Native challenge tests exercise both read workflows in both sample modes. They
verify every configured locus, agreement across all per-locus and fingerprint
exports (including MYOGA and long-read allele calls), zero evidence confidence
for absent/zero-coverage loci, retention of the 400-repeat long-read measurement,
and byte-identical core tables across repeated runs. The empty-FASTQ regression
also checks complete prior-only long-read exports.

## Accuracy, completion and remaining errors

All three clean layouts retain **11/11 correct numeric calls** in both versions.
For the default challenge, completion includes all 15 configured loci, while
exact accuracy uses all 14 truth-known loci, counting a missing baseline call as
not correct. The absent locus has no allele truth and remains `present=no` in
both versions despite its new prior prediction.

| Challenge metric | Short baseline | Short updated | Long baseline | Long updated |
| --- | ---: | ---: | ---: | ---: |
| Numeric completion / configured loci | 6/15 (40%) | 15/15 (100%) | 9/15 (60%) | 15/15 (100%) |
| Exact counts / all truth-known loci | 5/14 (35.7%) | 6/14 (42.9%) | 9/14 (64.3%) | 11/14 (78.6%) |
| Missing truth-known predictions | 8 | 0 | 5 | 0 |
| Accuracy among numeric predictions | 83.3% | 42.9% | 100% | 78.6% |
| MAE across all truth-known loci | Undefined: 8 missing | 29.544 | Undefined: 5 missing | 0.286 |
| MAE among numeric predictions | 0.167 | 29.544 | 0 | 0.286 |

The large increase in short-read numeric MAE is retained, not hidden by a
confidence filter. Baseline MAE excludes eight difficult loci; it is not an
all-locus accuracy benchmark. None of the previously emitted counts changed on
the five default layouts. The short-read 400-unit array is predicted as **17.5**
from an anchored partial run, contributing an error of 382.5 units; the 40-unit
expansion is also predicted as 17.5. Those fragments do not span the repeat, and
the finite reference prior cannot supply its unknown length. More completion
does not repair that lack of information. Long reads retain the measured 400-unit
count, including under both forced computational limits.

Other wrong short-read predictions are the existing homopolymer call 15 instead
of 16, a partial-bound 7.8889 instead of 8, a noisy partial-bound 12.5 instead of
16, and the three truth-known prior-only predictions. Both read types predict 4
instead of 6 at zero coverage and 4 instead of 5 at each indistinguishable locus.
These three prior-only errors have **zero sample-evidence confidence**. The
absent locus also receives prior 4, with zero evidence and no presence claim.

Updated accuracy by evidence class, retaining every truth-known locus:

| Evidence class | Short correct / total | Long correct / total |
| --- | ---: | ---: |
| EXACT | 4/4 | 4/4 |
| ESTIMATED | 1/2 | 6/6 |
| LOWER_BOUND | 0/4 | — |
| POSSIBLE_MIXTURE | 1/1 | 1/1 |
| PRIOR_ONLY | 0/1 | 0/1 |
| AMBIGUOUS_LOCUS (also prior-only) | 0/2 | 0/2 |

## Runtime and peak memory

Environment: macOS arm64, Python 3.13.15, minimap2 2.31-r1302, Sassy 0.2.1;
dependency versions are recorded in the JSON artifact. Each default value is a
median of three serial fresh-process trials.

| Workflow | Wall seconds baseline → updated | CPU seconds baseline → updated | Peak MiB baseline → updated |
| --- | ---: | ---: | ---: |
| Clean paired | 0.142 → 0.141 | 0.094 → 0.094 | 77.34 → 77.09 |
| Clean shotgun | 0.995 → 0.968 | 0.645 → 0.658 | 77.34 → 76.98 |
| Clean long | 0.679 → 0.680 | 0.565 → 0.567 | 80.08 → 80.42 |
| Challenge paired | 0.635 → 0.608 | 0.486 → 0.486 | 82.94 → 83.41 |
| Challenge long | 4.973 → 5.043 | 4.600 → 4.542 | 183.98 → 184.11 |

No major runtime or memory regression appeared on these workloads. Small wall-time
differences on these fixtures do not establish a production speedup. Core artifact
hashes are identical across all three default trials for
each version. All **80 benchmark workers succeeded**, and each updated source
hash matches the final implementation.

## Coverage, mixture mode and forced limits

The coverage sweep changes fixture multiplicity and retains independent seeded
error realizations. It is a sensitivity check, not a controlled real-sample
coverage calibration. Zero-coverage and absent loci remain at zero coverage.

| Copies (typical clean-locus molecules) | Short exact / 14 baseline → updated | Long exact / 14 baseline → updated | Updated short / long MAE |
| --- | ---: | ---: | ---: |
| 1 (12) | 6 → 7 | 7 → 11 | 30.339 / 0.286 |
| 2 (24) | 5 → 6 | 9 → 11 | 29.544 / 0.286 |
| 4 (48) | 6 → 7 | 10 → 11 | 29.944 / 0.286 |

All updated configurations complete 15/15 loci. The default JSON additionally
stratifies accuracy by coverage 0, 2, 20, 24 and 40, motif length, repeat length
and error rate. Higher coverage does not resolve indistinguishable loci or
unspanned arrays. At four copies, challenge long wall time is 10.849 → 9.678 s
and peak RSS is 299.14 → 300.34 MiB, illustrating the remaining depth-dependent
memory use in retained mapping/extraction evidence.

| Additional configuration | Short exact / 14 baseline → updated | Long exact / 14 baseline → updated | Updated short / long MAE |
| --- | ---: | ---: | ---: |
| Metagenome | 6 → 6 | 9 → 11 | 29.544 / 0.286 |
| Maximum candidates = 2 | 1 → 5 | 0 → 10 | 29.579 / 0.321 |
| Maximum likelihood cells = 1 | 1 → 7 | 0 → 11 | 29.472 / 0.286 |

Forced-limit updated completion remains 100%; the baseline supplies one short
count and no long counts. The 400-repeat long measurement remains 400 and carries
`candidate_limit` or `likelihood_cell_limit`. Candidate-limited long wall time
is 4.968 → 4.936 s and peak RSS 182.59 → 183.00 MiB; cell-limited long time is
4.963 → 5.037 s and RSS 181.06 → 182.55 MiB. These intentionally tiny budgets
alter candidate choice and can occasionally improve a synthetic answer by
chance; they are not recommended accuracy settings.

## Confidence diagnostics and limitations

This small synthetic set does **not** validate calibrated confidence. Updated
model-support Brier scores are 0.383 (short) and 0.242 (long), while evidence-score
Brier diagnostics are 0.128 and 0.239. Baseline support diagnostics cover only
its nonmissing predictions and therefore cannot be compared as if the evaluated
populations were identical. The baseline has no data-evidence confidence field;
its corresponding diagnostic is unavailable, not fabricated as zero.

The short-read evidence-score bins have accuracy 2/9 below .25, 0/1 from .25 to
.5, 2/2 from .5 to .75, and 2/2 above .75. Long-read bins have 3/6, 2/2, 2/2 and
4/4 respectively. This exposes a wrong homopolymer estimate with moderate
evidence support and deliberately conservative scores for some correct noisy
long-read measurements. Relative model support can be one for an incorrect
single-reference prior; it must never be read as validated accuracy.

Remaining limitations are the explicit partial-tail prior, catalog ascertainment,
uncalibrated pair-selection bias and ONT errors, incomplete-evidence mixture
fractions, and depth-dependent extraction memory. A bounded partial interval
can exclude the true 400-unit expansion because its upper endpoint is only the
end of the represented candidate grid. External consumers of wide fingerprints
must join evidence fields themselves; a wide numeric fingerprint alone cannot
distinguish observed from prior-only genotypes. Matched real-sample truth and
external validation are still required before any accuracy or calibration claim.

## Reproduction

With the existing development environment and minimap2/Sassy on `PATH`:

```bash
git archive 6851a68 | tar -x -C /path/to/empty-baseline-directory
python -m pytest -q
python scripts/benchmark_fastq_inference.py \
  --baseline-source /path/to/empty-baseline-directory \
  --workdir /path/to/new-empty-comparison-directory \
  --challenge --copies 2 --repetitions 3
```

Use new comparison directories for each additional run: `--sample-mode metagenome`,
`--max-candidates 2`, `--max-likelihood-cells 1`, or `--copies 1` / `--copies 4`.
Retained final benchmark directories are named `/private/tmp/mlvamaps-count-final-*`.
They contain identical FASTQs for each version, complete output tables, argument
files, logs, timings and truth. Configuration failures remain failures rather
than being counted as completed analyses.
