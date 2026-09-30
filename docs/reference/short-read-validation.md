# Short-read refactor validation

## Architecture and cleanup

SR FASTQ commands now require a database for competitive locus recruitment.
Observed products, local stitching and fragment geometry feed the shared length
interpreter. Reference counts never fill missing observations. Classification
remains a separate database-dependent step; assembly and default long-read
calling remain independent.

The historical sample-only timing results below describe lower-level recovery
helpers. New native-mapping regressions in `test_competitive_short_reads.py`
cover database calling, novel alleles, missing primers, unmapped mates, unknown
insert sizes and cross-locus ambiguity. The reported server sample is not
available locally, so these tests do not establish recovery on that sample.

## Important implementation files

- `mlvamaps/short_read_recruitment.py`: combined context index, competitive
  locus assignment, mate retention, bounded process batches.
- `mlvamaps/short_read_evidence.py`: native flank/primer localization,
  motif compatibility, explicit evidence classes, robust insert estimates.
- `mlvamaps/targeted_reconstruction.py`: direct products, unique pair merging,
  bounded overlap graph and conservative inferred bases.
- `mlvamaps/competitive_short_reads.py`: reference calling, vectorized
  censored-boundary and phase-corrected fragment likelihoods.
- `mlvamaps/locus_reconstruction.py`: sequence-first orchestration and QC.
- `mlvamaps/locus_products.py`: existing shared assembly calibration and SNP
  markers, plus additive observed repeat-sequence/motif-edit fields.
- `mlvamaps/short_reads.py`, `short_read_mapping.py`, `unified_fastq.py`,
  `pipeline.py`, `cli.py`, `report.py`: one entry path, cleanup, output and help.
- `scripts/benchmark_cross_mode.py`: call rates, incorrect-call rate, dropout,
  per-locus metrics, SNP/repeat-sequence concordance and joined QC features.
- `scripts/benchmark_reconstruction.py`: deterministic end-to-end resource fixture.

See [method details](../concepts/locus-reconstruction.md) for recruitment,
overlap-path constraints, candidate likelihood formulas and uncertainty rules.
Every input mode uses `LocusProduct` → `genotype_product` and the existing
assembly-equivalent repeat calibration. Existing masked SNP marker semantics
are preserved; motif edits add a separate observed-repeat view.

## Tests and scientific validation

Run `python -m pytest -q` for the complete regression suite. Native Sassy and
minimap2 checks run when their executables are available; optional accelerated
FASTQ backend checks are skipped when those backends are absent.

Coverage includes short/full spans, paired/overlapping reads, targeted
multi-read reconstruction, repeats longer than reads or fragment lengths,
anchored repeat-rich mates, both one-sided boundaries and soft clips,
imperfect/cyclic motifs, repeat and flank SNPs, physical flank indels, novel
alleles, competing natural contexts, low/high depth, sequencing errors,
mixtures, absent/unidentifiable loci, adjacent-candidate uncertainty,
thread equivalence, bounded graphs and failure cleanup. Direct-evidence tests
assert that assembly and candidate scoring are bypassed.

The three-mode integration test checks identical repeat calls and canonical
markers. The repository's one-locus synthetic paired-read/assembly example
also gave the same count (4), identical masked SNP sequence (72 comparable
bases), and identical observed repeat sequence (16 comparable bases). This is
a smoke check, not a clinical or population-level accuracy estimate. No real
matched Illumina/assembly dataset was available.

## Runtime and memory

Baseline revision: `e9dd281`, using the former optional likelihood engine.
The former default competitive mapper could not run without minimap2. Both
versions used identical generated FASTQs and panels, one thread, 10,000 pairs,
six loci and 25% off-target background. Three separate processes were run per
scenario; the table reports median wall time and median process peak RSS.
Input generation/import time is excluded; output writing is included.

| Synthetic workload | Before | After | Before RSS | After RSS |
| --- | ---: | ---: | ---: | ---: |
| Direct | 1.80 s | 1.77 s | 91.5 MiB | 90.9 MiB |
| Unresolved | 2.35 s | 1.61 s | 107.7 MiB | 90.9 MiB |

[Raw measurements and environment](short-read-validation.json) are retained.
Easy-locus time is essentially unchanged; unresolved-locus time improved by
about 32%, with lower peak memory. Profiling exposed background recruitment,
unnecessary external primer searches, and large per-molecule likelihood output
as avoidable costs. These were addressed rather than accepting a several-fold
slowdown. Real-library throughput, large reference panels and classification
runtime still need independent measurement.

## Pysam regression and reconstruction-limit validation

The `ce0800d` pysam migration regressed the high-depth reconstruction fixture.
Batching independent coordinate pools into native samtools consensus, retaining
weighted molecule-aware HTSlib pileups for conflicts, and reusing unchanged
sample-learning pools restores the earlier runtime while keeping pysam.

Two sequential runs per configuration used 12,000 pairs, six loci and four
threads, without profiling or concurrent test workloads. Median call wall time
includes outputs and excludes input generation/imports:

| Synthetic workload | NumPy `6cc7f66` | Pysam `ce0800d` | Corrected pysam |
| --- | ---: | ---: | ---: |
| Direct | 2.23 s | 1.91 s | 1.88 s |
| Pileup | 11.17 s | 31.69 s | 11.06 s |

Every run retained the six expected counts (5–10); corrected HTML reports were
also checked. Individual runs, stage timings, peak memory and environment are
under `pysam_regression_validation` in
[the reconstruction measurements](reconstruction-performance.json).
The benchmark now exports counts and checks them for sufficiently covered
direct, rescue and pileup fixtures, so speed cannot hide incorrect alleles.

The full suite passed 420 tests with 10 skips. New regression coverage includes
3,500 noisy 301-base reads with variable quality trimming: the previous code
hit `assembly_node_limit`, while the corrected caller reconstructs the 432-base
product and calibrated 8-repeat allele. Other checks cover unknown/conflicting
votes, native batching, coverage above 8,000, partial primer-arm consensus,
stale measurements after template changes, and node/path limit reporting.
The reported server sample was unavailable locally; these are synthetic
regression and performance checks, not validation of its final alleles.

## Primer-linked motif selection

The DRR125654 diagnostic graph had no ready repeat graphs. Its global motif
votes frequently favored homopolymers; Bams05 and Bams15 also had competing
substitution variants. The learner now selects families supported on both
primer-linked arms, preserves longer boundary observations through the existing
pysam consensus, and recognizes cyclic substitution variants. Split repeat runs
are joined only if the observed intervening sequence remains periodic.

Eight regression cases cover unlinked homopolymers, short primer fragments,
the supplied 39-base and 9-base motif variants in separate and interspersed
blocks, unrelated boundary repeats, and repeats outside the primers. An
end-to-end FASTQ regression checks the expected 30-repeat result in the call
tables, evidence TSV and HTML report. Inferred-count fixtures supply an
independently specified insert distribution; these do not establish that depth
alone determines length. The full suite passed 428 tests with 10 optional-tool
skips. DRR125654's FASTQs were unavailable locally, so its final alleles still
require a server rerun.

Learned-boundary classification reuses the existing locus process workers,
skips reads already classified against identical arms, and caches repeated
motif alignments. Sample graph entries include
`motif_selection: primer_linked_arms` to identify the corrected learner in a
fresh run. Reinstall the updated checkout in the server's active environment
before rerunning; use `--force` for directory/manifest batches to replace
existing successful results.

Two sequential runs per configuration compared this change against `1d22002`,
using 12,000 pairs, six loci and four threads with no concurrent test workload
or profiling. Median wall times include output generation:

| Synthetic workload | `1d22002` | Primer-linked learner |
| --- | ---: | ---: |
| Direct | 1.98 s | 1.77 s |
| Pileup | 10.98 s | 11.09 s |
| Rescue | 2.46 s | 2.34 s |

All runs retain the expected counts 5–10. The pileup workload now learns six
usable repeat graphs instead of two; its median runtime differs by about 1%.
These small synthetic comparisons do not establish server throughput. Raw
timings, stage times, memory, graph readiness and call checks are retained under
`primer_linked_motif_validation` in
[the reconstruction measurements](reconstruction-performance.json).

## Reproduction commands

```bash
python -m pytest -q -rs
python scripts/benchmark_reconstruction.py --output bench/direct --molecules 10000 --threads 1
python scripts/benchmark_reconstruction.py --output bench/unresolved --molecules 10000 --threads 1 --scenario unresolved
python scripts/benchmark_cross_mode.py manifest.tsv --output comparison.json
python scripts/benchmark_short_read_recruitment.py --reads1 R1.fq.gz --reads2 R2.fq.gz --primers panel.tsv --pairs 10000 --threads 1 4 --output recruitment.json
```

The cross-mode manifest has `sample_id`, `mode`, `outdir`; use the same panel
for independent assembly and FASTQ runs. No assembly-derived answers are used
by recruitment, reconstruction or likelihood scoring.

## Remaining limits

Overlap assembly tolerates 2% substitutions and remains conservative with
ambiguous repeat paths and large pools (2,048 compacted nodes/64 path expansions). Candidate
backgrounds are bounded at 64. Long pure repeats may remain ambiguous without
reliable empirical insert statistics. Short context flanks often provide too
few ordinary pairs for estimating those statistics; an independently measured
library override is supported.

Confidence scores and mixture thresholds are not externally calibrated.
Recruitment still uses short seeds and may become less selective on large
panels. Retained evidence memory grows with on-target depth. Inferred repeat
bases remain unknown, and partial flank SNPs do not provide fully phased
haplotypes. Motif-relative indels may have multiple equivalent placements;
nonidentical repeats above 10 kb retain sequence without global motif alignment.

Sub-1× controls are reproducible with `scripts/benchmark_low_coverage.py`.
See [the benchmark instructions](../../scripts/README.md) for paired 150-bp
sampling, truth comparisons and a server-side Shovill comparison protocol.
No real-data superiority claim follows from these synthetic controls.

## Competitive sample model validation (September 2026)

The full regression suite passes 466 tests, with 8 optional-backend skips.
Packaging and the single-path CLI were checked. Cross-mode regressions compare repeat counts,
product lengths, SNP sequences and combined markers for singleton reads,
overlapping pairs, surviving orphan mates and primer indels. Primer-only panels
and an unused, incompatible database candidate index are included: database
contents must not change allele measurement.

On the same six-locus, 12,000-pair direct-evidence fixture with two threads, the
previous reference-candidate mapper took 9.45 seconds; the competitive sample
model took 2.63 seconds, with the same 5–10 repeat calls. The harder primer-only
depth fixture dropped from 27.87 to 19.74 seconds with unchanged provisional
estimates. These are local synthetic timings, not whole-genome or server claims.

The sub-1× sweep used paired 150-bp reads, 100 seeds per coverage and identical
reads for both panels. Each cell below is exact definitive calls / all 600
locus evaluations, including missing and uncertain results:

| Coverage | Rich panel | Calibrated primer-only panel |
| --- | --- | --- |
| 0.125× | 34 / 600 | 13 / 600 |
| 0.25× | 57 / 600 | 24 / 600 |
| 0.5× | 95 / 600 | 44 / 600 |
| 0.75× | 123 / 600 | 72 / 600 |

All 462 definitive repeat calls matched assembly mode; 4,338 evaluations were
missing or uncertain. Provisional numeric estimates are reported separately
and are not counted as definitive matches. The two panels share reads, so
these evaluations are not independent biological samples. This supports sparse
observation handling, not high recall below 1× or superiority over Shovill.
Run the documented server comparison on real matched inputs before making
that claim.

SNP validation uses a neutral alignment score for `N` in both canonical output
and the comparator, so missing bases cannot become deletion calls. Partial
spans retain their measured repeat boundaries, including half units, instead
of allowing a motif-only mask to leak repeat bases into the SNP marker. The
benchmark compares each FASTQ run to assembly using the same panel, preserving
its SNP masking convention.
