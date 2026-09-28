# Short-read refactor validation

## Architecture and cleanup

Previously, the default short-read engine competitively mapped reads to candidate
alleles; an optional E/S/F/FRR engine fitted lengths before synthesizing products.
The new single pathway streams a combined context/flank index, retains mates,
classifies evidence, observes or reconstructs loci, then falls back to haploid
candidate likelihoods. Reference classification follows physical measurement.

Removed: `mlvamaps/repeat_likelihood.py`, its mixture-EM genotype fitting and
reference-filled repeat reconstruction, the old short-read candidate evidence
extractor, `--sr-engine`, unused short-read MAPQ/secondary-alignment controls,
E/S/F/FRR output columns, the per-molecule candidate-grid output, the obsolete
repeat-fitting benchmark and model-equivalence tests. Reusable native anchor,
quality and insert-statistics utilities were retained. No dependencies were
added; no dependency was exclusively required by the removed algorithm.

## Important implementation files

- `mlvamaps/short_read_recruitment.py`: combined context index, competitive
  locus assignment, mate retention, bounded process batches.
- `mlvamaps/short_read_evidence.py`: native flank/primer localization,
  motif compatibility, explicit evidence classes, robust insert estimates.
- `mlvamaps/targeted_reconstruction.py`: direct products, unique pair merging,
  bounded overlap graph, haploid candidate scoring, conservative inferred bases.
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

The final full suite passed **360 tests**, with **10 skips**: two require
minimap2 and eight require optional rapidgzip/isal backends unavailable here.
The pre-refactor baseline passed 336 tests with 11 skips. Obsolete likelihood
implementation tests were removed, useful biological cases were migrated, and
`tests/test_targeted_reconstruction.py` adds explicit acceptance checks.

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

Exact-overlap assembly is deliberately conservative with errors, ambiguous
repeat paths and large pools (256 distinct nodes/64 path expansions). Candidate
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
