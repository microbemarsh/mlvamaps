# Illumina short-read workflow

Illumina mode uses the shared FASTQ architecture. Paired reads are recruited by
unique flank seeds during QC, retaining both mates independently of competitive
minimap2 mapping. Direct sequence measurements enter shared allele inference;
uncertain loci can use a bounded local graph. Single-end reads and QC-surviving
orphan mates also participate in recruitment, direct measurement and local
assembly. Use `--short-no-assembly` to disable local graph and consensus
reconstruction while retaining read mapping and molecule likelihoods.

Recruitment uses both observed reference flanks and configured panel anchors.
This includes primer-panel loci without a database product; complete reads or
overlapping pairs can still supply a measured count. Every panel locus remains
in the results, with a blank accepted count and an explicit reason when the
sample evidence cannot resolve it.

Current FASTQ outputs prioritize supported calls over completeness. Accepted
`repeat_count`/`repeat_count_raw` values are otherwise blank;
`predicted_repeat_count`/`predicted_repeat_count_raw` retain diagnostic candidates,
with `call_accepted` and `call_rejection_reason` explaining the decision. Low-depth,
prior-only, partial/unbounded, conflicting and indirect estimates do not enter
fingerprints or repeat distances. Direct measurements or supported reconstruction
must meet the existing depth, model-support and margin criteria; those scores are
not calibrated accuracy guarantees. Presence support and assembly behavior are
unchanged. See the [output reference](../reference/outputs.md#complete-fastq-predictions-and-confidence).

## Mapping-only paired-read inference (2026-11-03)

The dated validation results below predate the current acceptance policy and
retain their original mandatory-count accuracy denominators. Unbounded and
indirect rescue predictions now remain diagnostic estimates with blank accepted
counts; the historical results do not establish current cohort concordance.

For an isolate, run:

```bash
mlvamaps call --database reference_build -i sr \
  --fq1 sample_R1.fastq.gz --fq2 sample_R2.fastq.gz \
  --sample-id sample --sample-mode isolate --short-no-assembly \
  -o results/sample -t 8
```

This path uses competitive mapping, direct read measurements, and a conditional
likelihood over paired-read placements. It builds no local graph or consensus
sequence. Direct measurements keep precedence; the new pair rescue applies only
to unresolved isolate loci without measured or mixed-allele support.

Both mates contribute to the fit, including weak opposite junctions excluded by
the earlier 21/12-base anchor rule. An exact primary flank anchor of at least
12 bases supplies locus identity, subject to competition across configured
loci. The model considers repeat-mate and boundary placements together, includes
one-substitution alternatives in the opposite mate, and conditions the fragment
likelihood on accepted placements. At least two informative junction molecules
on each side are required. Mates, duplicate reference contexts, and repeated
observations of one original molecule cannot multiply support.

External or valid unique-flank fragment calibration takes priority. When that
is unavailable, a bounded sample of 5,000 QC-passing molecules can estimate a
narrow overlapping fragment component from unique, exact, nonperiodic overlaps.
It requires more than 100 retained molecules and rejects broad, heterogeneous,
multimodal, or detection-edge-censored overlap distributions. The SD has a 1 bp
floor. This component is recorded separately in
`overlap_fragment_calibration.tsv`; it never fills the whole-library calibration
in `library_insert_distribution.tsv`. Fragments too long to overlap remain
unobserved, so even an accepted component cannot exclude a longer library tail.

The unchanged difficult fixture's expansion improves from **17.5 to 40 repeats**,
bringing individual-locus count agreement from **9/14 to 10/14**, without assembly.
The estimate assumes a pure repeat array and the calibrated fragment component.
Longer interrupted arrays can still explain these reads. The result therefore
remains low-confidence `ESTIMATED`, with `upper_unbounded=True`,
`finite_length_support=False`, and a blank `repeat_interval_high`. It generates
no assembled query product and no measured classification votes. The 400-repeat
array, zero-coverage locus, and individually indistinguishable loci remain
unresolved; real-sample agreement still needs matched server data.

`joint_pair_diagnostics.tsv` records the rescue's support, calibration source,
competition, fit outcome, and reasons for abstaining. The `pair_rescue_*` call
fields preserve its provenance, with method `joint_pair_read_likelihood`
(`conditional_joint_pair_likelihood` in the common inference table). Search and
likelihood-cell limits bound the additional work; these are computational limits,
not biological upper bounds.

Three serial fresh-process trials per workflow used identical FASTQ bytes and
alternated baseline/updated order. The baseline includes the preceding
shared-group and conditional-length changes; mapping-only is the updated
workflow with local assembly and consensus disabled.

| Fixture | Correct counts before → mapping-only | Before median seconds | Updated default | Mapping-only |
| --- | ---: | ---: | ---: | ---: |
| Clean paired reads | 11/11 → 11/11 | 0.215 | 0.199 | 0.209 |
| Clean shotgun reads | 11/11 → 11/11 | 1.596 | 1.746 | 1.754 |
| Original difficult paired reads | 9/14 → 10/14 | 1.037 | 1.087 | 1.007 |

All audited output tables were deterministic across the three trials. The full
native suite passed **647 tests**, with **13 pair-model tests** rerun after the
final resource-accounting correction. The native integration check raises if
either graph assembly or consensus is invoked in mapping-only mode. It also
checks that a stronger existing lower bound cannot be overwritten. These small
synthetic timings do not establish whole-genome throughput or server accuracy.

Appending 50,000 independent random 180 bp background fragments (100 bp mates)
to the unchanged challenge retained the **9/14 → 10/14** improvement. Three
alternating trials measured **2.660 → 3.877 seconds** median, an additional
1.217 seconds (46%). Maximum recorded RSS was **84.91 → 93.77 MiB**; this is the
maximum of the parent and largest child, not concurrent total memory. Output
tables were deterministic. This workload checks background-read overhead and
the bounded calibration pool; it is not a whole-genome accuracy benchmark.

The [validation record](../reference/short-read-pair-rescue-validation.json)
retains source/input hashes, individual trials, and remaining mismatches.

## Earlier shared-group and conditional-length validation (2026-11-03)

The preceding implementation added a shared **5-repeat** allele for
`ambiguous_a;ambiguous_b`, supported by 40 original molecules. Both member loci
remain ambiguous and retain their existing per-locus predictions. The original
fixture's individual-locus accuracy remained **9/14** at that stage: the new
conditional-length model had no usable independent library calibration for its
fixed 180 bp fragments. The mapping-only pair rescue above subsequently improves
that result to **10/14**. The fixture and truth were unchanged.

A separate native test generates 100 bp mates from a known Gaussian fragment
library and retains the two observation classes modeled by the caller. It uses
a true 40-repeat expansion with a 4-repeat reference. Supplying the same
260 bp mean and 30 bp SD to both versions changes the prediction from **19.5 to
40 repeats**, with likelihood support from 39 to 40.75. The result remains
`ESTIMATED` and supplies no measured classification votes. Without calibration,
the updated workflow retains the lower-bound prediction and unknown upper length.
This controlled test validates the length calculation; it is not a whole-genome
or real-sample accuracy benchmark.

Three alternating fresh-process trials on identical input files measured:

| Fixture | Correct counts before → after | Before median seconds | After median seconds |
| --- | ---: | ---: | ---: |
| Clean paired reads | 11/11 → 11/11 | 0.181 | 0.155 |
| Clean shotgun reads | 11/11 → 11/11 | 1.713 | 1.747 |
| Original difficult paired reads | 9/14 → 9/14 | 0.976 | 1.004 |
| Calibrated expansion, 1,007 fragments | 0/1 → 1/1 | 0.536 | 0.638 |

Outputs were deterministic across all three trials. The calibrated test added
about 0.10 seconds and used a maximum recorded RSS of 78.25 MiB versus 76.88 MiB
before; RSS is the maximum of the parent and largest child, not their sum. The
changes reuse the existing mapping and add no dependencies. The full native
suite passed **622 tests**, including duplicate-molecule, mixed-group, resource
limit, discrete-fragment, and unbounded-tail regressions. The
[validation record](../reference/short-read-joint-length-validation.json) retains
source/input hashes, timing trials, shared group calls, and remaining mismatches.
Matched server data are still needed to assess real-sample agreement.

## Earlier local recovery validation (2026-11-03)

The earlier bounded-graph change recovered the known product from a redundant
pool of 4,100 fragments that previously stopped at `fragment_limit`. It uses 82 graph
representatives and retains all 4,100 original support IDs. A separate regression
retains a competing branch supported by only two copies, preventing the reduction
from manufacturing a unique product. New partial-read checks cover right-flank
anchors, motif phase, calibrated fractional bounds, and reverse-oriented orphan
reads. Insert-calibration diagnostics distinguish unavailable and supplied
calibration without changing the acceptance rules.

Three alternating fresh-process trials on identical synthetic inputs, using two
copies and isolate mode, produced:

| Fixture | Before exact | After exact | Before median seconds | After median seconds |
| --- | ---: | ---: | ---: | ---: |
| Clean paired reads | 11/11 | 11/11 | 0.158 | 0.149 |
| Clean shotgun reads | 11/11 | 11/11 | 1.582 | 1.582 |
| Difficult paired reads | 9/14 | 9/14 | 1.002 | 0.936 |

Core output tables were deterministic across repeated trials in both versions.
The full native-enabled test suite passed 583 tests.
These small fixtures establish regression behavior, not server-cohort accuracy
or a general speedup. At that stage, the difficult fixture's unspanned expansions
and uninformative loci remained unresolved by sample length evidence. The
[validation record](../reference/short-read-recovery-validation.json) retains
source/input hashes, every timing trial, agreement metrics, and the deep-pool
recovery check.

## Command line

All FASTQ modes require a completed `--database` for the target taxon and panel.
The saved panel is used automatically when `-p` is omitted.

```bash
mlvamaps call -p panel.tsv -i sr \
  --fq1 SRR000001_1.fastq.gz \
  --fq2 SRR000001_2.fastq.gz \
  --profiles profiles.tsv \
  --database reference_build \
  --sample-metadata metadata.tsv \
  --sample-id SRR000001 \
  -o results/SRR000001 -t 8
```

For single-end data, omit `--fq2`. Mates are not inferred from filenames.
Interleaved FASTQ is not supported. Compressed files are read without whole-file
decompression. Pair files must have equal record counts and matching normalized
IDs at every record.

For multiple paired samples in one directory, filenames can supply the pairing:

```bash
mlvamaps call --database reference_build -i short_read_directory/ --short-reads \
  --sample-metadata metadata.tsv -o results -t 32
```

This recognizes exact `PREFIX_1.fastq.gz` and `PREFIX_2.fastq.gz` suffixes,
uses `PREFIX` as the sample ID, and fails before analysis if either mate is
missing. Discovery is non-recursive and ignores unrelated filenames.

## QC

The defaults require 40 post-trim bases and mean Q15. Three-prime trimming is
disabled unless `--short-trim-quality` is set. With the default
`--short-min-pair-retention 0.5`, a good mate remains as an orphan if the other
mate fails. IDs and mate association are not rewritten.

`short_read_qc_summary.tsv` reports input, rejected, retained, and orphan
counts. Insert-size fields remain empty unless external calibration is supplied:
candidate-relative template lengths do not provide an independent estimate of
a sample's library geometry.

## Flank recruitment and mate rescue

The existing `pysam.FastxFile` parser streams paired files once through QC.
`normalize_read_id` removes mate suffixes and validates synchronized names.
During this same pass, a compiled literal-prefix trie searches both mates for
21-base seeds from observed reference sequences outside their repeat intervals. Seeds
shared between loci, non-concrete/low-complexity seeds, and seeds found inside
configured repeat sequences are excluded. Both orientations are searched.

Either mate can identify a locus. Both QC-passing reads are retained, including
repeat-only, soft-clipped, or unmapped mates. A dictionary keyed by normalized
fragment name prevents duplicate seed hits, mates, or repeated records from
multiplying support. Conflicting flank assignments are not recruited. Existing
QC rejection thresholds still apply; mate rescue does not bypass quality QC.
Seed uniqueness is checked against configured locus contexts, not every genome
in a reference cohort. Background homology can therefore still recruit reads;
sequence measurements, a supported pair-likelihood estimate, or a covered
two-flank graph must support sizing.

The same observed-reference mapping is reused for reference classification. In short-read
allele inference, an exact flank assignment excludes competitors at other loci.
Mapper-only recovery requires at least 21 aligned bases outside the repeat and
90% alignment identity after excluding repeat-region length indels. Repeat-only mappings cannot independently identify a
locus. Fragment IDs also remain the unit of support in shared inference.

## Direct reference-locus mapping

Targets are actual observed sequences from the required database's per-locus
FASTAs. Duplicate sequences within a locus retain all source reference IDs.
Filtered mates map competitively to these targets with minimap2, without an early
taxon restriction. Both equivalent reference mappings and sequence mismatches
remain available to reference classification.

The sample's `candidate_mapping/` directory records the observed FASTA, metadata,
and provenance. Existing database candidate indexes remain compatibility assets;
they do not determine the sample's repeat count. Reference classification reuses
these alignments and the measured fragment counts without another FASTQ mapping.

## Direct VNTR inference

Mapped mates now use the same locus-window refinement as long reads. Mapping
supplies orientation and local search coordinates. Complete primer products take
the fast path; unresolved near-complete products receive bounded flank realignment.
Partial mates retain local measurement windows and, when enabled, the fragment
graph, avoiding unnecessary full-product realignment.

With consensus reconstruction enabled, supported complete segments containing
primer sequencing errors can use shared SPOARS correction before assembly PCR.
Different measured repeat alleles are clustered separately. A consensus that
conflicts with a confident molecule call remains unresolved; it does not add
independent molecule votes. `query_locus_reconstruction.tsv` records the support
and whether the final product was exported.

Direct measurements use both original mates, regardless of which candidates
retained their alignments. Overlapping pairs use the existing merger, with a
native prefix precheck and a requirement for a unique compatible overlap.
Repeat-only overlaps remain excluded. Native fuzzy-anchor prechecks avoid
launching Sassy for reads that cannot contain a complete primer pair or both
repeat flanks; the existing measurement function still makes the final call.

When at least the configured minimum number of directly measured fragments
agree on one count, weaker candidate-relative mappings remain presence evidence
and cannot move that count. The default minimum is three fragments.

Raw fragment lengths and measurement uncertainty enter a single-allele model
in isolate mode, or an allele-mixture EM model in metagenome mode. The selected
raw allele is then rounded using the assembly convention. Mates and duplicate
reference alignments never add independent molecule votes.

If library fragment lengths have been calibrated independently, supply
`--short-insert-mean BP --short-insert-sd BP` together. Mean and standard deviation
describe the full library fragment length, including both reads. In isolate mode,
unresolved loci can use a bounded joint length model combining opposite-flank
pairs, partial boundary pairs, anchored lower bounds, and repeat-only mates of
uniquely anchored reads. These remain diagnostic predictions, with
method `joint_short_read_length` (`short_read_joint_length_likelihood` in the
common inference table); they cannot replace direct
measurements or create measured support for classification.

Partial pairs require one exact 21-base anchor. The opposite boundary can use an
exact 12-base anchor solely for sizing. The model conditions the calibrated
Gaussian fragment distribution on the accepted observation classes, correcting
for fragment-length selection at the target. Boundary pairs and repeat-only
mates share one selection normalizer. The model compares zero and 1% background
components to allow anomalous pairs without imposing that fraction on clean
data. This choice is uncalibrated and exported as `length_null_fraction`.
Exact repeat tails are required, so interrupted or noisy repeats
may supply no such observation. Unassigned repeat-only reads do not supply locus
coverage or abundance. Targeted libraries and incorrect insert calibration can
still violate the model's assumptions. Metagenome mode retains its existing
paired-length fit; the new single-allele fit does not apply to mixtures.

Candidate and likelihood-cell limits bound this calculation, with repeated
observations scored together and every original molecule counted once. A finite
search limit is reported as `candidate_grid_high`, separately from length
uncertainty. Bounds alone, repeat-only mates alone, or a likelihood still supported
at the upper search edge leave `repeat_interval_high` blank and
`upper_unbounded=True`. A diagnostic point prediction can still depend on
the locus prior; it does not establish an upper length. Resource-limited fallback
retains a bound-constrained estimate with zero model support and unknown upper
length.

Without explicit settings, automatic
calibration requires at least 100 independent inward primary pairs wholly
within unique nonrepeat flanks, MAPQ at least 20, near-complete high-identity
alignments, nonzero trimmed dispersion, and sufficient distance from target
edges. Short amplicon catalogs commonly provide no eligible calibration.
Legacy paired evidence without selection geometry retains the original Gaussian
fallback and its target-selection limitation.

`library_insert_distribution.tsv` records the effective calibration and its
source. It also reports `calibration_status`, eligible and retained fragment
counts, and the required minimum. Failed automatic calibration distinguishes no
alignments, insufficient unique-flank pairs, insufficient pairs after trimming,
zero dispersion, and target-edge censoring. User-supplied calibration is labelled
separately. These diagnostics do not relax calibration requirements or create
an insert estimate when the targets lack suitable nonrepeat sequence.
The legacy QC and run-metadata insert fields still describe only
user-supplied settings. See [evidence validation](fastq-evidence-validation.md)
for current status definitions, limits, retained errors, and measured tradeoffs.

Reads containing a uniquely matched left or right repeat flank and an exact
terminal repeat run can supply a lower bound when the other boundary is absent.
Right-anchored reads may start or end within a motif unit; the visible repeat
bases are converted to the same calibrated allele units as left-anchored reads.
Both read orientations are checked, and mates still count as one molecule.
This evidence bounds an estimate; it does not establish the full repeat length.

## Shared alleles at indistinguishable loci

When loci have identical observed sequence catalogs, boundaries, primers, and
allele calibration, already ambiguous molecules can be measured as one group.
`locus_group_calls.tsv` and the HTML report show the shared allele and independent
molecule support. Duplicate references and mates do not multiply that support.
Different measured group alleles remain separate, since the member loci can
carry different repeat lengths even in an isolate.

Group calls leave individual loci ambiguous. They do not assign presence or an
allele to every member, enter the per-locus fingerprint, or supply classification
votes. Group fractions combine loci and do not establish strain fractions.
Molecules with competitors outside the group are excluded. Missing group evidence
produces an unresolved group row without a prior-derived allele.

## Conditional micro-assembly

This section describes the optional reconstruction path. `--short-no-assembly`
skips it and disables consensus reconstruction; mapping and the paired-read
likelihoods described above remain available.

The graph is attempted only when a locus has flank-recruited fragments and
lacks an agreeing direct count supported by the configured minimum fragments.
This includes absent direct spans, insufficient direct support, and conflicting
direct measurements. Strong direct calls skip it, including calls with some
soft-clipped or unmapped mates.

The in-process de Bruijn graph uses both mates, native overlapping k-mer
extraction, and a small Python adjacency dictionary. Each fragment contributes
at most once per edge; edges need two distinct fragments. Start with 41-base
words (shorter if needed), then increase by 20 bases and finally try the read
length minus ten, capped at 151 bases. Longer words can resolve repeats that
collapse at shorter word lengths. Stop at the first unique supported product.
Only a single covered, acyclic path between terminal 21-base reference seeds
is accepted. Interior flanks can contain sample SNPs or indels. Every
base in the reconstructed path comes from read-supported edges. No repeat
sequence or copy number is filled in from a candidate template.

The pool must contain at least three fragments (or the larger configured
minimum). Pools above 2,000 fragments first retain up to two original molecules
per identical sequence-pair pattern, selected deterministically by quality and
molecule ID. Both mates stay together. This reduction preserves every graph edge
supported by at least two fragments, including minority branches. Original
fragments still supply terminal-anchor evidence and support counts. If the
sequence diversity still needs more than 2,000 representatives, assembly reports
`fragment_diversity_limit` instead of discarding competing sequences.
`short_read_diagnostics.tsv` records the input/selected fragment counts and the
selection method. Limits of 50,000 edges and eight flank backgrounds also bound
work. Short anchors, gaps, competing branches, repeat cycles, or exceeded limits
yield no assembled sequence. Supported contigs pass through the assembly
workflow's Sassy PCR search (two primer edits), `pcr_rows_to_products`, and
`legacy_assembly_call_rows`, sharing product selection, primer-indel calibration
and rounding. Original fragment IDs supply support; an assembled product
is not reported as multiple independently observed full-spanning reads.

If assembly fails and no independent sequence measurement resolves the count,
usable mapping evidence can still supply a diagnostic prediction, with the
accepted count left blank. Loci without
usable sizing evidence receive explicit prior-only predictions with zero
sample-evidence confidence, provided a locus-specific prior is available.

Orphan mates contribute through flank recruitment even when their surviving
mate was originally R2. Their separate FASTQ file is mapped as unpaired reads
with the original molecule identities retained.

The current database schema is still required. Minimap2 output is streamed
through htslib rather than materialized as text SAM. `--keep-intermediates`
retains filtered reads and `candidate_mapping/candidate_alignments.bam`.

## Matching assembly fingerprints

Directly observed primer-bounded products use the assembly/MLVA_finder product
length convention, including configured primer lengths when a primer match
contains an insertion or deletion. The final allele states include measured
half alleles and counts outside the mapping candidate range. Exact ordered repeat
flanks can also measure expansions, contractions and zero-repeat alleles
outside the panel's expected range. Fuzzy flank matches retain the plausibility
check to avoid mistaking partial repeat sequence for a novel allele.

`--assembly-round-tolerance` applies to direct FASTQ measurements as well as
assembly calls. Its default of `0.25` matches MLVA_finder `-r 0.25`: raw values
exactly at `.25` and `.75` become half alleles. Set it to `0` to retain
unrounded values. Use the same panel calibration and rounding setting when
comparing `mlva_fingerprint.tsv` between modes.

The regression check uses matched synthetic FASTA and FASTQ molecules and
compares their final fingerprints against MLVA_finder-verified expected calls:

```bash
python -m pytest -q tests/test_fastq_assembly_concordance.py
```

This checks single reads, paired reads, overlapping pairs and fragmented
paired reads, including partial alleles and primer indels. It does not establish
concordance for a real dataset. On matched culture samples, compare both final
fingerprints using the same sample ID, locus panel and rounding option. For
example, from their common parent directory:

```bash
diff -u assembly/mlva_fingerprint.tsv fastq/mlva_fingerprint.tsv
```

Inspect discordant loci in `assembly_amplicons.tsv`, `common_locus_calls.tsv`
and `molecule_candidate_evidence.tsv`. Reads that do not span the informative
interval, collapsed or broken assembly repeats, and multiple alleles can still
produce different calls. A mixed sample's read-supported allele need not equal
MLVA_finder's selected assembly product. Custom MLVA_finder binning corrections
are not applied by this FASTQ caller.

The former short-read path measured mates only after candidate mapping retained
them. Unmapped mates were discarded, candidate-relative spans could favor short
templates, and repeat-only competitors could outweigh direct evidence. Direct
measurements already extended the shared allele state set beyond candidate
bounds; final rounding was not the source of this compression. Recovered direct
and assembled counts continue to extend that state set on demand. Candidate
libraries are not enlarged speculatively when partial reads cannot distinguish
longer alleles.

`short_read_diagnostics.tsv` separates recruitment, direct measurements, the
original candidate estimate, graph outcome, and final count. It includes rescued
mate, boundary, clipping, unmapped-mate and unique-fragment counts, plus the call
method. `--keep-intermediates` also retains successful products in
`short_read_microassembly.fasta`. Primary table columns and report layout are
unchanged.

## Short-read reconstruction regression checks

The new synthetic reconstruction cases compare final short-read fingerprints
with assembly calls on the identical underlying sequence. The database contains
allele 4, and none of the input fragments independently measures the complete
sample repeat interval. Before values use the three short-read inference modules
from `12fb524`; after values use the revised reconstruction.

| Synthetic case | Assembly allele | Before | After |
| --- | ---: | --- | --- |
| Single-end tiling over a perfect repeat | 5 | 5, estimated | 5, reconstructed |
| Only R1 survives paired-read QC | 5 | Not found | 5, reconstructed |
| Only R2 survives paired-read QC | 5 | Not found | 5, reconstructed |
| Paired reads, imperfect repeat and a sample flank SNP | 12 | 4, estimated | 12, reconstructed |

All four reconstructed sequences exactly match their input assembly sequence;
the final calls are `PASS`. Separate checks retain rejection of unresolved
repeat cycles, competing haplotypes and inadequate coverage. Run:

```bash
python -m pytest -q tests/test_short_read_fragments.py tests/test_fastq_assembly_concordance.py
```

These cases demonstrate specific recovery improvements. They do not measure
agreement with cgMLST or Parsnp, real-sample accuracy, or whole-genome throughput;
those require matched data from the server.

## Historical runtime checks

These measurements predate observed-reference FASTQ inference and do not
benchmark the current algorithm. Recheck runtime and accuracy on server data.

The paired-fragment changes were compared against `5b79e47` on 2026-10-02,
using two threads and three alternating fresh-process trials. No database
classification was included. Counts below were compared with assembly mode on
the identical synthetic products, without supplying assembly calls to inference.

| Workload | Exact agreement, before → after | Median seconds, before → after | Runtime change | Median peak RSS MiB, before → after | Graph attempts |
| --- | ---: | ---: | ---: | ---: | ---: |
| 11 easy loci, 1,160 pairs | 11/11 → 11/11 | 10.08 → 2.28 | −77.4% | 72.06 → 73.31 | 0/11 |
| Six difficult imperfect-repeat loci, 660 pairs | 0/6 → 6/6 | 2.50 → 0.46 | −81.5% | 76.23 → 75.25 | 5/6 |
| Easy loci plus 50,000 random background pairs | 11/11 → 11/11 | 11.78 → 4.06 | −65.5% | 73.05 → 74.86 | 0/11 |

On the difficult set, within-one agreement improved from 1/6 to 6/6,
within-two from 2/6 to 6/6, MAE from 8.67 to 0, median absolute error from
8.5 to 0, and Spearman correlation from 0.131 to 1. Counts were
5, 8, 12, 16, 20 and 24, with a candidate range ending at 6. Neither version
had missing calls on these measured datasets. Easy-set MAE stayed zero.

[Full trials, input hashes, generation recipes, and per-locus metrics](../reference/short-read-fragments-benchmark.json)
are retained. RSS is the larger of parent and largest child, not concurrent
total. These small synthetic workloads do not establish real-isolate accuracy,
whole-genome throughput, off-target homology behavior, or how often assembly
will be needed in routine samples. The matched real dataset remains on the
user's server. The [comparison and timing utility](../../scripts/README.md)
can run there unchanged.

The regression checks include recruitment/mate rescue, duplicate and ambiguous
seeds, repeat-only competitors, unique versus periodic pair overlaps, graph
coverage, branching/cycles, and complete paired inference through 24 repeats.
The existing full suite also exercises assembly, long-read inference, profiles,
classification, and reporting. Run `python -m pytest -q` with native tools on
`PATH`; tests requiring unavailable tools are otherwise skipped.

### Earlier runtime baseline

Anchor measurement rejects impossible matches before starting Sassy: a window
must be long enough, and a concrete primer with at most `k` edits must retain
at least one of `k+1` disjoint exact pieces. Possible matches still use the
existing native matcher. Short-read extraction also reuses up to 512 compact
measurement results per call, keyed by locus, sequence and quality. Independent
molecules retain their individual evidence and support counts.

A local benchmark on 2026-09-30 used 1,160 synthetic paired 100 bp reads from
11 loci, two threads and three fresh-process trials per version:

| Implementation | Median runtime | Median peak RSS |
| --- | ---: | ---: |
| Restored `16f786f8` | 40.67 s | 70.44 MiB |
| Concordance fixes with runtime optimizations | 12.16 s | 69.66 MiB |

All 18 TSV/CSV outputs matched the concordance implementation before runtime
optimization in every trial. This small workload emphasizes anchor measurement;
whole-genome I/O, database classification and server hardware can change the
gain. [Raw timings, input hashes and measurement details](../reference/short-read-runtime.json)
record the benchmark's scope. Peak RSS describes the larger of the Python
process and its largest native child, rather than combined concurrent memory.

The regression check verifies rejection of impossible searches, retention of
valid edit matches, and measurement reuse without merging molecule support:

```bash
python -m pytest -q tests/test_short_read_runtime.py
```

## Automatic taxon identification

With `--database`, classification reuses the original molecules’ alignments
against observed reference amplicons. Sequence mismatches and resolvable
repeat-length differences contribute to the shared alignment-likelihood model.
Isolate mode combines evidence for one reference source; metagenome mode uses
EM to estimate reference-component support. Taxon labels annotate reference
groups without pooling their support.

Results are written under `classification/`, appended to `profile_matches.tsv`,
and shown in `report.html`. Sequence evidence can remain informative when a
repeat count is unresolved. Low-support and indistinguishable references remain
explicit. See [mapping classification](../concepts/mapping-classification.md).

## Measured and presence-only evidence

Each assembled product, merged pair, and original read uses the same panel
anchors and assembly-calibrated repeat convention. Direct counts require an
independent product or flank measurement. Both boundaries on separate unmerged
mates do not establish a measured repeat interval. Adequately flank-supported
spanning CIGAR evidence or calibrated paired lengths can provide a diagnostic
prediction. Conflicting observations and partial reference matches retain a best
numeric estimate with uncertainty. Exact terminal motif runs with a unique flank
condition a locus prior on a lower bound. These distinctions appear in `genotype_evidence.tsv`.

`common_locus_calls.tsv` separates all detected molecules (`presence_support`)
from informative molecules (`molecule_support`). The compact call table uses
`PRESENT_COUNT_UNKNOWN` for detected loci without identifiable sample sizing.
Diagnostic predictions remain available when evidence or an eligible prior
supports them; withheld accepted counts stay blank. Prior-only predictions have
zero sample-evidence confidence. Presence is never inferred from a prediction.
Reported intervals may describe
bounded candidate envelopes rather than measured spans.

## Mixtures and confidence

Allele support counts only molecules with discriminating boundary evidence.
Repeat-internal locus reads remain in `uninformative_locus_reads`. Multiple
defensible alleles are preserved with primary/secondary support, informative
molecules, fractions, and `mixture_status`. Fractions are left empty when no
allele-discriminating molecule exists.

Acceptance reasons are textual and auditable. A primer/flank-bounded local
assembly with adequate molecule support can supply an accepted count. Direct
calls without an agreeing supported consensus also need enough independently
agreeing molecules and a selected-allele fraction meeting the existing minimum
probability threshold. Insufficient support and unresolved conflicts withhold
the count. The HTML report separates diagnostic predictions from accepted calls;
neither model support nor acceptance establishes calibrated accuracy.

## Metadata and MYOGA

Metadata can use `sample_id`, `run_accession`, `sra_run`, or `accession` as its
join key. Aliases for BioSample, dates, coordinates, location, country, host,
isolation source, and study are normalized without removing original columns.

`myoga_samples.csv` uses `genome_id = sample_id`. MYOGA recognizes
`genome_id`, `latitude`, `longitude`, `location`, and `collection_date`
directly. If a Newick tree is generated from the same sample IDs, its tip names
must remain exactly those values. In MYOGA, load the Newick/Parsnp tree and
`myoga_samples.csv`, then select `genome_id` as the ID column if it is not
selected automatically. `myoga_loci.csv` is a long-form companion for external
filtering; MYOGA does not require it.

## Manifest batches and HPC

```text
sample_id\treads1\treads2\tmetadata_id
SRR000001\t/path/SRR000001_1.fastq.gz\t/path/SRR000001_2.fastq.gz\tSAMN000001
SRR000002\t/path/SRR000002.fastq.gz\t.\tSAMN000002
```

```bash
mlvamaps call -p panel.tsv -i sr --manifest samples.tsv \
  --sample-metadata metadata.tsv \
  --profiles profiles.tsv --database reference_build \
  -o results -t 32
```

Samples run independently and write to `results/<sample_id>/`. One malformed or
missing sample is recorded in `results/batch_summary/batch_status.tsv` without
stopping the rest.
Successful sample directories resume by default; use a fresh output directory or
`--force` to recompute calls under the current acceptance policy. Re-exporting old
results alone does not apply it.
Combined standard, sample-summary, and MYOGA tables are written at the batch
root's clearly scoped `batch_summary/` directory.

Independent samples run concurrently under one global `--threads` budget. The
number of active samples is the minimum of sample count, available threads, and
a memory-aware cap of four. Each sample receives an equal integer share of the
budget, and the allocated total never exceeds `--threads`. Set a lower cap for
large samples or memory-constrained nodes:

```bash
export MLVAMAPS_MAX_CONCURRENT_SAMPLES=2
mlvamaps call -p panel.tsv -i short_read_directory/ --short-reads \
  --database reference_build -o results -t 32
```

Completion order does not change combined output order. Within each sample,
FASTQ/QC streams in bounded chunks and minimap2 performs the allocated
multithreaded alignment. Progress reports sample/worker allocation and stage
transitions unless `--quiet` is selected.

For Slurm arrays, split the manifest by row while preserving its header and run
one manifest shard per task into separate output roots. Merge the resulting TSV
or CSV files after all tasks complete. mlvamaps does not submit scheduler jobs.

## Synthetic worked example

Generate the tiny offline example:

```bash
python examples/make_illumina_example.py examples/illumina_demo
```

With a completed reference build for this synthetic panel at `reference_build`, run:

```bash
mlvamaps call --database reference_build -i sr \
  --fq1 examples/illumina_demo/SRR_DEMO_1.fastq.gz \
  --fq2 examples/illumina_demo/SRR_DEMO_2.fastq.gz \
  --sample-id SRR_DEMO \
  --sample-metadata examples/illumina_demo/metadata.tsv \
  -o examples/illumina_demo/results
```

Open `results/report.html`, inspect the exact-versus-unresolved evidence, and
load `results/myoga_samples.csv` into MYOGA. An assembly can be called
independently for scientific concordance assessment:

```bash
mlvamaps call -p examples/illumina_demo/panel.tsv \
  -i examples/illumina_demo/truth.fasta.gz -o examples/illumina_demo/truth
```

## Troubleshooting

- **Different FASTQ counts or IDs:** regenerate mates together; do not sort one
  file independently.
- **Many ambiguous pairs:** provide a reference build with longer, divergent
  locus flanks and review similar loci in the panel.
- **Presence-only locus:** this is expected when neither reads nor the local
  graph resolve both boundaries. Do not replace the blank call with the
  expected-range midpoint.
- **Database predates the context schema:** rebuild it with the current
  `build-reference` command, or omit `--database` and provide a rich panel.
- **MYOGA row does not attach to a tip:** make `genome_id` exactly equal to the
  Newick label, including suffixes and case.
