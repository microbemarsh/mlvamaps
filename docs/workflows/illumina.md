# Illumina short-read workflow

Illumina mode uses the shared FASTQ architecture. Paired reads are recruited by
unique flank seeds during QC, retaining both mates independently of competitive
minimap2 mapping. Direct sequence measurements enter shared allele inference;
uncertain loci can use a bounded local graph. Single-end reads and QC-surviving
orphan mates also participate in recruitment, direct measurement and local
assembly.

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
sequence measurements or a covered two-flank graph must support sizing.

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
Partial mates retain local measurement windows and the fragment graph, avoiding
unnecessary full-product realignment.

When supported complete segments have sequencing errors in their primers, the
shared SPOARS correction step can recover the product before assembly PCR.
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
`--short-insert-mean BP --short-insert-sd BP` together. The caller can then
estimate an unresolved repeat length from complete inward-facing reads on
opposite unique flanks. Mean and standard deviation describe the full library
fragment length, including both reads. Estimates use Gaussian insert-length
likelihoods and remain low-confidence `ESTIMATED`, with method
`paired_fragment_length`. They cannot replace direct measurements or create
measured support for classification. Calibration is not inferred from
overlapping pairs or locus-only template lengths, which can be length-biased.
The Gaussian fallback does not correct for target-selection bias; targeted or
amplicon libraries may violate its fragment-sampling assumptions.

## Conditional micro-assembly

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
minimum). Limits of 2,000 fragments, 50,000 edges, and eight flank backgrounds
bound work. Short anchors, gaps, competing branches, repeat cycles, or exceeded
limits yield no assembled sequence. Supported contigs pass through the assembly
workflow's Sassy PCR search (two primer edits), `pcr_rows_to_products`, and
`legacy_assembly_call_rows`, sharing product selection, primer-indel calibration
and rounding. Original fragment IDs supply support; an assembled product
is not reported as multiple independently observed full-spanning reads.

If assembly fails and no independent sequence measurement resolves the count,
usable mapping evidence can still supply a low-confidence `ESTIMATED` count.
This remains distinct from an independently measured product. Loci without
usable sizing evidence remain unavailable.

Orphan mates contribute through flank recruitment even when their surviving
mate was originally R2. Their separate FASTQ file is not mapped, so they do not
add reference-alignment SNP evidence.

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
mates do not establish a measured repeat interval. Unresolved calls can instead
carry an `ESTIMATED` count based on conflicting observations, repeat-spanning
CIGAR evidence, or reference-conditioned alignment support. These remain low
confidence and include their estimation method, support, and uncertainty reason.

`common_locus_calls.tsv` separates all detected molecules (`presence_support`)
from informative molecules (`molecule_support`). The compact call table uses
`PRESENT_COUNT_UNKNOWN` for detected loci without either a measurement or an
estimate and retains an empty `repeat_count`. Interval fields remain empty when
no interval has been measured. The report's visible per-locus table includes
all counts and estimates, with `Unavailable` for loci lacking sizing evidence.

## Mixtures and confidence

Allele support counts only molecules with discriminating boundary evidence.
Repeat-internal locus reads remain in `uninformative_locus_reads`. Multiple
defensible alleles are preserved with primary/secondary support, informative
molecules, fractions, and `mixture_status`. Fractions are left empty when no
allele-discriminating molecule exists.

Confidence reasons are textual and auditable. A primer/flank-bounded local
assembly with adequate molecule support is high confidence; direct molecule
evidence without depth is provisional; conflicts lower confidence; interval or
presence-only rows have no falsely precise probability. The HTML report uses a
dedicated Illumina table and labels unresolved rows explicitly.

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
Successful sample directories resume by default; use `--force` to recompute.
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
