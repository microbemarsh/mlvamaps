# Competitive FASTQ workflows

The descriptions below apply only to `--lr-engine competitive` and
`--sr-engine competitive`, the short-read default. Long reads default to
`--lr-engine spanning`; short reads can opt into `--sr-engine repeat-likelihood`.
Those engines use [locus reconstruction](../concepts/locus-reconstruction.md).

# FASTQ and amplicon sequencing workflow

Both long reads and [Illumina reads](illumina.md) use competitive minimap2
alignment against the same candidate MLVA allele contexts. The technologies
differ in evidence extraction, not in the meaning of the final locus call.

The FASTQ path accepts `.fastq`, `.fq`, and gzip-compressed equivalents. It
is intended for current high-accuracy long-read sequencing. Reads may span a
complete primer product or cover only part of a locus. The default
minimum mean Q score is 15, approximately 97% per-base accuracy.

```bash
mlvamaps call -p primers.tsv -i sample.fastq.gz -o results
```

## 1. Load the panel

mlvamaps reads a minimal primer table or a richer locus table. Optional profile
rows are loaded with `--profiles`.

The panel determines locus names, primer sequences, valid amplicon lengths,
repeat interpretation, and fingerprint column order.

## 2. Filter reads

For metagenomic input, an optional target-taxon screen can run before ordinary
read QC:

```bash
mlvamaps call -p primers.tsv -i metagenome.fastq.gz \
  --taxon-screen-index target_taxon.idx
```

The index is supplied rather than built by mlvamaps. See
[bede/deacon-indexes](https://github.com/bede/deacon-indexes) for
index-building information.

Deacon performs native SIMD minimizer matching and retains target-like reads.
mlvamaps records the original, retained, and rejected read totals in
`qc_summary.tsv`; the retained FASTQ and full native summary are written under
`taxon_screen/`. Screening is bypassed when no index is supplied, so pure
culture FASTQ behavior is unchanged.

Reads are parsed with their original qualities and filtered by:

- `--min-read-length`
- `--max-read-length`
- `--min-qscore`

Lower `--min-qscore` explicitly when working with older or noisier base calls.

This stage returns:

- `qc_summary.tsv`
- `filtered_reads.fastq.gz`

## 3. Recruit reads competitively to loci

The long-read algorithm maps retained reads against all panel loci and allowed
repeat alleles in one competitive minimap2 operation. With `--database`, it
reuses the completed database's candidate FASTA, metadata, provenance, and
`long.mmi` index. Candidate resources are not regenerated per sample. Complete
products from `--recruitment-database` (or `--database`) are preferred for the
additional recruitment path; the dedicated option does not activate
reference classification. Without a database, a rich panel can generate a
recorded, bounded synthetic candidate bank from primers, flanks, motif, and
repeat range.

minimap2 SAM is streamed to htslib instead of being written as a text
intermediate. htslib decodes flags, CIGAR operations, and required tags and
writes compressed BAM. The BAM is deleted in normal operation and retained as
`candidate_mapping/candidate_alignments.bam` only when the applicable
intermediate/debug option requests alignments.

Mapped evidence is separated into:

- `FULL_PRODUCT`: covers the complete product and can enter the normal
  assembly-equivalent calling path.
- `REPEAT_INFORMATIVE`: spans both repeat boundaries but not the complete
  product; supplies a provisional candidate allele.
- `PRESENCE_ONLY`: maps specifically to the locus but cannot measure the
  complete repeat.

This stage returns:

- `locus_recruited_reads.tsv`
- `locus_presence.tsv`
- `local_locus_products.fasta.gz`
- Native references and alignments under `recruitment/`

Presence statuses are `PRESENT_GENOTYPED`, `PRESENT_PROVISIONAL`,
`PRESENT_UNTYPED`, and `NO_EVIDENCE`. Presence is intentionally independent
from whether a repeat count can be reported.

## 4. Pair primers and apply the specificity fallback

`mlvamaps` creates a lossless FASTA projection of retained reads and passes it
through the same Sassy-backed in silico PCR engine used for assemblies. Sassy
discovers approximate primer locations; `mlvamaps` handles deterministic IUPAC
primer expansion, orientation, primer pairing, and product-length constraints.
The original FASTQ qualities remain connected to the assigned read. A valid
primer-pair assignment takes precedence when it and recruitment both recover
the same read. Primer pairing also supplies a fallback when a panel lacks a
usable complete recruitment product.

This stage returns:

- `filtered_reads.fasta.gz`
- `read_locus_assignments.tsv`
- Sassy-backed evidence under `in_silico_pcr/`

## 5. Extract repeat evidence and local products

Accepted reads are oriented to the forward-primer direction. mlvamaps defines
the inner primer-to-primer region and refines it with optional locus flanks.
For each usable read it calculates:

- Repeat-region coordinates and length.
- Raw and nearest-integer repeat count.
- Repeat pattern and motif k-mer count.
- Motif mismatches.
- Primer, flank, and mean read-quality scores.

This stage returns `read_repeat_features.tsv`.

Complete products in the EM-dominant cluster are globally aligned as a partial
order alignment with the `spoars` Python bindings. Its consensus is written to
`local_locus_products.fasta.gz`, then processed by the same Sassy-backed in
silico PCR,
legacy primer-coordinate product-size calculation, and repeat caller used for
whole assemblies. This corrects minority read insertions and deletions before
repeat counting. Repeat-spanning partial reads can produce provisional allele
evidence only when no complete dominant product can be assembled and resolved
by PCR. Presence-only reads never become allele calls.

`local_assembly_concordance.tsv` records the observed read-length range and
mode, POA consensus length, PCR-derived assembly product size, raw repeat
measurement, final allele, support depth, and fallback state for every
assembled locus.

## 6. Build mapping-derived repeat groups

The competitive recruitment mapping assigns each informative long read to a
locus and candidate repeat-product class. These mappings define the groups
used for abundance and dominant-product selection; sequence-ASV clustering is
not run.
Primer-only compatibility mode keeps all complete products at a locus together
so SPOARS—not an individual read measurement—defines the allele.

This stage returns:

- `mapped_variant_table.tsv`
- `mapped_variant_representatives.fasta.gz`
- `mapped_read_memberships.tsv`

The membership table retains the mapping group plus Parasail-derived
substitution and indel diagnostics for every complete product.

## 8. Estimate variant mixture abundance

mlvamaps fits an Emu-inspired expectation-maximization model to mapped-product
counts. Pairwise group-representative similarity supplies the
assignment likelihoods, while the abundance estimate from each iteration
becomes the prior for the next iteration. This separates meaningful secondary
variants from trace clusters and estimates the fraction of each component.

Secondary variants require both the configured abundance fraction and
`--min-secondary-reads` (default 2) to become confirmed. Singleton secondaries
remain visible as `CANDIDATE` evidence but do not alter the primary signature.

This stage returns `vntr_mixture_abundance.tsv`. Control the meaningful/trace
boundary with `--min-mixture-fraction`.

## 9. Map to dominant POA products

The assembly-PCR-resolved SPOARS product at each locus supplies the minimap2
reference. All usable reads assigned to that locus are mapped back to the POA
product. Raw reads and mapping groups therefore annotate confidence and
mixtures but do not replace the primary assembly-derived sequence.

This stage returns:

- `locus_mapping_references.fasta.gz`
- `locus_read_alignments.sam`
- `locus_mapping_summary.tsv`
- `locus_snps.tsv`
- Native mapping inputs under `minimap2/locus_mapping/`

See [representative mapping and SNP evidence](../concepts/representative-mapping.md)
for thresholds and interpretation.

## 10. Predict read alleles

Reads in the EM-dominant cluster contribute primary repeat-count confidence.
Representative-relative edits reduce evidence weight so clean observations
contribute more than heavily edited sequences. Concordant likelihoods are
combined multiplicatively, allowing confidence to increase with support.
`--max-confidence-depth` caps the effective evidence at 25 by default to limit
overconfidence from correlated reads.

Each sufficiently informative long read contributes an explicit molecule-level
allele observation before consensus generation. Complete primer-bounded reads
use the same calibrated product-length conversion as assemblies; complete VNTR
spans, repeat-sized indels, one-boundary reads, and flank-only reads occupy
successively weaker evidence tiers. SPOARS supplies corrected representative
sequences, confirmatory measurements, and sequence-variant characterization,
but is not a prerequisite for a repeat-count call.

Singleton clusters are retained by default (`--min-cluster-size 1` and
`--min-depth 1`). A single spanning read can therefore contribute a
provisional allele and remain in the fingerprint. The
`SINGLE_MOLECULE_PROVISIONAL` evidence label preserves the distinction between
detection and replicated support.

This stage returns `read_level_allele_predictions.tsv`.

## 11. Call each locus

The caller fixes the primary allele from the dominant complete local product,
then combines primary-cluster read probabilities to report its confidence and
the alternatives. It also reports total and primary depth, retained candidate
and confirmed variant counts, dominant fraction, and secondary alleles.
Secondary variants are interpreted independently and never averaged into the
primary allele posterior.

This stage returns:

- `allele_calls.tsv`
- `calls.tsv`

Statuses include `PASS`, `LOW_DEPTH`, `AMBIGUOUS`, `OUT_OF_RANGE`,
`MULTIPLE_VARIANTS`, and `LOCUS_DROPOUT`.

The default is `--sample-mode metagenome`, where any meaningful secondary allele causes
`MULTIPLE_VARIANTS`, even when one allele exceeds 80%, so the dominant
per-locus signature is not mistaken for an unqualified single-strain result.
Use `--sample-mode isolate` explicitly for cultured material. Both modes keep
the assembly-equivalent dominant allele, posterior probability, and dominant
variant fraction, allowing a rapid metagenomic detection to be compared with a
later cultured assembly without erasing the original uncertainty.
Alleles at loci that are not linked by the same reads cannot be phased into
organism-specific metagenomic signatures.

## 12. Fingerprint, profiles, and report

mlvamaps converts the locus calls to wide and probabilistic fingerprints. If a
profile database is present, conventional repeat-count distance and
matched-locus count are the primary ranking keys in both FASTQ and assembly
modes. FASTQ allele probabilities break otherwise equal matches. The output
also reports missing-locus-aware comparison counts and confidence.

This stage returns:

- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `report.html`

The report leads with panel completeness, locus-specific review findings, and
the closest profile/reference interpretation. It only shows mixture and mapping
sections when those data exist. Exact allele, mixture, mapping, SNP, and
distance-component tables remain available in collapsed detail sections.


# Illumina short-read workflow

Illumina mode uses the shared FASTQ architecture. Pair identity is retained
through competitive minimap2 mapping against candidate MLVA allele contexts.
Illumina-specific molecule, overlap, boundary, repeat-indel, and pair-geometry
evidence then enters the same locus-level inference used for long reads.

## Command line

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
mlvamaps call -p panel.tsv -i short_read_directory/ --short-reads \
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
counts. When enough exact opposite-orientation mappings exist, it also records
the empirical fragment-span median, median absolute deviation, and pair count.

## Competitive locus-context mapping

Contexts come from one of two explicit sources:

1. the versioned `competitive_mapping/candidate_metadata.tsv` and
   `competitive_mapping/candidate_contexts.fasta` files in a
   current `--database`; or
2. complete products synthesized from a rich panel's primers, flanks, motif,
   and expected repeat range when no database is supplied.

A primer-only panel cannot define the repeat boundaries required by this
algorithm and is rejected with guidance to build a reference database or enrich
the panel.

Filtered mates are mapped once with minimap2 against all candidate MLVA contexts,
without an early taxon restriction. Contexts retain locus, reference, taxon,
repeat interval, expected allele, and flank provenance. Equivalent alignments
remain available to inference rather than being forced to one reference.

Mate scores are resolved together. A confident mate can rescue its unaligned
mate, equal best scores across loci are ambiguous, and confident mates assigned
to different loci are discordant. Neither category is counted as unique support
for multiple loci. Conventional MAPQ contributes to locus assignment but does
not define allele confidence among intentionally similar repeat states.

When database resources are available, the candidate FASTA, metadata,
provenance, and `short.mmi` index are reused directly. They are not copied or
regenerated for every sample. Without a database, locally synthesized resources
are written under the sample's `candidate_mapping/` directory.

## Direct VNTR inference

Candidate alleles vary only in whole repeat units. Evidence includes unique
flank mappings, boundary junctions, full VNTR spans, opposite-flank proper pairs,
forced. Candidate competition, not the primary alignment label alone, supplies
repeat-state evidence; no per-locus de novo assembly is required.

When `--database` is supplied, it must contain the versioned
`competitive_mapping/candidate_metadata.tsv`, `candidate_contexts.fasta`, and
technology-specific minimap2 indexes produced by the
current reference builder. Older databases must be rebuilt rather than being
silently reinterpreted. minimap2 output is streamed through htslib rather than
materialized as text SAM. `--keep-intermediates` retains the filtered reads and
compressed `candidate_mapping/candidate_alignments.bam`; with a database, the
candidate bank and index continue to reside in the database.

## Automatic taxon identification

With `--database`, the original retained molecules are aligned competitively
against observed reference amplicons. Sequence mismatches and resolvable
repeat-length differences contribute to the shared alignment-likelihood model.
Isolate mode combines evidence for one reference source; metagenome mode uses
EM to estimate reference-component support. Taxon labels annotate reference
groups without pooling their support.

Results are written under `classification/`, appended to `profile_matches.tsv`,
and shown in `report.html`. Sequence evidence can remain informative when a
repeat count is unresolved. Low-support and indistinguishable references remain
explicit. See [mapping classification](../concepts/mapping-classification.md).

## Exact, interval, and presence evidence

Each contig, merged pair, and original read is evaluated with the same panel
anchors and assembly-calibrated repeat convention used elsewhere in mlvamaps.
Evidence is classified as:

- `COMPLETE_ASSEMBLED_PRODUCT`
- `BOUNDARY_SPANNING_READ_PAIR`
- `BOUNDARY_SPANNING_SINGLE_READ`
- `PARTIAL_REPEAT_EVIDENCE`
- `PRESENCE_ONLY`
- `AMBIGUOUS_ASSEMBLY`
- `MULTIPLE_ALLELES`
- `LOW_DEPTH`
- `NOT_FOUND`

An exact `repeat_count` requires a contig, merged pair, or original read that
directly resolves both repeat boundaries. The two boundaries may be the rich
panel flanks or, when flanks are absent, the product primers. A read inside the
repeat or covering one boundary cannot create an exact value.

Opposite boundaries on separate mates may produce `repeat_count_min` and
`repeat_count_max`. Empirical insert size can narrow that interval when at
least two concordant spans support it. The midpoint is never copied into
`repeat_count`. Without an adequate insert estimate, the panel's expected range
is retained and the reason explains why.

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

Then run:

```bash
mlvamaps call -p examples/illumina_demo/panel.tsv -i sr \
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
