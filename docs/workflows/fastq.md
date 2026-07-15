# FASTQ and amplicon sequencing workflow

The FASTQ path accepts `.fastq`, `.fq`, and gzip-compressed equivalents. It
is intended for accurate amplicon or long-read sequencing in which individual
reads contain a valid forward-primer/reverse-primer product.

```bash
mlvamaps call primers.tsv sample.fastq.gz -o results
```

## 1. Load the panel

MLVAMaps reads a minimal primer table or a richer locus table. Optional profile
rows are loaded with `--profiles`.

The panel determines locus names, primer sequences, valid amplicon lengths,
repeat interpretation, and fingerprint column order.

## 2. Filter reads

Reads are parsed with their original qualities and filtered by:

- `--min-read-length`
- `--max-read-length`
- `--min-qscore`

This stage returns:

- `qc_summary.tsv`
- `filtered_reads.fastq.gz`

## 3. Pair primers and assign loci

MLVAMaps creates a lossless FASTA projection of retained reads for Amplirust.
Amplirust handles IUPAC primer bases, both orientations, primer alignment, and
valid product-length constraints. The original FASTQ qualities remain connected
to the assigned read.

This stage returns:

- `filtered_reads.fasta`
- `read_locus_assignments.tsv`
- Native evidence under `amplirust/`

Assignment requires a valid paired-primer product. A shotgun read covering only
part of a target is not called by this workflow.

## 4. Extract repeat evidence

Accepted reads are oriented to the forward-primer direction. MLVAMaps defines
the inner primer-to-primer region and refines it with optional locus flanks.
For each usable read it calculates:

- Repeat-region coordinates and length.
- Raw and nearest-integer repeat count.
- Repeat pattern and motif k-mer count.
- Motif mismatches.
- Primer, flank, and mean read-quality scores.

This stage returns `read_repeat_features.tsv`.

## 5. Dereplicate and cluster with VSEARCH

Each locus is processed independently. VSEARCH performs exact amplicon
dereplication followed by abundance-sorted global clustering. The identity
definition includes gaps. Low-complexity masking is disabled because VNTR
sequences are expected to be repetitive, and sensitive short-word seeding is
used so indels do not disappear behind long-word requirements.

Defaults:

- Global identity: `0.97`
- Minimum retained cluster size: `2`

Controls:

- `--cluster-min-identity`
- `--min-cluster-size`
- `--vsearch-bin`

This stage returns `vntr_asv_table.tsv` and diagnostic files under
`vsearch/`.

## 6. Preserve observed representatives

The VSEARCH centroid is an actual observed read. MLVAMaps never replaces it
with a generated consensus. Parasail globally aligns each unique cluster repeat
sequence to its selected observed representative and annotates substitutions,
insertions, deletions, and edit distance without clipping sequence ends.

This stage returns:

- `vntr_asv_representatives.fasta`
- `vntr_asv_memberships.tsv`

The FASTA contains representative repeat regions. The membership table retains
the raw and aligned sequence evidence for every read in a retained cluster.

## 7. Map to dominant representatives

The most supported retained VSEARCH variant at each locus supplies its complete
observed amplicon as a minimap2 reference. All usable reads assigned to that
locus are mapped back to the representative.

This stage returns:

- `locus_mapping_references.fasta`
- `locus_read_alignments.sam`
- `locus_mapping_summary.tsv`
- `locus_snps.tsv`
- Native mapping inputs under `minimap2/locus_mapping/`

See [representative mapping and SNP evidence](../concepts/representative-mapping.md)
for thresholds and interpretation.

## 8. Predict read alleles

Every read in a retained cluster contributes a repeat-count probability.
Representative-relative edits reduce evidence weight so clean observations
contribute more than heavily edited sequences.

This stage returns `read_level_allele_predictions.tsv`.

## 9. Call each locus

The Bayesian caller combines read probabilities and reports the best and
second-best repeat count, posterior values, raw and effective depth, retained
variant count, and dominant variant.

This stage returns:

- `allele_calls.tsv`
- `calls.tsv`

Statuses include `PASS`, `LOW_DEPTH`, `AMBIGUOUS`, `OUT_OF_RANGE`,
`MULTIPLE_VARIANTS`, and `LOCUS_DROPOUT`.

## 10. Fingerprint, profiles, and report

MLVAMaps converts the locus calls to wide and probabilistic fingerprints. If a
profile database is present, it ranks known rows by repeat-count distance and
reports mismatched loci and comparison confidence. It also summarizes novelty
from distance and call uncertainty.

This stage returns:

- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `novelty_scores.tsv`
- `report.html`

The report includes call status, a generated gel comparison, representative
mapping coverage, and SNP evidence.
