# FASTQ and amplicon sequencing workflow

The FASTQ path accepts `.fastq`, `.fq`, and gzip-compressed equivalents. It
is intended for current high-accuracy long-read sequencing. Reads may span a
complete primer product or cover only part of a locus. The default
minimum mean Q score is 17, approximately 98% per-base accuracy.

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

Lower `--min-qscore` explicitly when working with older or noisier base calls.

This stage returns:

- `qc_summary.tsv`
- `filtered_reads.fastq.gz`

## 3. Recruit reads competitively to loci

The default `--fastq-strategy recruit` builds a competitive reference bank for
every panel locus and allowed repeat allele, then maps all retained reads
against that bank in one minimap2 operation. Complete products from
`--recruitment-database` (or `--database`) are preferred. The dedicated option
does not activate phylogenetic placement. When no database product exists, a
rich panel can generate a recorded synthetic product from its primers, flanks,
motif, and repeat range.

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
- `local_locus_products.fasta`
- Native references and alignments under `recruitment/`

Presence statuses are `PRESENT_GENOTYPED`, `PRESENT_PROVISIONAL`,
`PRESENT_UNTYPED`, and `NO_EVIDENCE`. Presence is intentionally independent
from whether a repeat count can be reported.

## 4. Pair primers and apply the specificity fallback

MLVAMaps creates a lossless FASTA projection of retained reads for Amplirust.
Amplirust handles IUPAC primer bases, both orientations, primer alignment, and
valid product-length constraints. The original FASTQ qualities remain connected
to the assigned read. A valid primer-pair assignment takes precedence when it
and recruitment both recover the same read. Primer pairing also supplies a
fallback when a panel lacks a usable complete recruitment product.

This stage returns:

- `filtered_reads.fasta`
- `read_locus_assignments.tsv`
- Native evidence under `amplirust/`

Use `--fastq-strategy primer` to disable recruitment and retain the historical
primer-only workflow.

## 5. Extract repeat evidence and local products

Accepted reads are oriented to the forward-primer direction. MLVAMaps defines
the inner primer-to-primer region and refines it with optional locus flanks.
For each usable read it calculates:

- Repeat-region coordinates and length.
- Raw and nearest-integer repeat count.
- Repeat pattern and motif k-mer count.
- Motif mismatches.
- Primer, flank, and mean read-quality scores.

This stage returns `read_repeat_features.tsv`.

For recruited full products, reads of the modal product length are combined
into a local majority product in `local_locus_products.fasta`. Repeat-spanning
partial reads can produce provisional allele evidence when no full product is
available. Presence-only reads never become allele calls.

## 6. Dereplicate and cluster with VSEARCH

Each locus is processed independently. VSEARCH performs exact amplicon
dereplication followed by abundance-sorted global clustering. The identity
definition includes gaps. Low-complexity masking is disabled because VNTR
sequences are expected to be repetitive, and sensitive short-word seeding is
used so indels do not disappear behind long-word requirements.

Defaults:

- Global identity: `0.97`
- Minimum retained cluster size: `1`

Controls:

- `--cluster-min-identity`
- `--min-cluster-size`
- `--vsearch-bin`

This stage returns `vntr_asv_table.tsv` and diagnostic files under
`vsearch/`.

## 7. Preserve observed representatives

The VSEARCH centroid is an actual observed read. MLVAMaps never replaces it
with a generated consensus. Parasail globally aligns each unique cluster repeat
sequence to its selected observed representative and annotates substitutions,
insertions, deletions, and edit distance without clipping sequence ends.

This stage returns:

- `vntr_asv_representatives.fasta`
- `vntr_asv_memberships.tsv`

The FASTA contains representative repeat regions. The membership table retains
the raw and aligned sequence evidence for every read in a retained cluster.

## 8. Estimate variant mixture abundance

MLVAMaps fits an Emu-inspired expectation-maximization model to retained
VSEARCH count evidence. Pairwise representative similarity supplies the
assignment likelihoods, while the abundance estimate from each iteration
becomes the prior for the next iteration. This separates meaningful secondary
variants from trace clusters and estimates the fraction of each component.

Secondary variants require both the configured abundance fraction and
`--min-secondary-reads` (default 2) to become confirmed. Singleton secondaries
remain visible as `CANDIDATE` evidence but do not alter the primary signature.

This stage returns `vntr_mixture_abundance.tsv`. Control the meaningful/trace
boundary with `--min-mixture-fraction`.

## 9. Map to dominant representatives

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

## 10. Predict read alleles

Reads in the EM-dominant cluster contribute primary repeat-count likelihoods.
Representative-relative edits reduce evidence weight so clean observations
contribute more than heavily edited sequences. Concordant likelihoods are
combined multiplicatively, allowing confidence to increase with support.
`--max-confidence-depth` caps the effective evidence at 25 by default to limit
overconfidence from correlated reads.

By default, primer-spanning reads use the same calibrated product-length
conversion and historical rounding convention as assembly calls. This makes
the FASTQ fingerprint comparable to the fingerprint obtained after assembling
the same sample. The unrounded measurement remains in the read-level output.
Use `--read-calling-convention probabilistic` to retain direct half-unit
inference instead.

Singleton clusters are retained by default (`--min-cluster-size 1`). A single
spanning read can therefore contribute a provisional allele and remain in the
fingerprint; `--min-depth` still marks insufficient support as `LOW_DEPTH`.

This stage returns `read_level_allele_predictions.tsv`.

## 11. Call each locus

The Bayesian caller combines primary-cluster read probabilities and reports
the best and second-best repeat count, posterior values, total and primary
depth, retained candidate and confirmed variant counts, dominant fraction, and
secondary alleles. Secondary variants are interpreted independently and never
averaged into the primary allele posterior.

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

MLVAMaps converts the locus calls to wide and probabilistic fingerprints. If a
profile database is present, it ranks known rows using the full allele
probability distributions at observed loci, while also reporting conventional
repeat-count distance, missing-locus-aware comparison counts, and confidence.

This stage returns:

- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `report.html`

The report leads with panel completeness, locus-specific review findings, and
the closest profile/reference interpretation. It only shows mixture and mapping
sections when those data exist. Exact allele, mixture, mapping, SNP, and
distance-component tables remain available in collapsed detail sections.
