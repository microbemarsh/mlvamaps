# FASTQ and amplicon sequencing workflow

The FASTQ path accepts `.fastq`, `.fq`, and gzip-compressed equivalents. It
is intended for current high-accuracy long-read sequencing. Reads may span a
complete primer product or cover only part of a locus. The default
minimum mean Q score is 15, approximately 97% per-base accuracy.

```bash
mlvamaps call primers.tsv sample.fastq.gz -o results
```

## 1. Load the panel

MLVAMaps reads a minimal primer table or a richer locus table. Optional profile
rows are loaded with `--profiles`.

The panel determines locus names, primer sequences, valid amplicon lengths,
repeat interpretation, and fingerprint column order.

## 2. Filter reads

For metagenomic input, an optional target-taxon screen can run before ordinary
read QC:

```bash
mlvamaps call primers.tsv metagenome.fastq.gz \
  --taxon-screen-index target_taxon.idx
```

The index is supplied rather than built by MLVAMaps. See
[bede/deacon-indexes](https://github.com/bede/deacon-indexes) for
index-building information.

Deacon performs native SIMD minimizer matching and retains target-like reads.
MLVAMaps records the original, retained, and rejected read totals in
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

Complete products in the EM-dominant cluster are globally aligned as a partial
order alignment with the `spoars` Python bindings. Its consensus is written to
`local_locus_products.fasta`, then processed by the same Sassy in-silico PCR,
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
- `mapped_variant_representatives.fasta`
- `mapped_read_memberships.tsv`

The membership table retains the mapping group plus Parasail-derived
substitution and indel diagnostics for every complete product.

## 8. Estimate variant mixture abundance

MLVAMaps fits an Emu-inspired expectation-maximization model to mapped-product
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

- `locus_mapping_references.fasta`
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

By default, the dominant cluster's SPOARS consensus goes through the complete
assembly calling path. This includes primer matching, legacy-compatible product
boundary calculation, product selection, repeat conversion, and historical
rounding. Per-read likelihoods increase or decrease confidence around that
allele but cannot move the primary call to a different repeat count. The
unrounded product measurement is retained in `allele_calls.tsv` and
`calls.tsv`. Use `--read-calling-convention probabilistic` to retain direct
per-read half-unit inference instead.

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

MLVAMaps converts the locus calls to wide and probabilistic fingerprints. If a
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
