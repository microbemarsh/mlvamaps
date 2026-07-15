# Assembly workflow

The assembly path accepts FASTA assemblies and can optionally add read support.

```bash
mlvamaps call primers.tsv assembly.fasta -o results
```

## 1. Find paired-primer products

Amplirust searches both assembly orientations with degenerate-primer and
primer-error support. MLVAMaps filters circular-wrap records and products
outside the locus-specific amplicon range.

This stage returns:

- `assembly_amplicons.tsv`: locus, contig, 1-based coordinates, orientation,
  product size, and primer mismatches.
- `assembly_amplicons.fasta`: extracted products.
- Native evidence under `amplirust/`.

If several products are found for a locus, the compact call uses the product
with the fewest total primer mismatches and then the shortest product.
All accepted products remain available in the amplicon table.

## 2. Estimate repeat count

MLVAMaps converts product size to repeat count when the panel provides
repeat-unit length, nominal repeat units, and expected product size. It retains
the raw estimate and reports the nearest integer as the conventional MLVA call.

Assembly statuses:

- `PASS`: product found and repeat count could be calculated.
- `PRESENT_COUNT_UNKNOWN`: product found, but panel metadata was insufficient
  for a repeat-count calculation.
- `NOT_FOUND`: no accepted paired-primer product was found.

This stage returns `calls.tsv`.

## 3. Add optional FASTQ support

```bash
mlvamaps call primers.tsv assembly.fasta --reads sample.fastq.gz
```

minimap2 uses `assembly_amplicons.fasta` as its reference and maps accurate
reads back to the extracted products. MLVAMaps counts primary alignments and
aligned reference bases for each product.

This stage returns:

- `read_support.sam`
- `read_support.tsv`

## 4. Add optional SAM/BAM support

```bash
mlvamaps call primers.tsv assembly.fasta --bam assembly_reads.bam
```

For an existing assembly-aligned SAM or BAM, MLVAMaps measures alignment-block
overlap with each extracted product. Secondary, supplementary, and unmapped
records are ignored.

This stage returns `read_support.tsv`.

FASTQ and SAM/BAM support are mutually exclusive for one call. Read support
adds mapped-read and mean-coverage evidence; it does not replace the
assembly-product sequence or alter its size-derived repeat count.

## 5. Fingerprint and report

Assembly calls use the same fingerprint, profile matching, and novelty formats
as FASTQ calls.

This stage returns:

- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `novelty_scores.tsv`
- `report.html`

The report displays product calls, coordinates, profile matches, and a generated
gel. Band position follows product size. Band intensity follows read support
when available and otherwise uses a uniform default intensity.
