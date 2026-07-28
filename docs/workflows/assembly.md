# Assembly workflow

The assembly path accepts FASTA assemblies and can optionally add read support.

```bash
mlvamaps call primers.tsv assembly.fasta -o results
```

## 1. Find paired-primer products

Amplirust searches both assembly orientations with degenerate-primer and
primer-error support. For compatibility, assembly gaps represented by `N`
bases are permitted inside products; MLVA_finder did not reject them. MLVAMaps
filters circular-wrap records and uses the historical raw-allele limit for
discovery. Novel calling still enforces each locus's configured amplicon range
before selecting a product.

This stage returns:

- `assembly_amplicons.tsv`: locus, contig, 1-based coordinates, orientation,
  product size, and primer mismatches.
- `assembly_amplicons.fasta`: extracted products.
- Native evidence under `amplirust/`.

The default `--algorithm legacy` path reproduces MLVA_finder's decision rule:
it uses the lowest successful primer-error round, applies the original
forward-strand and equal-length-match precedence, and retains the smallest raw
allele on the last FASTA record with a valid product in that round. Product
sizes also follow the original configured-primer-length formula when a primer
match contains an indel. `--algorithm novel` selects the probabilistic caller
described below. All accepted products remain available in the amplicon table.
The regression suite includes a self-contained MLVA_finder oracle covering
perfect and mismatched primers, strict half-unit rounding, multiple FASTA
records, and parallel PCR execution.

## 2. Estimate repeat count

Both algorithms convert product size to repeat count when the panel provides
repeat-unit length, nominal repeat units, and expected product size. It retains
the raw estimate. The default legacy algorithm applies MLVA_finder's strict
integer tolerance (configured with `--assembly-round-tolerance`) and otherwise
uses the intervening half allele.

With `--algorithm novel`, MLVAMaps evaluates the raw estimate against explicit
integer and half-unit allele states. The reported probability reflects distance
from those states at the locus's repeat-unit resolution. Exact midpoint ties
remain `AMBIGUOUS` instead of being resolved by an arbitrary rounding rule.
`--min-posterior` controls the novel caller's required confidence.

In novel mode, when FASTQ or BAM support is supplied and the assembly contains
multiple accepted products for a locus, mapped-read counts weight the
product-length distribution. This lets the supported allele win rather than
defaulting to the smallest legacy allele. Legacy mode records depth for its
selected product but does not allow depth to change the historical call.

Assembly statuses:

- `PASS`: product found and repeat count could be calculated.
- `AMBIGUOUS`: the leading allele is below `--min-posterior` or leads the
  runner-up by less than 0.2.
- `PRESENT_COUNT_UNKNOWN`: product found, but panel metadata was insufficient
  for a repeat-count calculation.
- `NOT_FOUND`: no accepted paired-primer product was found.

This stage returns `calls.tsv`.

For direct migration checks, every assembly run also writes CSV compatibility
views: `legacy_output.csv`, `legacy_mlva_analysis.csv`,
`legacy_predicted_pcr_sizes.csv`, and `legacy_primer_mismatches.txt`. These sit
alongside the richer TSV, FASTA, evidence, and HTML outputs; they do not replace
them.

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

Assembly calls use the same fingerprint and profile-matching formats as FASTQ
calls.

This stage returns:

- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `report.html`

The report leads with panel recovery, locus-specific review findings, closest
profile/reference interpretation, and the generated gel. Product coordinates,
primer mismatches, and alternative-call probabilities remain in expandable
detail tables. Band position follows product size; band intensity follows read
support when available and otherwise uses a uniform default intensity.
