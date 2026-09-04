# Assembly workflow

The assembly path accepts FASTA assemblies and can optionally add read support.
It uses the same Sassy-backed batch in-silico PCR implementation used during
reference extraction and for FASTQ-derived local consensus products; there is
no separate assembly-only primer-matching algorithm.

```bash
mlvamaps call -p primers.tsv -i assembly.fasta -o results
```

## 1. Find paired-primer products

`mlvamaps` uses the `sassy search` command-line tool to discover approximate
primer matches. It then applies its own MLVA_finder-compatible IUPAC expansion,
strand logic, primer pairing, length filtering, and product selection:

- Both primers are supplied 5-prime to 3-prime. In each oriented sequence,
  `mlvamaps` searches for the forward primer and the reverse complement of the
  reverse primer, then pairs downstream reverse-primer matches with each
  forward-primer match.
- Configured IUPAC primer codes are expanded deterministically before Sassy is
  called. Assembly `N` bases are not treated as primer wildcards; overlap with
  an ambiguous target base consumes edit distance. Assembly gaps represented by
  `N` are nevertheless permitted elsewhere inside a recovered product.
- Sassy's own reverse-complement search is disabled because `mlvamaps` controls
  orientation. For each concrete forward-primer expansion it searches the input
  contig first and searches the reverse-complemented contig only as a fallback
  when no input-strand forward match is found. The reported orientation is
  relative to the original contig.
- `--max-primer-mismatches` defaults to 2 and is applied independently to each
  primer as edit distance (substitutions plus insertions and deletions).
  Searches run cumulatively from error round 0 through the configured maximum.
- A candidate must have the primers in order and satisfy the PCR engine's
  inclusive product-length bounds. Assembly calls use global bounds broad enough
  to cover every locus's historical MLVA_finder-valid raw allele (`0 <= raw
  repeat count < 100`); this is deliberately not the same as enforcing every
  rich panel's `expected_amplicon_min_bp`/`expected_amplicon_max_bp` per locus.
- Sassy discovers candidate locations. The `regex` package is used only in
  windows around those locations to reproduce legacy equal-cost indel traceback
  choices; it is not the genome-wide primer search engine.

This stage returns:

- `assembly_amplicons.tsv`: locus, contig, 1-based coordinates, orientation,
  product size, and primer mismatches.
- `assembly_amplicons.fasta.gz`: extracted products.
- Sassy-backed evidence under `in_silico_pcr/` (`primers.csv`, `matches.tsv`,
  and `products.fasta.gz`).

The assembly caller reproduces MLVA_finder's decision rule. It uses the lowest
successful primer-error round, applies input-strand and equal-length-match
precedence, and retains the smallest raw repeat count on the last FASTA record
with a valid product in that round. Product sizes also follow the original
configured-primer-length formula when a primer match contains an indel. All
recovered candidate products remain available in `assembly_amplicons.tsv`, even
though only the selected product supplies the locus repeat-count call.
The regression suite includes a self-contained MLVA_finder oracle covering
perfect and mismatched primers, strict half-unit rounding, multiple FASTA
records, and parallel PCR execution.

When a directory of assemblies is supplied, samples have isolated output
directories and are combined in deterministic filename order. See
[workflow architecture](architecture.md) for the global thread-budget and
sample-identity rules.

## 2. Estimate repeat count

The caller converts product size to repeat count when the panel provides
repeat-unit length, nominal repeat units, and expected product size. It retains
the raw estimate and applies MLVA_finder's strict
integer tolerance (configured with `--assembly-round-tolerance`) and otherwise
uses the intervening half allele.

FASTQ or BAM support records depth for the selected product but does not alter
the historical product selection or repeat count.

Assembly statuses:

- `PASS`: product found and repeat count could be calculated.
- `NOT_FOUND`: no recovered product was eligible for the historical
  repeat-count selection rule. A primer product can therefore remain visible in
  `assembly_amplicons.tsv` while the locus has no repeat-count call when the
  panel cannot convert product length to repeat count or the raw count is 100 or
  greater.

This stage returns `calls.tsv`.

For direct migration checks, every assembly run also writes CSV compatibility
views: `legacy_output.csv`, `legacy_mlva_analysis.csv`,
`legacy_predicted_pcr_sizes.csv`, and `legacy_primer_mismatches.txt`. These sit
alongside the richer TSV, FASTA, evidence, and HTML outputs; they do not replace
them.

## 3. Add optional FASTQ support

```bash
mlvamaps call -p primers.tsv -i assembly.fasta --reads sample.fastq.gz
```

minimap2 uses `assembly_amplicons.fasta.gz` as its reference and maps accurate
reads back to the extracted products. mlvamaps counts primary alignments and
aligned reference bases for each product.

This stage returns:

- `read_support.sam`
- `read_support.tsv`

## 4. Add optional SAM/BAM support

```bash
mlvamaps call -p primers.tsv -i assembly.fasta --bam assembly_reads.bam
```

For an existing assembly-aligned SAM or BAM, mlvamaps measures alignment-block
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
