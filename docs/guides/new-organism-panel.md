# Adapting mlvamaps to a new organism or MLVA scheme

mlvamaps is organism-agnostic. Supporting another bacterium, fungus, parasite,
or other microbe is primarily a panel-definition and validation task rather
than a code change.

## 1. Start from an established scheme

Record the source of the scheme, primer sequences, primer orientation, repeat
unit definition, product-size convention, and any laboratory-specific allele
numbering rules. Two schemes that reuse a locus name can still define its
repeat count differently.

## 2. Build the minimal panel

Create one row per locus with:

```text
locus_id
forward_primer
reverse_primer
```

Use stable, unique locus identifiers. Keep primer sequences in their published
5-prime to 3-prime orientation and retain IUPAC degenerate bases.

## 3. Add interpretation metadata

Recommended fields:

- `repeat_unit_length_bp`
- `expected_product_size_bp`
- `nominal_repeat_units`
- `repeat_motif`
- `expected_min_repeats`
- `expected_max_repeats`
- `expected_amplicon_min_bp`
- `expected_amplicon_max_bp`

Optional left/right repeat flanks improve repeat-boundary localization in
FASTQ reads. Optional reference coordinates document the scheme but do not
convert representative-mapping positions into genomic coordinates.

## 4. Validate with assemblies

Run representative assemblies whose expected MLVA results are already known:

```bash
mlvamaps call -p new_panel.tsv -i known_assembly.fasta -o validation/assembly
```

Check:

- Expected products are present once at the expected loci.
- Product orientations and coordinates are sensible.
- Product sizes fall within panel bounds.
- Size-to-repeat conversion reproduces the scheme's published convention.
- Primer mismatch allowances do not create off-target products.

## 5. Validate with amplicon reads

```bash
mlvamaps call -p new_panel.tsv -i known_amplicons.fastq.gz -o validation/reads
```

Inspect assignment, repeat features, mapped product groups, SPOARS assemblies,
mapping coverage, and SNP evidence. Include negative controls,
mixed alleles when relevant, low-depth samples, and reads with realistic
sequencing errors.

## 6. Build a profile database

Use a TSV containing `profile_id`, optional `strain_id`, and one column per
locus. Confirm that every profile uses the same repeat-count convention as the
panel.

Do not silently combine databases from different laboratories when their
schemes, primer locations, repeat units, or allele rounding differ. Document
conversions and keep original source identifiers.

## 7. Establish acceptance criteria

Before operational use, define:

- Required loci and allowed dropout.
- Minimum read depth.
- Posterior and multi-variant handling.
- Accepted primer mismatches.
- Competitive mapping specificity and locus-score separation.
- Representative-mapping and SNP evidence thresholds.
- How ambiguous, out-of-range, mixed, and novel-looking results are reviewed.

## 8. Preserve validation evidence

Keep the panel version, source publications, known sample results, command
lines, mlvamaps version, and complete output directories. A profile match is
only meaningful when the typing scheme and panel version are traceable.

## Current scope

mlvamaps can use any microbial MLVA scheme expressible as paired primers plus
repeat/product metadata. It does not yet discover VNTR panels automatically,
translate between incompatible allele-numbering schemes, or validate a new
scheme without known truth data.
