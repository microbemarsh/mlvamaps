# Example input files

`mlvamaps_primers.example.tsv` is the recommended minimal primer format for
MLVAMaps. It is a cleaned, headered version of the legacy MLVA_finder primer
file and contains one row per VNTR locus.

Required primer columns:

- `locus_id`: stable locus name used in all outputs.
- `forward_primer`: forward primer sequence, 5' to 3'.
- `reverse_primer`: reverse primer oligo sequence, 5' to 3'. MLVAMaps searches
  for its reverse complement in oriented reads.

Important VNTR columns:

- `repeat_unit_length_bp`: length of the repeat unit, parsed from legacy names
  like `vrrA_12bp_314bp_10U`.
- `expected_product_size_bp`: nominal PCR product size including primers.
- `nominal_repeat_units`: nominal repeat count from the legacy primer name.

`insilicoMLVAprimers_all.raw.example.csv` is a direct copy of the legacy-style
file. Despite the `.csv` suffix, it is whitespace-delimited:

```text
locus_id forward_primer reverse_primer
```

MLVAMaps can ingest either file.

For assembly calls:

```bash
mlvamaps call examples/seer_lab_Ba/mlvamaps_primers.example.tsv assembly.fasta
```

For FASTQ reads:

```bash
mlvamaps call examples/seer_lab_Ba/mlvamaps_primers.example.tsv reads.fastq.gz
```

`mlva_loci.example.tsv` is the richer optional format for MLVAMaps. It includes
primer sequences plus repeat motif, flank, coordinate, and amplicon metadata.

Optional columns:

- `chrom_or_contig`, `start`, `end`: reference coordinates when known.
- `left_flank_sequence`, `right_flank_sequence`: sequence immediately outside
  the repeat region, used to improve repeat extraction.
- `pool_id`: primer-pool or multiplex identifier.

`amplirust_primers.example.csv` shows the comma-delimited CSV shape that
MLVAMaps can write for `amplirust`. Users normally provide
`mlvamaps_primers.example.tsv`, the raw legacy primer file, or `mlva_loci.tsv`.
