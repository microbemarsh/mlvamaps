# Example input files

`mlva_seer_primers.example.tsv` is the recommended minimal primer format for
MLVA Seer. It is a cleaned, headered version of the legacy MLVA_finder primer
file and contains one row per VNTR locus.

Required primer columns:

- `locus_id`: stable locus name used in all outputs.
- `forward_primer`: forward primer sequence, 5' to 3'.
- `reverse_primer`: reverse primer oligo sequence, 5' to 3'. MLVA Seer searches
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

MLVA Seer can ingest either file with `--primers`.

For assembly/contig/reference extraction, only `--input` and `--primers` are
required:

```bash
mlva-nanopore extract-amplicons \
  --input assembly.fasta \
  --primers examples/mlva_seer_primers.example.tsv
```

For FASTQ reads:

```bash
mlva-nanopore call \
  --input reads.fastq.gz \
  --primers examples/mlva_seer_primers.example.tsv
```

`mlva_loci.example.tsv` is the richer optional format for MLVA Seer. It includes
primer sequences plus repeat motif, flank, coordinate, and amplicon metadata.

Optional columns:

- `chrom_or_contig`, `start`, `end`: reference coordinates when known.
- `left_flank_sequence`, `right_flank_sequence`: sequence immediately outside
  the repeat region, used to improve repeat extraction.
- `pool_id`: primer-pool or multiplex identifier.

`amplirust_primers.example.csv` shows the comma-delimited CSV shape that MLVA
Seer writes internally for `amplirust`. Users normally provide
`mlva_seer_primers.example.tsv`, the raw legacy primer file, or `mlva_loci.tsv`;
the `extract-amplicons` command converts it to the amplirust primer CSV.
