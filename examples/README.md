# Example input files

`make_illumina_example.py` creates a tiny offline paired-end dataset, rich
panel, truth assembly, and SRA-style metadata. See the
[Illumina workflow](../docs/workflows/illumina.md#synthetic-worked-example) for
the complete call, validation, and MYOGA commands.

`mlva_loci.example.tsv` is the richer optional format for mlvamaps. It includes
primer sequences plus repeat motif, flank, coordinate, and amplicon metadata.

Optional columns:

- `chrom_or_contig`, `start`, `end`: reference coordinates when known.
- `left_flank_sequence`, `right_flank_sequence`: sequence immediately outside
  the repeat region, used to improve repeat extraction.
- `pool_id`: primer-pool or multiplex identifier.

`amplirust_primers.example.csv` shows the comma-delimited CSV shape that
mlvamaps can write for `amplirust`. Users normally provide
`mlvamaps_primers.example.tsv`, the raw legacy primer file, or `mlva_loci.tsv`.
