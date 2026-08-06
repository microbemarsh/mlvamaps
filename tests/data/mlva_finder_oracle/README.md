# MLVA_finder compatibility oracle

`expected.csv` was generated from `assembly.fasta` and `primers.tsv` with
i2bc/MLVA_finder commit `15dacfc413b7fcb8c432c1cc8700e27e643ce0c3`:

```bash
python MLVA_finder.py \
  --input INPUT_DIRECTORY \
  --output OUTPUT_DIRECTORY \
  --primer primers.tsv \
  --contig
```

The fixture exercises exact products, a one-error primer product, strict
half-unit rounding, multiple FASTA records, and MLVA_finder's last-matching
record selection rule. mlvamaps runs the same fixture with multiple PCR workers
to ensure parallel discovery does not alter result order.
