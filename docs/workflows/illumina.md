# Illumina short-read workflow

The default `--sr-engine repeat-likelihood` recruits molecules by non-repeat
flanks, then combines E (enclosing), S (spanning pair), F (flanking), and anchored
FRR evidence. Repeat states are synthetic and can be novel or outside the panel's
expected range. Haploid likelihood and the existing microbial EM model retain
distinct components, including supported SNP haplotypes with the same count.

```bash
mlvamaps call --loci panel.tsv -i sr --fq1 reads_1.fastq.gz --fq2 reads_2.fastq.gz \
  --insert-mean 400 --insert-sd 35 --sample-mode metagenome -o results
```

Insert overrides must describe independently measured library statistics. They
are optional when repeat-independent flank pairs can estimate them or enclosing
reads determine the repeat directly. Single-end data can use E and F evidence;
FRRs without a locus anchor are excluded. Unidentifiable lengths remain blank.

Read-supported non-repeat sequence is reconstructed after repeat inference;
uncovered bases become N. Products enter the same genotype and SNP-mask path as
assembly and LR products. See [formulas, diagnostics, CLI semantics and
limitations](../concepts/locus-reconstruction.md). For regression comparison,
`--sr-engine competitive` retains the [legacy workflow](competitive-fastq.md).
