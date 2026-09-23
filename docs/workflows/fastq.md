# Accurate long-read and amplicon workflow

The default `--lr-engine spanning` filters reads, identifies both primer anchors
using Sassy, orients and trims complete molecules, retains distinct sequence
clusters, and applies the existing sequence mixture EM. Partial molecules retain
presence evidence. Novel repeat counts do not require database representation.

```bash
mlvamaps call --loci panel.tsv -i accurate_reads.fastq.gz -o results --lr-engine spanning
```

Products enter the same calibrated repeat and masked SNP typing functions used
for assembly and reconstructed SR products. See [reconstruction details and
outputs](../concepts/locus-reconstruction.md). Accurate molecules are required;
the default does not perform noisy-read correction. `--lr-engine competitive`
retains the [legacy minimap2/SPOARS workflow](competitive-fastq.md).
