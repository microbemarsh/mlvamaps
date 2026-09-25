# Workflow architecture

mlvamaps has one database-construction workflow and two calling entry paths:
assembly FASTA and sequencing FASTQ. Both calling paths produce compatible
locus calls, fingerprints, profile matches, and reports, but they preserve the
different evidence available from assembled contigs and individual molecules.

Short reads default to [competitive candidate mapping](competitive-fastq.md).
Long reads default to spanning-molecule recovery. The optional SR
repeat-likelihood engine and default LR engine are documented in
[locus reconstruction](../concepts/locus-reconstruction.md).

## Resource ownership

Reusable resources belong to the completed database, not to an individual
taxon or sample:

- `database/competitive_mapping/` owns the candidate FASTA, candidate metadata,
  provenance, and short-/long-read minimap2 indexes.
- `database/deacon/` owns the broad target-group recruitment reference and its
  Deacon index.
- Per-locus FASTAs contain observed reference amplicons for classification.
  Classification aligns Sassy-derived assembly products or FASTQ molecules to
  these sequences and combines sequence and repeat-length evidence.

A taxid build first retains each taxon's real amplicons and extraction QC in an
isolated work directory. It then merges all successful taxon outputs and builds
each reusable resource once. This is also the path for a one-row taxid build;
there is no separate per-taxon finalization followed by a duplicate combined
finalization.

Local-assembly `build-reference` input already represents one complete cohort,
so it proceeds directly from extraction to one finalization.

## Locus reconstruction and candidate mapping

Assembly products, complete LR molecules, and SR likelihood/sequence
reconstructions enter `locus_products.genotype_product`. The optional SR model uses
panel flanks or borrowed database context to identify loci and synthesizes
repeat states independently of known alleles. Its candidate scores and per-pair
evidence are retained for inspection.

`--sr-engine competitive` (the default) and `--lr-engine competitive` use the
`unified_fastq` candidate-mapping workflow and cached minimap2 indexes.
See [competitive workflows](competitive-fastq.md) for those details.

## Global thread budget

`--threads` is the process-wide CPU budget:

- reference extraction uses at most the resolved budget;
- native tools receive the threads allocated to their active stage;
- Illumina directory and manifest batches run independent samples concurrently;
- active sample count is bounded by sample count, CPU budget, and a memory-aware
  concurrency cap;
- each active sample receives `floor(total_threads / active_samples)` threads,
  so allocated native threads do not exceed the global budget.

The default batch cap is four active samples. Set it lower for large datasets or
memory-constrained nodes:

```bash
export MLVAMAPS_MAX_CONCURRENT_SAMPLES=2
mlvamaps call -p panel.tsv -i reads/ --short-reads -o results -t 32
```

Input discovery and final combined tables retain deterministic input order even
when samples finish out of order. Each sample has an isolated output directory,
and one failed sample is recorded without being interpreted as a biological
negative result.

## Progress and auditability

Unless `--quiet` is used, elapsed-time messages report phase boundaries and
counts such as taxa, assemblies, loci, candidates, and workers where available.
These messages are operational diagnostics, not biological QC. Biological
retention, dropout, ambiguity, and exclusion reasons remain in TSV and JSON
outputs and must be reviewed separately.

## Stable result boundary

Implementation details may differ between input technologies, but these public
outputs remain available:

- `calls.tsv`
- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `profile_match_loci.tsv`
- `report.html`

Do not infer biological absence from a blank repeat count. Consult the locus
status and evidence columns to distinguish not found, insufficient depth,
presence without sizing, ambiguity, filtering, and mixture evidence.
