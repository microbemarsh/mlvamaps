# Workflow architecture

mlvamaps has one database-construction workflow and two calling entry paths:
assembly FASTA and sequencing FASTQ. Both calling paths produce compatible
locus calls, fingerprints, profile matches, and reports, but they preserve the
different evidence available from assembled contigs and individual molecules.

```text
Reference construction
  taxon genomes ──> per-taxon amplicons + QC ──┐
                                                ├─> merged real references
  local genomes ────────────────────────────────┘          │
                                                          ├─> candidate bank + minimap2 indexes (once)
                                                          ├─> Deacon index (once)
                                                          └─> locus phylogenies (once)

Calling
  assembly FASTA ──> shared in-silico PCR ──> assembly locus calls ──┐
                                                                    ├─> fingerprints, taxonomy, report
  FASTQ ──> QC ──> cached candidates/index ──> molecule evidence ───┘
```

## Resource ownership

Reusable resources belong to the completed database, not to an individual
taxon or sample:

- `database/competitive_mapping/` owns the candidate FASTA, candidate metadata,
  provenance, and short-/long-read minimap2 indexes.
- `database/deacon/` owns the broad target-group recruitment reference and its
  Deacon index.
- `phylogeny/` owns alignments and trees built only from observed reference
  amplicons. Synthetic repeat-state candidates never enter these trees.

A taxid build first retains each taxon's real amplicons and extraction QC in an
isolated work directory. It then merges all successful taxon outputs and builds
each reusable resource once. This is also the path for a one-row taxid build;
there is no separate per-taxon finalization followed by a duplicate combined
finalization.

Local-assembly `build-reference` input already represents one complete cohort,
so it proceeds directly from extraction to one finalization.

## Shared candidate mapping

Illumina, ONT, ONT-HQ, and HiFi calls use the same versioned candidate contexts.
With `--database`, candidates and the appropriate minimap2 index are loaded from
the completed database. They are not regenerated in each sample directory.
Without a database, a rich panel may synthesize a bounded local candidate bank;
a primer-only panel does not contain enough information to do so.

minimap2 emits SAM on a pipe. htslib, through `pysam`, decodes it and writes a
compressed temporary BAM while the candidate alignments are interpreted. Normal
operation removes that BAM. `--keep-intermediates` retains
`candidate_mapping/candidate_alignments.bam`; this option does not produce a
text `candidate_alignments.sam`.

The current Python-facing evidence API preserves competing taxa and repeat
states after htslib decoding. Sequence and quality evidence retain molecule
identity, including Illumina mate identity. Missing or presence-only evidence is
not converted into an allele value.

## Global thread budget

`--threads` is the process-wide CPU budget:

- reference extraction uses at most the resolved budget;
- native tools receive the threads allocated to their active stage;
- every RAxML-NG process is forced to `--threads 1`;
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
outputs retain their established schemas:

- `calls.tsv`
- `mlva_fingerprint.tsv`
- `mlva_fingerprint_probabilistic.tsv`
- `profile_matches.tsv`
- `profile_match_loci.tsv`
- `report.html`

Do not infer biological absence from a blank repeat count. Consult the locus
status and evidence columns to distinguish not found, insufficient depth,
presence without sizing, ambiguity, filtering, and mixture evidence.