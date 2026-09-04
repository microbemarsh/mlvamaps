# Reference building

`build-reference` is the complete database-construction workflow for mlvamaps.
It accepts one NCBI taxonomy ID, a `taxids.csv` cohort, or local assemblies plus
metadata. No separate preparation or locus-extraction command is required.

```bash
mlvamaps build-reference \
  --taxids taxids.csv \
  --panel primer.csv \
  --output mlvamaps_db \
  --threads 32
```

`taxids.csv` requires `taxid` (aliases `taxon_id` and `ncbi_taxid`) and may
include `name`:

```csv
taxid,name
86661,bacillus_cereus_group
1280,staphylococcus_aureus
```

For local material:

```bash
mlvamaps build-reference -i assemblies/ --metadata metadata.tsv \
  --panel primer.csv --output mlvamaps_db --threads 32
```

## Extract, merge, finalize

Reference construction deliberately separates taxon-local extraction from
database-wide finalization.

1. **Prepare each taxon.** NCBI packages and normalized metadata are retained in
   the taxon's work directory.
2. **Extract each taxon.** The shared Sassy-backed PCR engine examines assemblies
   and retains real primer-bounded amplicons, product-selection evidence, and
   amplifiability QC. It does not build taxon-local candidate indexes, a Deacon
   index, or taxon-local phylogenies.
3. **Merge real observations.** Amplicons, metadata, and source-assembly records
   are merged by locus. A reference identifier present in multiple taxon cohorts
   is an error rather than a silent many-to-one merge. Metadata for references
   without any usable locus is excluded from the callable database while the
   taxon's failure remains visible in summary QC.
4. **Finalize once.** The merged database generates one deduplicated candidate
   bank, one short-read index, one long-read index, one broad Deacon index, and
   one set of real-reference phylogenies.

A single `--taxid` follows exactly the same coordinator path and is finalized
once. Local `-i ASSEMBLIES --metadata ...` input is already one cohort and is
also finalized once.

Unless `--quiet` is selected, progress output identifies the current taxon,
assembly extraction counts, phase, locus counts, worker counts, and elapsed
phase time. This is intended to make a cluster run diagnosable without a
profiler.

`--threads` is the overall build CPU budget. Extraction workers do not exceed
that budget. RAxML-NG is always invoked with one internal thread; users should
not reduce the whole build to one thread merely to protect short-locus RAxML-NG
runs.

Candidate generation uses `expected_min_repeats` and `expected_max_repeats`
when explicitly supplied. Otherwise it conservatively uses observed calibrated
alleles plus one state on each side, capped to a 20-state window and an absolute
maximum of 100. Synthetic editing changes only the defined repeat interval.

## Database layout

```text
mlvamaps_db/
├── manifest.json
├── database/
│   ├── reference_panel.tsv
│   ├── reference_metadata.tsv
│   ├── reference_assemblies.tsv
│   ├── LOCUS.fasta.gz
│   ├── competitive_mapping/
│   │   ├── candidate_contexts.fasta
│   │   ├── candidate_metadata.tsv
│   │   ├── candidate_provenance.tsv
│   │   ├── short.mmi
│   │   └── long.mmi
│   └── deacon/
│       ├── reference_genomes.fasta
│       └── target_recruitment.idx
├── phylogeny/
├── reference_build_manifest.tsv
└── reference_locus_amplifiability.tsv
```

Multi-taxon builds retain isolated taxon work directories and write
`taxon_reference_summary.tsv` plus `taxon_locus_amplifiability.tsv`. The
combined Deacon index includes all requested real genomes, so recruitment does
not preselect one taxon. Real locus FASTAs and phylogenies never contain
synthetic candidate alleles. Candidate prevalence is metadata rather than
duplicated mapping evidence.

The isolated taxon directories are extraction and QC artifacts, not standalone
fully indexed databases. The supported call target is the completed top-level
database.

## Candidate and recruitment strategy

Candidate generation starts from distinct observed locus backgrounds before
expanding bounded whole-repeat states. Identical locus/sequence/state hypotheses
are collapsed to one candidate while `candidate_provenance.tsv` retains all
source reference, taxon, sequence, and background links. Candidate identifiers
are stable within the completed build, and calls reuse these files rather than
regenerating them.

The Deacon reference contains the requested cohort's real source genomes and
excludes synthetic candidate alleles. This broad combined strategy avoids
preselecting one taxon before metagenomic recruitment. The source assembly table
and build manifest provide the retained source identifiers and checksums.

Use the completed build directly:

```bash
mlvamaps call -i sample.fasta --database mlvamaps_db -o results/sample
```

Databases without schema 2.0 competitive assets must be rebuilt. NCBI package
downloads can be reused with `--resume` after an interrupted build.