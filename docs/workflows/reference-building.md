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

## Internal stages

The command fetches references, extracts real primer-bounded loci, summarizes
amplifiability, builds real-reference MAFFT/RAxML-NG assets, generates bounded
candidate repeat states, deduplicates identical sequence hypotheses, builds
short- and long-read minimap2 indexes, builds a Deacon index from complete real
genomes, and finally writes `manifest.json`.

`--threads` sets the overall build budget and remains available to parallel
extraction, MAFFT, and other threaded stages. RAxML-NG reference-tree searches
are safely run with one internal thread per process; users do not need to reduce
the entire build to `--threads 1`.

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

Use the completed build directly:

```bash
mlvamaps call -i sample.fasta --database mlvamaps_db -o results/sample
```

Databases without schema 2.0 competitive assets must be rebuilt. NCBI package
downloads can be reused with `--resume` after an interrupted build.