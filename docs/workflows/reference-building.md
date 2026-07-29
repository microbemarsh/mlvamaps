# Reference building

MLVAMaps separates reference acquisition from reference construction:

1. `prepare-reference` downloads a portable NCBI Datasets package, extracts its
   assemblies, normalizes metadata, and records a checksum and tool versions.
2. `build-reference` extracts each MLVA locus and builds its fixed reference
   alignment and phylogeny.

`build-reference` can run both stages as one pipeline, or it can retain its
original local-assembly behavior.

## One taxid

```bash
mlvamaps build-reference \
  --taxid 86661 \
  --loci panels/b_cereus_group.csv \
  --assembly-source refseq \
  --output references \
  --threads 32
```

The resulting database is `references/taxid_86661/reference/database`.

## A CSV of taxids

The input must have a `taxid`, `taxon_id`, or `ncbi_taxid` column. An optional
`name`, `reference_name`, or `database_name` column controls the output
directory:

```csv
taxid,name
86661,bacillus_cereus_group
1280,staphylococcus_aureus
```

```bash
mlvamaps build-reference \
  --taxids-csv taxa.csv \
  --loci panels/mlva_loci.tsv \
  --output references \
  --threads 32
```

This creates independent databases:

```text
references/
├── bacillus_cereus_group/
│   ├── prepared/
│   │   ├── package/
│   │   ├── metadata.tsv
│   │   ├── ncbi_dataset.zip
│   │   └── download_manifest.json
│   └── reference/
│       ├── database/
│       ├── phylogeny/
│       └── reference_build_manifest.tsv
├── staphylococcus_aureus/
│   ├── prepared/
│   └── reference/
└── reference_pipeline_manifest.json
```

All taxids in one invocation use the same `--primers` or `--loci` panel. Run
separate commands when taxa require different MLVA schemes.

`--loci` accepts the same rich comma- or tab-delimited locus table used by
`mlvamaps call`. Minimal three-column primer CSV/TSV files can instead be
passed with `--primers`.

## Prepare now, build later

To acquire inputs without performing the computational build:

```bash
mlvamaps prepare-reference \
  --taxids-csv taxa.csv \
  --assembly-source refseq \
  --output references
```

The prepared metadata includes an `assembly_file` path relative to the
extracted package, so each cohort can later be built through the local input
interface:

```bash
mlvamaps build-reference \
  --assemblies references/bacillus_cereus_group/prepared/package \
  --metadata references/bacillus_cereus_group/prepared/metadata.tsv \
  --loci panels/b_cereus_group.tsv \
  --output references/bacillus_cereus_group/reference
```

Use `--resume` only to reuse an existing
`prepared/ncbi_dataset.zip` after an interrupted acquisition. Published
reference outputs should otherwise be treated as immutable.

## Requirements and selection

Taxid acquisition requires NCBI's `datasets` and `dataformat` executables. They
are included in the conda environment through `ncbi-datasets-cli`.

`--assembly-source refseq` is the default. `genbank` selects GenBank assemblies,
and `all` asks NCBI for both sources. Extra NCBI Datasets options can be passed
by repeating `--datasets-arg`; use the `--datasets-arg=--option` form for values
that begin with a dash.

NCBI downloads are retried three times by default when a transient network
error or incomplete ZIP is detected. Change this with `--download-retries N`.
An incomplete archive is removed before retrying, while a completed,
ZIP-validated archive can be reused with `--resume`.
