# Reference building

mlvamaps separates reference acquisition from reference construction:

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
  -p panels/b_cereus_group.csv \
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
  -p panels/mlva_loci.tsv \
  --output references \
  --threads 32
```

This creates independent databases:

```text
references/
├── taxon_reference_summary.tsv
├── taxon_locus_amplifiability.tsv
├── bacillus_cereus_group/
│   ├── prepared/
│   │   ├── package/
│   │   ├── metadata.tsv
│   │   ├── ncbi_dataset.zip
│   │   └── download_manifest.json
│   └── reference/
│       ├── database/
│       ├── phylogeny/
│       ├── reference_locus_amplifiability.tsv
│       └── reference_build_manifest.tsv
├── staphylococcus_aureus/
│   ├── prepared/
│   └── reference/
└── reference_pipeline_manifest.json
```

All taxids in one invocation use the same `-p` panel. Run
separate commands when taxa require different MLVA schemes.

The two top-level TSVs summarize panel compatibility across taxa.
`taxon_reference_summary.tsv` has one row per taxon, while
`taxon_locus_amplifiability.tsv` has one row per taxon/locus and can be pivoted
into a heatmap using `percent_genomes_amplifiable`. A locus is amplifiable when
at least one examined genome has a product retained by the normal in-silico PCR
and multiple-product policy. Thus `NO_AMPLICONS` is distinct from
`INSUFFICIENT_REFERENCES`: the latter is amplifiable but lacks enough retained
references for `--min-references-per-tree`.

Taxon status is `BUILT` only when every panel locus has a tree, `PARTIAL` when
at least one locus is amplifiable but the full panel was not built, and
`NO_USABLE_LOCI` when no locus is amplifiable. A no-usable-loci taxon skips
phylogeny construction without stopping later taxon rows. Unless `--quiet` is
used, the command also prints this summary after each taxon.

`-p` accepts either the rich comma- or tab-delimited locus table used by
`mlvamaps call` or a minimal three-column primer CSV/TSV. Header aliases are
normalized in both cases, including
`name`/`locus`/`id`, `forward`/`fwd`, and `reverse`/`rev`.

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
  -i references/bacillus_cereus_group/prepared/package \
  --metadata references/bacillus_cereus_group/prepared/metadata.tsv \
  -p panels/b_cereus_group.tsv \
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
