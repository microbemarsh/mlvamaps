# mlvamaps

`mlvamaps` calls microbial MLVA/VNTR loci from Illumina reads, accurate long or
amplicon reads, and genome assemblies. It uses a user-supplied primer panel, so
no organism or typing scheme is hard-coded.

The main outputs are an MLVA fingerprint, per-locus calls and evidence, and a
self-contained HTML report. Reference databases add alignment-based reference support using nucleotide
variation and VNTR repeat-length evidence.

```mermaid
flowchart TD
    CLI[mlvamaps call] --> Route{Input selector}
    Route -->|sr or short-reads directory| SR[Short-read recovery]
    Route -->|FASTQ file| LR[Long-read recovery]
    Route -->|FASTA file| AS[Assembly PCR and product selection]
    DB[Database panel and observed locus sequences] --> SR
    DB --> LR
    SR --> Shared[Shared molecule measurement and allele inference]
    LR --> Shared
    Shared --> Calls[Calls and MLVA fingerprint]
    AS --> Calls
    Shared -->|Reuse observed-reference alignments| Typing[Reference classification]
    AS -->|With database: align selected products| Typing
    Calls --> Report[Profile comparison and HTML report]
    Typing --> Report
```

Assembly and FASTQ inputs retain their own evidence-recovery paths:

- **Genome assemblies:** Sassy-backed, primer-directed in silico PCR recovers
  candidate amplicons; the assembly caller applies MLVA_finder-compatible
  product selection and size-to-repeat conversion.
- **Illumina reads:** minimap2 maps molecules directly to observed database
  loci. Read-based repeat measurements supply the shared allele caller. Paired reads also use flank recruitment with mate
  rescue and, for uncertain loci, bounded local assembly through the shared
  product-length parser.
- **Long-read or amplicon FASTQ:** Sassy primer recovery supplies complete
  products, followed by mapping to observed database loci and shared molecule
  measurement. Primer-recovered products remain usable when mapping misses them. SPOARS
  remains available for representative sequence correction and confirmation;
  it is not required before a repeat-count call.

Locus detection does not by itself imply a confident repeat-count call.
Unresolved and missing calls remain explicit and are never converted to zero.

See [how each workflow works](workflow.md) for source-traced diagrams, processing
steps, and evidence outputs for short reads, long reads, and assemblies.

FASTQ genotyping and reference typing share one competitive mapping to the
actual observed database amplicons. Repeat counts require a measured product,
resolved repeat flanks, or a supported local assembly; a matching reference
alone does not supply the sample's allele. Novel measured alleles remain callable.
The same alignments and molecule-level repeat measurements inform reference
classification, preserving equivalent references and unresolved evidence.

Reference typing uses one Emu-inspired alignment-likelihood framework:
Sassy-recovered assembly amplicons are aligned to observed reference loci with
Parasail; FASTQ molecules are competitively aligned to observed reference
amplicons with minimap2. Nucleotide mismatches and repeat-length disagreement
both contribute. Assemblies use isolate support; mixed read samples use EM.
This is an independent implementation; the Emu package is not required.

Sassy-backed assembly PCR, MLVA_finder-compatible repeat counts and tables,
recovered amplicons, and read SNP evidence remain available. Assembly runs with
a database also retain reference-relative mismatch counts, repeat differences,
and alignment CIGARs in `classification/assembly_reference_evidence.tsv`.

Optional repeat-profile trees visualize sample similarity and do not determine
reference classification. Taxon annotations describe the reference catalog;
they do not turn a marker-similarity tree into a complete taxonomic phylogeny.
The per-sample tree currently uses repeat counts only, even though reference
classification also uses nucleotide variation. See [alignment-based classification](docs/concepts/mapping-classification.md)
and [migration to the emu workflow](docs/workflows/emu-migration.md).

## Install

### Conda/Miniforge (recommended)

[Miniforge](https://github.com/conda-forge/miniforge) provides `conda` on Linux
and macOS. Until the Bioconda package is published, install from this checkout:

```bash
git clone https://github.com/microbemarsh/mlvamaps.git
cd mlvamaps
conda env create -f environment.yml
conda activate mlvamaps
python -m pip install --no-deps .
```

Verify the installation:

```bash
mlvamaps --version
mlvamaps --help
```

After the Bioconda recipe is accepted, installation will be:

```bash
conda create -n mlvamaps -c conda-forge -c bioconda mlvamaps
conda activate mlvamaps
```

> **Bioconda status:** the `mlvamaps` recipe is staged in
> [`packaging/bioconda/meta.yaml`](packaging/bioconda/meta.yaml). Its native
> dependencies include the `sassy>=0.2.2` command-line tool used for primer
> matching. If Sassy is installed outside the active environment, set
> `SASSY_BIN` to its executable path.

## Run your first sample

You need:

1. a FASTA/FASTQ input; and
2. for FASTQ, a completed reference database built for the target taxon and
   primer panel (`--database`); or, for assembly-only calling, a CSV/TSV panel
   with at least `locus_id`, `forward_primer`, and `reverse_primer` columns.

FASTQ commands load the saved database panel automatically. Build the database
once using the [reference-building guide](docs/workflows/reference-building.md).

A richer panel can also describe repeat motifs, flanks, expected repeat ranges,
and accepted amplicon sizes. See the [input format reference](docs/reference/input-formats.md).

### Genome assembly

```bash
mlvamaps call \
  -p examples/mlva_loci.example.tsv \
  -i sample.fasta \
  -o results/sample \
  -t 8
```

### Paired-end Illumina reads

```bash
mlvamaps call \
  --database references \
  -i sr \
  --fq1 sample_R1.fastq.gz \
  --fq2 sample_R2.fastq.gz \
  --sample-id sample \
  -o results/sample \
  -t 8
```

For a directory containing exact `SAMPLE_1.fastq.gz` / `SAMPLE_2.fastq.gz`
pairs:

```bash
mlvamaps call --database references -i reads/ --short-reads -o results -t 8
```

### Accurate long or amplicon reads

```bash
mlvamaps call --database references -i sample.fastq.gz -o results/sample -t 8
```

## Find the results

Start with:

| Output | Purpose |
| --- | --- |
| `report.html` | Human-readable calls, QC, evidence, and matches. |
| `mlva_fingerprint.tsv` | Sample-by-locus repeat-copy-number profile. |
| `calls.tsv` | Tidy per-locus calls and statuses. |
| `locus_repeat_counts.tsv` | Compact individual-locus repeat counts. |

Failed or unresolved loci are reported explicitly rather than silently changed
to zero. See the complete [output reference](docs/reference/outputs.md) for all
evidence and diagnostic files.

Single-input calls write directly to the requested output directory. Directory
and manifest calls write each sample under `OUTDIR/<sample_id>/` and place only
batch-wide aggregate tables and status information under
`OUTDIR/batch_summary/`.

## Build a reference database

Build directly from NCBI assemblies for one taxon:

```bash
mlvamaps build-reference \
  --taxid 86661 \
  -p panel.tsv \
  -o references \
  -t 16
```

The panel may be a minimal `primer.csv` (locus, forward primer, reverse primer)
or a rich CSV/TSV panel. To compare one panel across taxa, provide `taxids.csv`:

```csv
taxid,name
86661,bacillus_cereus_group
1280,staphylococcus_aureus
```

```bash
mlvamaps build-reference \
  --taxids-csv taxids.csv \
  -p panel.tsv \
  -o references \
  -t 16
```

Multi-taxon builds write:

- `taxon_reference_summary.tsv`: one row per taxon;
- `taxon_locus_amplifiability.tsv`: one row per taxon and locus, suitable for
  compatibility heatmaps;
- one isolated extraction/QC work area per taxon; and
- a combined top-level database containing all taxa for automatic taxon
  identification.

Candidate indexes and the Deacon index are built once from the
merged cohort, including when `taxids.csv` contains only one taxon.

A locus is amplifiable when at least one examined genome produces an amplicon
retained by the normal primer-matching and filtering rules. Any retained reference amplicon remains available for comparison, including
loci represented by a single genome. Taxa with no usable loci are recorded
and do not stop later taxa.

Use a built database during calling:

```bash
mlvamaps call \
  -i sample.fasta \
  --database references \
  -o results/sample
```

Current databases retain `competitive_mapping/` assets and a broad real-genome
Deacon recruitment index for database compatibility. FASTQ calling uses the
observed per-locus FASTAs directly, rather than their synthetic allele expansion.
The observed mapping FASTA and metadata are recorded under each sample’s
`candidate_mapping/` directory. All FASTQ commands require `--database`, including
single-end, paired-end, manifest, and directory inputs. A standalone panel or
`--recruitment-database` alone does not satisfy this requirement. Assembly-only
calling can still use a panel without a database.

For a multi-taxid build, that command automatically loads the saved panel and
taxon metadata, then identifies the closest reference IDs with taxon annotations. No separate panel,
target taxid, calibration artifact, or taxon-identification flag is required.
Databases that predate schema 2.0 must be rebuilt with the current
`mlvamaps build-reference` so all competitive-mapping and recruitment assets
have reproducible provenance.

See the [reference-building guide](docs/workflows/reference-building.md) for
local assemblies, metadata, resuming downloads, and output interpretation.
See the [workflow architecture](docs/workflows/architecture.md) for resource
reuse, streamed candidate mapping, and batch thread allocation.

## Common next steps

- Run `mlvamaps COMMAND --help` for command-specific options.
- [CLI options and thresholds](docs/reference/cli.md)
- [Input and panel formats](docs/reference/input-formats.md)
- [Output file reference](docs/reference/outputs.md)
- [Calling and profiles](docs/concepts/calling-and-profiles.md)
- [Dataset aggregation and MYOGA export](docs/workflows/myoga-export.md)

`mlvamaps` uses 32 threads by default. Pass `-t N` to set a limit or `-t 0` to
use all detected CPUs. Use `--quiet` to suppress progress messages.

## Development

```bash
conda env create -f environment.yml
conda activate mlvamaps
python -m pip install --no-deps -e .
pytest -q
```

The software is licensed under GPL-3.0-only. Please report problems through
[GitHub Issues](https://github.com/microbemarsh/mlvamaps/issues).
