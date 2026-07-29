# Input and panel formats

## Sequencing reads

Supported read suffixes:

- `.fastq`
- `.fq`
- `.fastq.gz`
- `.fq.gz`

FASTQ mode expects individual reads to contain a valid paired-primer product.
This directly supports amplicon sequencing. Accurate primer-spanning reads from
other assays can also be used, while non-spanning shotgun reads should be
supplied as assembly support instead.

## Input directories

`mlvamaps call` accepts a directory in place of one input file:

```bash
mlvamaps call primers.tsv sequence_files/ -o results
```

MLVAMaps processes supported FASTA and FASTQ files directly inside the
directory in filename order. Discovery is non-recursive, unrelated files are
ignored, and FASTA and FASTQ inputs may be mixed. Each input filename supplies
the sample ID and gets a separate `results/<sample-id>/` directory. Filenames
that collapse to the same sample ID, such as `sample.fasta` and
`sample.fastq`, are rejected to prevent output collisions.

`--sample-id`, `--reads`, and `--bam`/`--alignments` apply to single-file calls
and cannot be combined with directory input.

## Assemblies

Supported assembly suffixes include:

- `.fasta`
- `.fa`
- `.fna`
- `.fas`
- Gzip-compressed FASTA equivalents

## Minimal primer panel

A CSV, TSV, or legacy whitespace-delimited table can provide:

```text
locus_id
forward_primer
reverse_primer
repeat_unit_length_bp
expected_product_size_bp
nominal_repeat_units
```

The first three fields are required.

Example:

```text
locus_id	forward_primer	reverse_primer	repeat_unit_length_bp	expected_product_size_bp	nominal_repeat_units
vrrA_12bp_314bp_10U	CACAACTACCACCGATGGCACA	GCGCGTTTCGTTTGATTCATAC	12	314	10
```

Primer sequences are supplied 5-prime to 3-prime. MLVAMaps searches for the
reverse complement of the reverse primer in an oriented product. IUPAC
degenerate bases are handled by Amplirust.

## Rich locus table

A richer TSV can additionally define:

| Column | Purpose |
| --- | --- |
| `repeat_motif` | Expected VNTR motif used for read-pattern evidence. |
| `left_flank_sequence` | Sequence immediately left of the repeat, used to refine the boundary. |
| `right_flank_sequence` | Sequence immediately right of the repeat. |
| `expected_min_repeats` | Lower expected repeat-count bound. |
| `expected_max_repeats` | Upper expected repeat-count bound. |
| `expected_amplicon_min_bp` | Smallest accepted paired-primer product. |
| `expected_amplicon_max_bp` | Largest accepted paired-primer product. |
| `chrom_or_contig`, `start`, `end` | Optional reference annotation. |
| `pool_id` | Optional multiplex or primer-pool identifier. |

Examples are under `examples/`.

## Profile database

A profile database is a TSV with:

- Required `profile_id`.
- Optional `strain_id` and `metadata`.
- One column per comparable locus.

```text
profile_id	strain_id	VNTR_01	VNTR_02
P1	STRAIN_1	5	4
```

Profile locus names must match the panel. Empty values are skipped during
comparison. Profile rows from different MLVA schemes should not be mixed unless
their locus definitions and repeat-number conventions are compatible.

## Phylogenetic sequence database

`--database PATH` enables EPA-ng fixed-tree query placement independently of
the MLVA profile table. Prefer the top-level directory produced by
`mlvamaps build-reference` (or its `database/` subdirectory): MLVAMaps then
reuses the saved reference alignment, RAxML-NG tree, and selected model for
each locus. It runs MAFFT only to add the query without changing reference
coordinates, followed by EPA-ng placement.

Because each locus contains only one short query, MLVAMaps parallelizes EPA-ng
across callable loci. The `--threads` CPU budget is divided among the concurrent
placement processes; with at least as many loci as CPUs, each process gets one
native thread.

Sequence-only databases remain supported and build missing reference trees on
demand. Their recommended layout is one FASTA per locus in a directory;
each filename stem must exactly match the panel locus and each FASTA header is
the reference identifier:

```text
reference_sequences/
  VNTR_01.fasta
  VNTR_02.fasta
```

A long-form TSV with `reference_id`, `locus_id`, and `sequence` columns is also
accepted. A single FASTA can be used when every header contains exactly one
panel locus, for example `>REFERENCE_1|VNTR_01`.

Use the same reference identifiers across locus files. Only references present
at every successfully placed query locus are included in the summed-distance
ranking, preventing incomplete references from receiving artificially small
totals.

For repeat-aware phylogenetic analysis, rich panels should provide both flank
sequences and repeat-unit length or motif. MLVAMaps uses the flanks to remove
the tandem-repeat tract from the SNP-tree alignment while retaining its repeat
count and unit haplotype as explicit marker features. Motif-run detection is a
fallback when flanks are unavailable; sequences that cannot be bounded remain
unmasked and are labelled as such in the output.

An optional `--reference-metadata` TSV or CSV can associate references with
space, time, and source information:

```text
reference_id	collection_date	latitude	longitude	location	source
R1	2024-01-02	40.0	-75.0	Site A	environment
```

When the sequence database is a directory, a file named
`reference_metadata.tsv` is detected automatically.

Existing metadata used or produced by `geonome-flow reference-build` can be
converted directly:

```bash
python scripts/convert_geonome_metadata.py /path/to/geonome/metadata.tsv \
  --output reference_metadata.tsv
```

The input can be the canonical `metadata.tsv` passed to the workflow or its
surveillance `inputs/metadata.csv`. A Python-built reference directory,
`reference_manifest.json`, or `normalized_metadata.tsv` is also accepted. The
converter maps Geonome `genome_id`, `sample`, or `accession` identifiers to
`reference_id`, retains dates and coordinates, and maps location and sample
source fields to the columns above.

## Existing alignment support

Assembly mode accepts SAM or BAM aligned to the supplied assembly through
`--bam` or `--alignments`. Contig names must match the assembly contig names
reported by Amplirust.
