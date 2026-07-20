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

`--database PATH` enables MAFFT alignment, RAxML-NG maximum-likelihood reference
trees, and EPA-ng fixed-tree query placement independently of the MLVA profile
table. The recommended layout is one FASTA per locus in a directory;
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

## Existing alignment support

Assembly mode accepts SAM or BAM aligned to the supplied assembly through
`--bam` or `--alignments`. Contig names must match the assembly contig names
reported by Amplirust.
