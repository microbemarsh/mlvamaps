# FASTQ and amplicon sequencing workflow

FASTQ calling requires a completed `--database`. Short and long reads map
directly to its observed reference loci and share molecule-level repeat
inference and reference classification. A reference match alone never supplies
the sample's repeat count.

```bash
mlvamaps call --database reference_build -i sample.fastq.gz -o results
```

For paired reads, see the [Illumina workflow](illumina.md).

## 1. Load the database panel and filter reads

The saved `reference_panel.tsv` supplies primers, locus names, repeat
calibration, and fingerprint order unless `-p` supplies an explicit panel.
The database's per-locus FASTAs supply the observed mapping targets.
Current-schema databases remain usable; their older `competitive_mapping/`
assets are retained for compatibility but do not supply the primary FASTQ
mapping targets or numeric allele estimates.

An optional Deacon target screen can precede QC. In metagenome mode, an available
database recruitment index is used automatically. Reads then pass the configured
length and mean-quality filters; original qualities remain attached.
The outputs include `qc_summary.tsv` and `filtered_reads.fastq.gz`.

## 2. Recover primer-bounded long-read products

Retained long reads are projected losslessly to FASTA and searched with the
existing Sassy-backed in silico PCR engine. Primer orientation, pairing,
IUPAC handling, and product-length constraints remain in use. Tied assignments
to different loci or product sizes do not become accepted primer assignments.
PCR searches share a Sassy process across batches of up to 64 records or 65,536
input bases; a larger record occupies its own batch. Each record and orientation
remains a separate search target, preserving molecule identities and coordinates.

Recovered products contribute direct molecule measurements even if minimap2
misses the read. This scan is not a gate that discards reads without both
primers: partial reads continue to observed-reference mapping.

Outputs include `filtered_reads.fasta.gz`, the `in_silico_pcr/` files, and
`read_locus_assignments.tsv` from the subsequent combined assignment stage.

## 3. Map to observed reference loci

Minimap2 maps retained reads competitively against observed database amplicons.
Identical sequences within a locus share a target while retaining all source
reference IDs. Alternative reference mappings are preserved. Synthetic repeat
expansions do not define the sample's allele hypotheses in this path.

For short reads, flank recruitment and mate rescue run during QC and feed the
same measurement stage. The appropriate technology preset is used for mapping.
The observed target FASTA, metadata, and provenance are recorded under
`candidate_mapping/`; minimap2 builds its target index in memory.

SAM is streamed through htslib to temporary BAM. The applicable intermediate or
debug flag retains `candidate_mapping/candidate_alignments.bam`; normal operation
removes it. Decoded alignments remain available to reference classification.

## 4. Measure molecules and call repeats

Complete primer products use the same calibrated product-length and rounding
functions as assemblies. Independently resolved repeat flanks can also measure
a repeat interval. Paired reads can contribute a uniquely merged fragment or a
supported local assembly, as described in the Illumina workflow.
For concrete DNA anchors in windows of at most 512 bases, an in-process edit-distance
existence check rejects impossible matches before launching Sassy. Accepted
approximate matches still use Sassy's coordinates and traceback.

Only measured evidence supplies a numeric allele. Flank-only mappings,
reference-relative alignment scores, and paired geometry establish presence
without filling in a count. Conflicting measurements of the same molecule remain
unresolved. Reference copies and mate records do not multiply molecule support.
Novel and fractional measured alleles remain callable outside the catalog range.

`common_locus_calls.tsv` distinguishes `presence_support` from informative
`molecule_support`. `molecule_candidate_evidence.tsv` retains reference IDs,
alignment diagnostics, evidence types, and inferred molecule counts.

Low-depth calls, mixtures, ambiguous measurements, detected-but-unresolved loci,
and missing loci retain distinct statuses. `PRESENT_COUNT_UNKNOWN` has no numeric
repeat count; missing values are never converted to zero.

## 5. Confirm representative sequences

The existing long-read recruitment, mapped-variant, SPOARS consensus, and SNP
workflow remains available for sequence confirmation and reporting. It retains
its own recruitment and representative mappings. Consensus is not required for
a directly measured molecule call, and a nearest-reference allele does not
overwrite the measured count.

Relevant outputs include `local_locus_products.fasta.gz`,
`local_assembly_concordance.tsv`, mapped-variant tables, and locus SNP evidence.

## 6. Classify against references using the shared evidence

Reference classification reuses the observed-reference alignments from step 3;
it does not map the FASTQ a second time. Sequence mismatches and measured repeat
length differences both inform reference likelihoods. A measured paired-fragment
count contributes once per molecule, even when both mates align. Where no measured
count exists, the existing alignment-derived repeat evidence remains available
to classification without becoming a numeric locus call.

Equivalent references remain grouped. Isolate support and metagenome mixture EM
retain an unclassified component, and weak or ambiguous evidence remains explicit.
Model support is conditional on the available reference catalog.

Outputs include `classification/molecule_reference_likelihoods.tsv.gz`, reference
matches, and taxonomic identification tables when applicable. Primer-only recovery
without an observed-reference alignment can support a locus count but does not
by itself establish a reference identity.

## 7. Export results

Start with `report.html`, `calls.tsv`, `mlva_fingerprint.tsv`, and
`locus_repeat_counts.tsv`. Profile matching and reporting use the measured locus
calls. Use the same panel calibration and `--assembly-round-tolerance` when
comparing FASTQ and assembly fingerprints.

Synthetic regression coverage includes single-end, paired, overlapping,
shotgun, and long reads, novel and fractional alleles, zero-repeat products,
primer indels, and primer rescue when mapping misses reads. Performance and
accuracy on a particular sample still require comparison with its known truth.

## Synthetic runtime check

On 2026-10-19, three alternating runs of the pre-optimization and updated workflows
used two threads and fresh processes on the same machine. Median calling times,
including output generation but excluding Python imports, were:

| Input | Before | After | Time reduction |
| --- | ---: | ---: | ---: |
| 11,100 pairs of 100-base reads | 7.60 s | 5.05 s | 34% |
| 132 long-read amplicons | 2.69 s | 1.16 s | 57% |

Both inputs used an 11-locus synthetic panel. The paired input contained 1,100
locus-associated pairs and 10,000 random background pairs; the long-read input
contained 12 copies of each product. The baseline already used observed-reference
mapping. All 28 short-read and 37 long-read output tables matched across versions
in each trial after normalizing output-directory paths and decompressing gzip
tables. These small synthetic timings measure reduced Sassy process overhead;
they do not establish throughput or accuracy on real server datasets.
