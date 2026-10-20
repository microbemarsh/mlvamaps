# FASTQ and amplicon sequencing workflow

FASTQ calling requires a completed `--database`. Short and long reads map
directly to its observed reference loci and share molecule-level repeat
inference and reference classification. Directly measured counts take precedence;
unresolved sizing can fall back to an explicitly labelled low-confidence estimate.

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

Complete primer products use the same calibrated product-length functions as
assemblies. Raw lengths and measurement uncertainty are retained through
inference; assembly rounding is applied to the final allele. Isolate mode fits
a single-allele likelihood with a broad error component for occasional outliers.
Metagenome mode uses the existing weighted EM solver to estimate per-locus
allele fractions. A secondary allele requires sufficient fraction and effective
molecule support, convergence, and a BIC-penalized likelihood improvement over
the single-allele model. EM convergence
and fractions are reported separately from single-allele model support.

Identical observations share likelihood calculations while retaining molecule
multiplicity. Long-read anchor measurements are cached by locus, sequence and
quality, including failed measurements, rather than repeated for each reference
allele. These caches are bounded and local to one inference run.

Independently resolved repeat flanks can also measure
a repeat interval. Paired reads can contribute a uniquely merged fragment or a
supported local assembly, as described in the Illumina workflow.
For concrete DNA anchors in windows of at most 512 bases, an in-process edit-distance
existence check rejects impossible matches before launching Sassy. Accepted
approximate matches still use Sassy's coordinates and traceback.

Measured evidence supplies direct alleles. If direct sizing remains unresolved,
report the best supported approximation as `ESTIMATED`: prefer competing measured
observations, then repeat-spanning CIGAR projections, then reference-conditioned
partial-alignment support. Alignment estimates require locus-flank support and
remain low confidence even when only one catalog allele is available. Reference
copies and mate records do not multiply molecule support. Estimates are kept
separate from measured counts used in reference classification.
Novel and fractional measured alleles remain callable outside the catalog range.
Optional independently calibrated paired-fragment lengths provide another
labelled estimate when no direct measurement is available; see the
[Illumina guide](illumina.md).

`common_locus_calls.tsv` distinguishes `presence_support` from informative
`molecule_support`. `molecule_candidate_evidence.tsv` retains reference IDs,
alignment diagnostics, evidence types, and inferred molecule counts.

Low-depth calls, mixtures, ambiguous measurements, detected-but-unresolved loci,
and missing loci retain distinct statuses. `PRESENT_COUNT_UNKNOWN` has no numeric
repeat count when neither sizing nor estimation is possible. The report shows
every locus in an always-visible table; missing values are never converted to zero.

## 5. Confirm representative sequences

The existing long-read recruitment, mapped-variant, SPOARS consensus, and SNP
workflow remains available for sequence confirmation and reporting. It retains
its own recruitment and representative mappings. Consensus is not required for
a directly measured molecule call, and a nearest-reference allele does not
overwrite the measured count.

A dominant-cluster consensus with sufficient input support can resolve an
uncertain call when its allele has molecule model support within fourfold of
the best allele. The result remains `ESTIMATED`, with method
`consensus_supported`. Consensus proposes a candidate or resolves ambiguity;
it never adds independent depth or changes the molecule measurements supplied
to reference classification. Strong calls and supported mixtures are retained.

Relevant outputs include `local_locus_products.fasta.gz`,
`local_assembly_concordance.tsv`, mapped-variant tables, and locus SNP evidence.

## 6. Classify against references using the shared evidence

Reference classification reuses the observed-reference alignments from step 3;
it does not map the FASTQ a second time. Sequence mismatches and measured repeat
length differences both inform reference likelihoods. A measured paired-fragment
count contributes once per molecule, even when both mates align. Where no measured
count exists, the existing alignment-derived repeat evidence remains available
to classification. A reported estimate does not become an independent measured
count in that classification model.

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

## Raw-length inference runtime check

On 2026-10-20, the previous caller and the raw-length/EM implementation were
compared on the same machine with two threads and five runs per workflow.
The benchmark times FASTQ calling and output generation, excluding database
construction and assembly truth generation.

| Synthetic input | Previous median | Updated median | Time reduction |
| --- | ---: | ---: | ---: |
| 4,400 full-product pairs | 0.723 s | 0.722 s | Essentially unchanged |
| 29,000 shotgun pairs | 3.751 s | 3.507 s | 6.5% |
| 4,400 long-read products | 14.643 s | 10.196 s | 30.4% |

All resulting fingerprints matched the assembly truth. Inputs use an 11-locus
panel covering novel, fractional, zero-repeat and primer-indel cases, with
replicated high-quality reads. These timings establish performance for this
synthetic isolate workload, not real noisy reads or mixtures.

Reproduce with dev dependencies and minimap2/Sassy on PATH:

```bash
python scripts/benchmark_fastq_inference.py --repetitions 5 --copies 100
```

The [benchmark record](../reference/fastq-likelihood-benchmark.json) contains
individual timings, the baseline commit, environment and input conditions.
Regression tests additionally cover noisy raw-length aggregation, allele
mixtures, duplicate-reference invariance, uncertain consensus and calibrated
paired-fragment estimates.
