# The emu typing workflow

This branch uses alignment-based reference support for all `call --database`
inputs. It independently implements an Emu-inspired likelihood/mixture model;
it does not invoke or require the Emu package.

## Preserved outputs and assembly extraction

Assembly inputs still use Sassy-backed in silico PCR, the established product
selection and repeat-count conventions, and MLVA_finder-compatible tables.
Recovered products are aligned to corresponding observed reference loci using
Parasail. `classification/query_amplicons.fasta` and
`classification/assembly_reference_evidence.tsv` retain the query sequences,
reported repeat counts, reference counts, mismatch counts, and extended CIGARs.
Mismatch counts describe aligned differences; they are not a quality-filtered
SNP call set. Read SNP outputs retain their existing coverage/quality filters.

FASTQ allele inference retains competitive candidate mapping. Reference typing
uses a separate alignment pass against observed reference amplicons. Synthetic
candidate alleles and representative consensus sequences do not become
reference-classification observations.

## Removed typing options

The CLI no longer accepts `--phylogenetics`, `--phylogeny-snp-weight`,
`--phylogeny-repeat-weight`, the legacy target-calibration/bootstrap/placement
options, or `calibrate-taxa`. Reference building no longer creates locus
phylogenies or accepts `--min-references-per-tree`. MAFFT, RAxML-NG, EPA-ng and
MUMmer are no longer runtime dependencies. `export-myoga` retains repeat and
categorical profile comparisons; its combined-marker options are removed from
both the CLI and Python API. The legacy `phylogeny`, `taxon_assignment`, and
`combined_marker_export` modules have been deleted. Shared reference-loading
helpers now live in `reference_database`; repeat-profile Newick generation
lives in `profile_tree`.

Unused short-read schema-1 contexts, Amplirust aliases and executable options,
and the ignored `--min-cluster-size`, `--cluster-min-identity` and `--vsearch-bin`
options are also removed. Use the current Sassy PCR functions and versioned
competitive candidate contexts. MLVA_finder-compatible products and output
tables remain supported.

Use `--classification-repeat-scale` for the repeat discrepancy term and
`--taxon-min-loci` for the minimum observed loci needed for supported reference
classification. These controls are not equivalent to the old distance weights.
See [the model description](../concepts/mapping-classification.md).

Use `mapping_reference_matches` in Python result dictionaries. The old
`combined_marker_matches` alias is removed. Reference summaries use
`reference_status` instead of `tree_status`; taxon summaries no longer report
`trees_built`. A single observed reference can be retained, but a catalog with
only one reference cannot establish a supported unique identification.

## Existing databases and results

Schema 2.0 databases with the competitive mapping assets remain usable. Trees
in an older database are ignored. Earlier schemas must be rebuilt.
Use a fresh results directory for branch comparisons. Existing results are not
deleted or migrated. Manifest samples with a completed mapping classification
are skipped unless `--force` is supplied; use it when rerunning with new settings.
Reports never fall back to an older phylogenetic classification.

Repeat-profile trees remain downstream similarity visualizations. Neither EM
fractions nor query-to-reference likelihoods are pairwise evolutionary distances.
