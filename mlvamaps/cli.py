from __future__ import annotations

import argparse
from pathlib import Path

from .assembly_call import ASSEMBLY_ALGORITHMS, run_assembly_call
from .concurrency import DEFAULT_THREADS
from .in_silico_pcr import run_in_silico_pcr
from .io import open_text
from .pipeline import run_call
from .reference_builder import build_reference_database
from .reference_pipeline import (
    build_taxon_references,
    prepare_taxon_references,
    read_taxon_references,
)
from .simulation import simulate_reads


FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
FASTA_SUFFIXES = (
    ".fasta",
    ".fa",
    ".fna",
    ".fas",
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fas.gz",
)


def _sample_id_from_path(path: str) -> str:
    sample = Path(path).name
    for suffix in FASTQ_SUFFIXES + FASTA_SUFFIXES:
        if sample.lower().endswith(suffix):
            return sample[: -len(suffix)]
    return Path(sample).stem


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be at least 0")
    return parsed


def _fraction(value: str) -> float:
    parsed = float(value)
    if not 0.0 <= parsed <= 1.0:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return parsed


def _round_tolerance(value: str) -> float:
    parsed = float(value)
    if not 0.0 <= parsed <= 0.5:
        raise argparse.ArgumentTypeError("must be between 0 and 0.5")
    return parsed


def _nonnegative_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be at least 0")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def _input_kind(path: str) -> str:
    lower = path.lower()
    if lower.endswith(FASTQ_SUFFIXES):
        return "fastq"
    if lower.endswith(FASTA_SUFFIXES):
        return "fasta"
    with open(path) as handle:
        first = handle.read(1)
    if first == "@":
        return "fastq"
    if first == ">":
        return "fasta"
    raise ValueError(f"Could not tell whether {path!r} is FASTQ reads or FASTA assembly.")


def _input_files(path: str) -> list[Path]:
    """Resolve one input file or the supported sequence files in a directory."""
    input_path = Path(path)
    if input_path.is_file():
        return [input_path]
    if not input_path.exists():
        raise ValueError(f"Input path does not exist: {path}")
    if not input_path.is_dir():
        raise ValueError(f"Input path is not a file or directory: {path}")
    supported = FASTQ_SUFFIXES + FASTA_SUFFIXES
    files = sorted(
        (
            candidate
            for candidate in input_path.iterdir()
            if candidate.is_file()
            and candidate.name.lower().endswith(supported)
        ),
        key=lambda candidate: candidate.name.lower(),
    )
    if not files:
        suffixes = ", ".join(supported)
        raise ValueError(
            f"Input directory {path!r} contains no supported FASTA or FASTQ "
            f"files ({suffixes})."
        )
    return files


def _looks_like_loci_file(path: str) -> bool:
    with open_text(path, "rt") as handle:
        header = handle.readline().strip().lower().replace(",", "\t").split("\t")
    return "locus_id" in header and (
        "repeat_motif" in header
        or "left_flank_sequence" in header
        or "expected_min_repeats" in header
        or "chrom_or_contig" in header
    )


def _set_panel_path(args: argparse.Namespace, path: str) -> None:
    if _looks_like_loci_file(path):
        args.loci = path
    else:
        args.primers = path


def _resolve_call_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    paths = list(args.paths)
    if len(paths) > 2:
        parser.error(
            "call accepts at most two positional paths: primers.tsv and a "
            "FASTQ/FASTA file or input directory"
        )
    if len(paths) == 2:
        if not args.loci and not args.primers:
            _set_panel_path(args, paths[0])
        elif not args.input_path:
            parser.error("call got two positional paths plus --primers/--loci; pass only the sample path positionally")
        if not args.input_path:
            args.input_path = paths[1]
    elif len(paths) == 1:
        if args.loci or args.primers:
            if not args.input_path:
                args.input_path = paths[0]
        elif args.input_path:
            _set_panel_path(args, paths[0])
        else:
            parser.error("call needs both a primer file and an input file")

    if not args.input_path and args.reads_path:
        args.input_path = args.reads_path
        args.reads_path = None
    if not args.loci and not args.primers:
        parser.error(
            "call requires a primer file, for example: "
            "mlvamaps call primers.tsv sample.fastq.gz"
        )
    if not args.input_path:
        parser.error("call requires an input FASTQ/FASTA file or a directory containing them")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mlvamaps",
        description="Simple MLVA/VNTR calling from primers plus FASTQ, FASTA, or an input directory",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    call = subparsers.add_parser(
        "call",
        help="Call VNTRs from primers plus FASTQ/FASTA files or a directory",
        epilog=(
            "Examples:\n"
            "  mlvamaps call primers.tsv sample.fastq.gz\n"
            "  mlvamaps call primers.tsv assembly.fasta\n"
            "  mlvamaps call primers.tsv sequence_files/\n"
            "  mlvamaps call primers.tsv assembly.fasta --reads sample.fastq.gz\n"
            "  mlvamaps call primers.tsv assembly.fasta --bam assembly_reads.bam"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    call.add_argument(
        "paths",
        nargs="*",
        metavar="PATH",
        help="primers.tsv plus a FASTQ/FASTA file or directory",
    )
    call.add_argument(
        "--input",
        dest="input_path",
        metavar="PATH",
        help="FASTQ/FASTA file or directory containing supported sequence files",
    )
    call.add_argument("--reads", dest="reads_path", metavar="FASTQ", help="Reads to map for assembly depth support")
    call.add_argument("--bam", "--alignments", dest="alignments_path", metavar="BAM/SAM", help="Assembly-aligned BAM/SAM for assembly depth support")
    call.add_argument("--loci")
    call.add_argument("--primers", help="Primer-pair CSV/TSV/whitespace file with locus, forward, reverse columns")
    call.add_argument("--profiles")
    call.add_argument(
        "--database",
        help="Reference-build directory or per-locus sequence database for phylogenetic placement",
    )
    call.add_argument(
        "--reference-metadata",
        help="TSV/CSV with reference_id and optional date, coordinates, location, and source",
    )
    call.add_argument(
        "-o",
        "--output",
        "--outdir",
        dest="outdir",
        default="results",
        metavar="DIR",
        help="Output directory (default: %(default)s)",
    )
    call.add_argument("--sample-id")
    call.add_argument("--min-read-length", type=int, default=50)
    call.add_argument("--max-read-length", type=int, default=100000)
    call.add_argument(
        "--min-qscore",
        type=float,
        default=17.0,
        help=(
            "Minimum mean read Phred score (default: %(default)s, approximately "
            "98%% per-base accuracy)"
        ),
    )
    call.add_argument(
        "--recruitment-database",
        help=(
            "Reference-build directory supplying canonical locus products for "
            "FASTQ recruitment without enabling phylogenetic placement"
        ),
    )
    call.add_argument("--max-primer-mismatches", type=int, default=2)
    call.add_argument(
        "--algorithm",
        "--assembly-algorithm",
        choices=ASSEMBLY_ALGORITHMS,
        default="legacy",
        help=(
            "FASTA allele caller: exact historical MLVA_finder rules or the "
            "depth-aware probabilistic caller (default: %(default)s)"
        ),
    )
    call.add_argument(
        "--assembly-round-tolerance",
        type=_round_tolerance,
        default=0.25,
        metavar="FRACTION",
        help="Legacy integer-rounding tolerance (default: %(default)s)",
    )
    call.add_argument(
        "--read-calling-convention",
        choices=("assembly", "probabilistic"),
        default="assembly",
        help=(
            "FASTQ allele convention: assembly uses the same calibrated "
            "product-length rounding as FASTA calls; probabilistic retains "
            "raw half-unit inference (default: %(default)s)"
        ),
    )
    call.add_argument(
        "--sample-mode",
        choices=("isolate", "metagenome"),
        default="metagenome",
        help="Interpret FASTQ evidence as an isolate or metagenome (default: %(default)s)",
    )
    call.add_argument(
        "--fastq-strategy",
        choices=("recruit", "primer"),
        default="recruit",
        help=(
            "FASTQ discovery strategy: competitive locus recruitment with "
            "primer fallback, or legacy primer-only assignment "
            "(default: %(default)s)"
        ),
    )
    call.add_argument(
        "--recruitment-preset",
        help="Optional minimap2 -x preset for competitive long-read recruitment",
    )
    call.add_argument(
        "--recruitment-min-identity",
        type=_fraction,
        default=0.9,
        help="Minimum recruitment alignment identity (default: %(default)s)",
    )
    call.add_argument(
        "--recruitment-min-aligned-bp",
        type=_positive_int,
        default=100,
        help="Minimum aligned bases for locus presence (default: %(default)s)",
    )
    call.add_argument(
        "--recruitment-min-locus-margin",
        type=_nonnegative_int,
        default=10,
        help=(
            "Minimum alignment-score lead over the next-best locus "
            "(default: %(default)s)"
        ),
    )
    call.add_argument(
        "--taxon-screen-index",
        metavar="DEACON_IDX",
        help=(
            "Before MLVA analysis, retain only reads matching this target-taxon "
            "Deacon pangenome index (see github.com/bede/deacon-indexes)"
        ),
    )
    call.add_argument(
        "--taxon-screen-abs-threshold",
        type=_positive_int,
        default=2,
        help="Minimum shared Deacon minimizers for target retention (default: %(default)s)",
    )
    call.add_argument(
        "--taxon-screen-rel-threshold",
        type=_fraction,
        default=0.01,
        help="Minimum relative Deacon minimizer match fraction (default: %(default)s)",
    )
    call.add_argument(
        "--deacon-bin",
        default="deacon",
        metavar="PATH",
        help="Deacon executable used by --taxon-screen-index (default: %(default)s)",
    )
    call.add_argument(
        "--min-depth",
        type=_positive_int,
        default=1,
        help=(
            "Minimum informative reads required to avoid LOW_DEPTH "
            "(default: %(default)s)"
        ),
    )
    call.add_argument("--min-posterior", type=float, default=0.75)
    call.add_argument(
        "--repeat-range-tolerance",
        type=_nonnegative_float,
        default=1.0,
        metavar="REPEATS",
        help=(
            "Allowed repeat-count distance beyond expected locus bounds before "
            "setting OUT_OF_RANGE (default: %(default)s)"
        ),
    )
    call.add_argument(
        "--max-confidence-depth",
        type=_positive_float,
        default=25.0,
        help=(
            "Cap on effective primary-read evidence used to sharpen allele "
            "confidence (default: %(default)s)"
        ),
    )
    call.add_argument(
        "--min-cluster-size",
        type=_positive_int,
        default=1,
        help="Deprecated compatibility option; mapping groups retain low-depth evidence",
    )
    call.add_argument(
        "--cluster-min-identity",
        type=_fraction,
        default=0.97,
        help="Deprecated compatibility option; sequence clustering is no longer used",
    )
    call.add_argument(
        "--min-mixture-fraction",
        type=_fraction,
        default=0.01,
        help="Minimum EM-estimated fraction for a meaningful variant (default: %(default)s)",
    )
    call.add_argument(
        "--min-secondary-reads",
        type=_positive_int,
        default=2,
        help=(
            "Minimum reads required to promote a secondary variant from "
            "candidate to confirmed (default: %(default)s)"
        ),
    )
    call.add_argument(
        "--vsearch-bin",
        default="vsearch",
        metavar="PATH",
        help="Deprecated compatibility option; FASTQ grouping now uses read mapping",
    )
    call.add_argument(
        "--amplirust-bin",
        default="amplirust",
        metavar="PATH",
        help=argparse.SUPPRESS,
    )
    call.add_argument(
        "--minimap2-bin",
        default="minimap2",
        metavar="PATH",
        help="minimap2 executable for representative and assembly-support mapping (default: %(default)s)",
    )
    call.add_argument(
        "--mafft-bin",
        default="mafft",
        metavar="PATH",
        help="MAFFT executable for optional per-locus phylogenetic placement (default: %(default)s)",
    )
    call.add_argument(
        "--raxml-ng-bin",
        default="raxml-ng",
        metavar="PATH",
        help="RAxML-NG executable for maximum-likelihood locus trees (default: %(default)s)",
    )
    call.add_argument(
        "--epa-ng-bin",
        default="epa-ng",
        metavar="PATH",
        help="EPA-ng executable for fixed-tree query placement (default: %(default)s)",
    )
    call.add_argument(
        "--dnadiff-bin",
        default="dnadiff",
        metavar="PATH",
        help="MUMmer dnadiff executable for assembly whole-genome SNP tie breaking (default: %(default)s)",
    )
    call.add_argument(
        "--raxml-model",
        default="DNA",
        metavar="MODEL",
        help="RAxML-NG nucleotide model for locus trees (default: %(default)s)",
    )
    call.add_argument(
        "--phylogeny-snp-weight",
        type=_nonnegative_float,
        default=1.0,
        help="Weight for normalized SNP-tree distance in combined marker ranking (default: %(default)s)",
    )
    call.add_argument(
        "--phylogeny-repeat-weight",
        type=_nonnegative_float,
        default=1.0,
        help="Weight for normalized tandem-repeat distance in combined marker ranking (default: %(default)s)",
    )
    call.add_argument(
        "--no-locus-mapping",
        action="store_true",
        help="Skip minimap2 representative mapping and SNP evidence generation",
    )
    call.add_argument(
        "--min-mapping-quality",
        type=_nonnegative_int,
        default=0,
        help="Minimum minimap2 MAPQ used for locus mapping evidence (default: %(default)s)",
    )
    call.add_argument(
        "--min-base-quality",
        type=_nonnegative_int,
        default=20,
        help="Minimum base quality used for coverage and SNP evidence (default: %(default)s)",
    )
    call.add_argument(
        "--min-snp-depth",
        type=_positive_int,
        default=3,
        help="Minimum quality-filtered depth for a SNP call (default: %(default)s)",
    )
    call.add_argument(
        "--min-snp-alternate-reads",
        type=_positive_int,
        default=2,
        help="Minimum reads supporting a non-reference SNP allele (default: %(default)s)",
    )
    call.add_argument(
        "--min-snp-frequency",
        type=_fraction,
        default=0.2,
        help="Minimum non-reference allele frequency for a SNP call (default: %(default)s)",
    )
    call.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Worker threads (default: %(default)s; 0 uses all available CPUs)",
    )
    call.add_argument(
        "--minimap2-preset",
        help="Optional minimap2 -x preset for assembly read-depth mapping",
    )
    call.add_argument("--quiet", action="store_true", help="Suppress live progress updates")
    call.add_argument(
        "--debug-disagreements",
        action="store_true",
        help="Write read- and locus-level mapping versus measurement audit TSVs",
    )

    simulate = subparsers.add_parser("simulate", help="Simulate amplicon reads for a VNTR panel")
    simulate.add_argument("--loci", required=True)
    simulate.add_argument("--profile", dest="profiles")
    simulate.add_argument("--profile-id")
    simulate.add_argument("--sample-id", required=True)
    simulate.add_argument("--depth", type=int, default=200)
    simulate.add_argument("--error-rate", type=float, default=0.03)
    simulate.add_argument("--seed", type=int, default=13)
    simulate.add_argument(
        "-o", "--output", "--outdir", dest="outdir", required=True, metavar="DIR", help="Output directory"
    )

    extract = subparsers.add_parser(
        "extract-amplicons",
        help="Extract MLVA_finder-compatible amplicons from FASTA with Sassy",
    )
    extract.add_argument("--input", required=True, help="Input FASTA, optionally gzip-compressed")
    extract.add_argument("--loci")
    extract.add_argument("--primers", help="Primer-pair CSV/TSV/whitespace file with locus, forward, reverse columns")
    extract.add_argument(
        "-o",
        "--output",
        "--outdir",
        dest="outdir",
        default="assembly_amplicons",
        metavar="DIR",
        help="Output directory (default: %(default)s)",
    )
    extract.add_argument("--max-errors", type=int, default=2)
    extract.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="PCR search threads (default: %(default)s; 0 auto-detects CPUs)",
    )
    extract.add_argument("--circular", action="store_true")
    extract.add_argument("--no-search-rc", action="store_true")
    extract.add_argument("--trim-primers", action="store_true")
    extract.add_argument("--amplirust-bin", default="amplirust", help=argparse.SUPPRESS)

    prepare_reference = subparsers.add_parser(
        "prepare-reference",
        help="Download reproducible NCBI assembly inputs for one taxid or a CSV of taxids",
    )
    prepare_source = prepare_reference.add_mutually_exclusive_group(required=True)
    prepare_source.add_argument("--taxid", help="Single NCBI taxonomy identifier")
    prepare_source.add_argument(
        "--taxids-csv",
        help="CSV/TSV with a taxid column and optional name column",
    )
    prepare_reference.add_argument(
        "-o", "--output", "--outdir", dest="outdir", default="reference_builds"
    )
    prepare_reference.add_argument(
        "--assembly-source",
        choices=("refseq", "genbank", "all"),
        default="refseq",
        help="NCBI assembly source (default: %(default)s)",
    )
    prepare_reference.add_argument(
        "--datasets-arg",
        action="append",
        default=[],
        help="Extra NCBI Datasets argument; repeat as needed",
    )
    prepare_reference.add_argument(
        "--resume",
        action="store_true",
        help="Reuse each existing prepared/ncbi_dataset.zip",
    )
    prepare_reference.add_argument("--datasets-bin", default="datasets", help=argparse.SUPPRESS)
    prepare_reference.add_argument("--dataformat-bin", default="dataformat", help=argparse.SUPPRESS)

    reference = subparsers.add_parser(
        "build-reference",
        help="Build reference databases from local assemblies, one taxid, or a CSV of taxids",
    )
    reference_source = reference.add_mutually_exclusive_group(required=True)
    reference_source.add_argument(
        "--assemblies", help="Directory containing reference FASTA assemblies"
    )
    reference_source.add_argument("--taxid", help="Single NCBI taxonomy identifier")
    reference_source.add_argument(
        "--taxids-csv",
        help="CSV/TSV with a taxid column and optional name column",
    )
    panel = reference.add_mutually_exclusive_group(required=True)
    panel.add_argument("--primers", help="Primer-pair CSV/TSV with locus, forward, and reverse columns")
    panel.add_argument("--loci", help="Rich loci TSV (recommended when repeat motif/flanks are known)")
    reference.add_argument(
        "--metadata",
        help="Reference metadata CSV/TSV; required with --assemblies",
    )
    reference.add_argument("-o", "--output", "--outdir", dest="outdir", default="reference_build")
    reference.add_argument(
        "--assembly-source",
        choices=("refseq", "genbank", "all"),
        default="refseq",
        help="NCBI assembly source for taxid inputs (default: %(default)s)",
    )
    reference.add_argument(
        "--datasets-arg",
        action="append",
        default=[],
        help="Extra NCBI Datasets argument for taxid inputs; repeat as needed",
    )
    reference.add_argument(
        "--resume",
        action="store_true",
        help="Reuse each existing prepared/ncbi_dataset.zip for taxid inputs",
    )
    reference.add_argument("--datasets-bin", default="datasets", help=argparse.SUPPRESS)
    reference.add_argument("--dataformat-bin", default="dataformat", help=argparse.SUPPRESS)
    reference.add_argument("--max-primer-mismatches", type=_nonnegative_int, default=2)
    reference.add_argument(
        "--multiple-products",
        choices=("exclude", "best", "error"),
        default="exclude",
        help="Policy for assembly/locus pairs with multiple products (default: %(default)s)",
    )
    reference.add_argument(
        "--min-references-per-tree",
        type=_positive_int,
        default=3,
        help="Minimum extracted references required to infer a locus tree (default: %(default)s)",
    )
    reference.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Parallel extraction workers and tree-tool threads (default: %(default)s; 0 auto-detects CPUs)",
    )
    reference.add_argument("--quiet", action="store_true", help="Suppress live progress updates")
    reference.add_argument("--amplirust-bin", default="amplirust", help=argparse.SUPPRESS)
    reference.add_argument("--mafft-bin", default="mafft")
    reference.add_argument("--raxml-ng-bin", default="raxml-ng")
    reference.add_argument(
        "--raxml-model",
        default="DNA",
        help="RAxML-NG nucleotide model or model-selection set (default: %(default)s)",
    )
    return parser


def _run_single_input(
    args: argparse.Namespace,
    input_path: Path,
    outdir: Path,
    sample_id: str,
) -> None:
    if _input_kind(str(input_path)) == "fastq":
        result = run_call(
            reads_path=str(input_path),
            loci_path=args.loci,
            primers_path=args.primers,
            profiles_path=args.profiles,
            database_path=args.database,
            recruitment_database_path=args.recruitment_database,
            outdir=str(outdir),
            sample_id=sample_id,
            min_read_length=args.min_read_length,
            max_read_length=args.max_read_length,
            min_qscore=args.min_qscore,
            max_primer_mismatches=args.max_primer_mismatches,
            min_depth=args.min_depth,
            min_posterior=args.min_posterior,
            repeat_range_tolerance=args.repeat_range_tolerance,
            min_cluster_size=args.min_cluster_size,
            cluster_min_identity=args.cluster_min_identity,
            min_mixture_fraction=args.min_mixture_fraction,
            min_secondary_reads=args.min_secondary_reads,
            vsearch_bin=args.vsearch_bin,
            amplirust_bin=args.amplirust_bin,
            minimap2_bin=args.minimap2_bin,
            mafft_bin=args.mafft_bin,
            raxml_ng_bin=args.raxml_ng_bin,
            epa_ng_bin=args.epa_ng_bin,
            raxml_model=args.raxml_model,
            phylogeny_snp_weight=args.phylogeny_snp_weight,
            phylogeny_repeat_weight=args.phylogeny_repeat_weight,
            reference_metadata_path=args.reference_metadata,
            locus_mapping=not args.no_locus_mapping,
            min_mapping_quality=args.min_mapping_quality,
            min_base_quality=args.min_base_quality,
            min_snp_depth=args.min_snp_depth,
            min_snp_alternate_reads=args.min_snp_alternate_reads,
            min_snp_frequency=args.min_snp_frequency,
            threads=args.threads,
            show_progress=not args.quiet,
            sample_mode=args.sample_mode,
            assembly_equivalent_reads=args.read_calling_convention == "assembly",
            assembly_round_tolerance=args.assembly_round_tolerance,
            max_confidence_depth=args.max_confidence_depth,
            fastq_strategy=args.fastq_strategy,
            recruitment_preset=args.recruitment_preset,
            recruitment_min_identity=args.recruitment_min_identity,
            recruitment_min_aligned_bp=args.recruitment_min_aligned_bp,
            recruitment_min_locus_margin=args.recruitment_min_locus_margin,
            debug_disagreements=args.debug_disagreements,
            taxon_screen_index=args.taxon_screen_index,
            taxon_screen_abs_threshold=args.taxon_screen_abs_threshold,
            taxon_screen_rel_threshold=args.taxon_screen_rel_threshold,
            deacon_bin=args.deacon_bin,
        )
        print(f"Wrote easy MLVA calls to {result['calls']}")
        print(f"Wrote detailed allele evidence to {result['allele_calls']}")
        print(f"Wrote individual locus repeat counts to {result['repeat_counts']}")
        print(f"Wrote mapped VNTR variant groups to {result['mapped_variant_table']}")
        print(f"Wrote EM variant abundance estimates to {result['mixture_abundance']}")
        print(f"Wrote mapped read-group evidence to {result['mapped_read_memberships']}")
        if args.profiles or args.database:
            print(f"Wrote ranked profile matches to {result['profile_matches']}")
            print(f"Wrote per-locus profile comparisons to {result['profile_match_loci']}")
        if not args.no_locus_mapping:
            print(f"Wrote locus mapping summaries to {result['mapping_summary']}")
            print(f"Wrote locus SNP evidence to {result['mapping_snps']}")
        print(f"Wrote report to {result['report']}")
        if args.database:
            print(f"Wrote per-locus trees to {result['phylogeny']}")
            print(f"Wrote phylogenetic matches to {result['phylogenetic_matches']}")
            print(f"Wrote combined repeat/SNP matches to {result['combined_marker_matches']}")
            print(f"Wrote MYOGA-compatible tree to {result['combined_marker_tree']}")
        return

    result = run_assembly_call(
        assembly_path=str(input_path),
        loci_path=args.loci,
        primers_path=args.primers,
        outdir=str(outdir),
        sample_id=sample_id,
        reads_path=args.reads_path,
        alignments_path=args.alignments_path,
        profiles_path=args.profiles,
        database_path=args.database,
        max_primer_mismatches=args.max_primer_mismatches,
        assembly_round_tolerance=args.assembly_round_tolerance,
        algorithm=args.algorithm,
        min_posterior=args.min_posterior,
        threads=args.threads,
        minimap2_preset=args.minimap2_preset,
        minimap2_bin=args.minimap2_bin,
        amplirust_bin=args.amplirust_bin,
        mafft_bin=args.mafft_bin,
        raxml_ng_bin=args.raxml_ng_bin,
        epa_ng_bin=args.epa_ng_bin,
        dnadiff_bin=args.dnadiff_bin,
        raxml_model=args.raxml_model,
        phylogeny_snp_weight=args.phylogeny_snp_weight,
        phylogeny_repeat_weight=args.phylogeny_repeat_weight,
        reference_metadata_path=args.reference_metadata,
        show_progress=not args.quiet,
    )
    print(f"Wrote easy MLVA calls to {result['calls']}")
    print(f"Wrote individual locus repeat counts to {result['repeat_counts']}")
    print(f"Wrote assembly amplicons to {result['amplicons']}")
    if args.profiles or args.database:
        print(f"Wrote ranked profile matches to {result['profile_matches']}")
        print(f"Wrote per-locus profile comparisons to {result['profile_match_loci']}")
    print(f"Wrote report to {result['report']}")
    if args.database:
        print(f"Wrote per-locus trees to {result['phylogeny']}")
        print(f"Wrote phylogenetic matches to {result['phylogenetic_matches']}")
        print(f"Wrote combined repeat/SNP matches to {result['combined_marker_matches']}")
        print(f"Wrote MYOGA-compatible tree to {result['combined_marker_tree']}")
    if args.reads_path or args.alignments_path:
        print(f"Wrote read-depth support to {result['read_support']}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "call":
        _resolve_call_args(parser, args)
        if args.reads_path and args.alignments_path:
            parser.error("call accepts either --reads or --bam/--alignments for assembly depth support, not both")
        try:
            input_files = _input_files(args.input_path)
        except ValueError as exc:
            parser.error(str(exc))
        batch = Path(args.input_path).is_dir()
        if batch and args.sample_id:
            parser.error("--sample-id cannot be used with an input directory; sample IDs come from filenames")
        if batch and (args.reads_path or args.alignments_path):
            parser.error("--reads and --bam/--alignments can only be used with a single input file")
        sample_ids = [_sample_id_from_path(str(path)) for path in input_files]
        duplicate_ids = sorted(
            sample_id for sample_id in set(sample_ids) if sample_ids.count(sample_id) > 1
        )
        if duplicate_ids:
            parser.error(
                "input filenames produce duplicate sample IDs: "
                + ", ".join(duplicate_ids)
            )
        for input_path, derived_sample_id in zip(input_files, sample_ids):
            sample_id = args.sample_id or derived_sample_id
            outdir = Path(args.outdir) / sample_id if batch else Path(args.outdir)
            if batch:
                print(f"Processing {input_path} as sample {sample_id}")
            _run_single_input(args, input_path, outdir, sample_id)
        return 0
    if args.command == "simulate":
        result = simulate_reads(
            loci_path=args.loci,
            profiles_path=args.profiles,
            profile_id=args.profile_id,
            outdir=args.outdir,
            sample_id=args.sample_id,
            depth=args.depth,
            error_rate=args.error_rate,
            seed=args.seed,
        )
        print(f"Wrote simulated reads to {result['reads']}")
        print(f"Wrote truth profile to {result['truth']}")
        return 0
    if args.command == "extract-amplicons":
        if not args.loci and not args.primers:
            parser.error("extract-amplicons requires either --loci or --primers")
        result = run_in_silico_pcr(
            input_path=args.input,
            loci_path=args.loci,
            primers_path=args.primers,
            outdir=args.outdir,
            max_errors=args.max_errors,
            threads=args.threads,
            circular=args.circular,
            search_rc=not args.no_search_rc,
            trim_primers=args.trim_primers,
        )
        print(f"Wrote normalized primer CSV to {result['primers']}")
        print(f"Wrote extracted amplicons to {result['products']}")
        print(f"Wrote primer-match stats to {result['stats']}")
        return 0
    if args.command == "prepare-reference":
        references = read_taxon_references(
            taxid=args.taxid, taxids_csv=args.taxids_csv
        )
        results = prepare_taxon_references(
            references,
            args.outdir,
            assembly_source=args.assembly_source,
            datasets_args=args.datasets_arg,
            datasets_bin=args.datasets_bin,
            dataformat_bin=args.dataformat_bin,
            resume=args.resume,
        )
        for result in results:
            print(
                f"Prepared taxid {result['taxid']} ({result['name']}) in "
                f"{result['outdir']}"
            )
        return 0
    if args.command == "build-reference":
        if args.assemblies and not args.metadata:
            parser.error("build-reference with --assemblies also requires --metadata")
        if not args.assemblies and args.metadata:
            parser.error("--metadata is only accepted with --assemblies")
        if not args.assemblies:
            references = read_taxon_references(
                taxid=args.taxid, taxids_csv=args.taxids_csv
            )
            result = build_taxon_references(
                references,
                primers_path=args.primers or args.loci,
                loci_path=args.loci,
                outdir=args.outdir,
                assembly_source=args.assembly_source,
                datasets_args=args.datasets_arg,
                datasets_bin=args.datasets_bin,
                dataformat_bin=args.dataformat_bin,
                resume=args.resume,
                multiple_products=args.multiple_products,
                max_primer_mismatches=args.max_primer_mismatches,
                min_references_per_tree=args.min_references_per_tree,
                threads=args.threads,
                amplirust_bin=args.amplirust_bin,
                mafft_bin=args.mafft_bin,
                raxml_ng_bin=args.raxml_ng_bin,
                raxml_model=args.raxml_model,
                show_progress=not args.quiet,
            )
            for reference in result["references"]:
                print(
                    f"Built taxid {reference['taxid']} ({reference['name']}) "
                    f"database at {reference['database']}"
                )
            print(f"Wrote taxid pipeline manifest to {result['manifest']}")
            return 0
        result = build_reference_database(
            assemblies_dir=args.assemblies,
            primers_path=args.primers or args.loci,
            loci_path=args.loci,
            metadata_path=args.metadata,
            outdir=args.outdir,
            multiple_products=args.multiple_products,
            max_primer_mismatches=args.max_primer_mismatches,
            min_references_per_tree=args.min_references_per_tree,
            threads=args.threads,
            amplirust_bin=args.amplirust_bin,
            mafft_bin=args.mafft_bin,
            raxml_ng_bin=args.raxml_ng_bin,
            raxml_model=args.raxml_model,
            show_progress=not args.quiet,
        )
        print(f"Wrote per-locus reference database to {result['database']}")
        print(f"Wrote reference build QC to {result['manifest']}")
        print(f"Wrote per-locus reference trees to {result['phylogeny']}")
        print(f"Wrote MYOGA metadata to {result['myoga_metadata']}")
        return 0
    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
