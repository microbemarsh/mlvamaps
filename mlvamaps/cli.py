from __future__ import annotations

import argparse
import csv
import copy
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from . import __version__
from .assembly_call import run_assembly_call
from .concurrency import DEFAULT_THREADS, resolve_threads
from .io import open_text, write_tsv
from .myoga_export import export_myoga
from .pipeline import run_call
from .reference_builder import build_reference_database
from .reference_pipeline import (
    build_taxon_references,
    read_taxon_references,
)
from .sample_metadata import MYOGA_SAMPLE_FIELDS, metadata_by_sample, read_sample_metadata, write_csv
from .short_reads import SAMPLE_SUMMARY_FIELDS, run_short_read_call


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
BATCH_OUTPUT_DIRECTORY = "batch_summary"


def _batch_allocation(total_threads: int, job_count: int) -> tuple[int, int]:
    """Return bounded sample workers and native threads per active sample."""
    total = resolve_threads(total_threads)
    memory_limit = int(os.environ.get("MLVAMAPS_MAX_CONCURRENT_SAMPLES", "4"))
    workers = max(1, min(job_count, total, max(1, memory_limit)))
    return workers, max(1, total // workers)


def _args_with_threads(args: argparse.Namespace, threads: int) -> argparse.Namespace:
    cloned = copy.copy(args)
    cloned.threads = threads
    return cloned


def _sample_output_dir(output_root: Path, sample_id: str) -> Path:
    """Return the output directory for one sample in a batch."""
    return output_root / sample_id


def _batch_output_dir(output_root: Path) -> Path:
    """Return the directory reserved for multi-sample aggregate outputs."""
    return output_root / BATCH_OUTPUT_DIRECTORY


def _check_batch_sample_ids(sample_ids: list[str]) -> None:
    """Prevent sample directories from colliding with batch-level outputs."""
    reserved = [
        sample_id
        for sample_id in sample_ids
        if sample_id.casefold() == BATCH_OUTPUT_DIRECTORY.casefold()
    ]
    if reserved:
        raise ValueError(
            f"sample ID {reserved[0]!r} is reserved for batch aggregate outputs"
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


def _short_read_directory_rows(path: str | Path) -> list[dict[str, str]]:
    """Discover exact PREFIX_1/2.fastq.gz mate pairs in one directory."""
    input_dir = Path(path)
    if not input_dir.exists():
        raise ValueError(f"Input path does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise ValueError("--short-reads requires -i to name a directory")
    suffixes = {"1": "_1.fastq.gz", "2": "_2.fastq.gz"}
    mates: dict[str, dict[str, Path]] = {}
    for candidate in input_dir.iterdir():
        if not candidate.is_file():
            continue
        for mate, suffix in suffixes.items():
            if candidate.name.endswith(suffix):
                prefix = candidate.name[: -len(suffix)]
                if prefix in {"", ".", ".."}:
                    raise ValueError(
                        f"Invalid short-read sample prefix in {candidate.name!r}"
                    )
                mates.setdefault(prefix, {})[mate] = candidate.resolve()
                break
    if not mates:
        raise ValueError(
            f"Input directory {str(input_dir)!r} contains no "
            "PREFIX_1.fastq.gz/PREFIX_2.fastq.gz pairs"
        )
    incomplete = [
        prefix for prefix, pair in sorted(mates.items()) if set(pair) != {"1", "2"}
    ]
    if incomplete:
        details = ", ".join(
            f"{prefix} (missing mate {'2' if '1' in mates[prefix] else '1'})"
            for prefix in incomplete
        )
        raise ValueError(f"Unpaired short-read files: {details}")
    return [
        {
            "sample_id": prefix,
            "reads1": str(mates[prefix]["1"]),
            "reads2": str(mates[prefix]["2"]),
        }
        for prefix in sorted(mates, key=str.casefold)
    ]


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


def _resolve_panel_option(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    panel_path = getattr(args, "panel_path", None)
    if not panel_path:
        return
    if getattr(args, "loci", None) or getattr(args, "primers", None):
        parser.error("provide exactly one panel with -p/--panel")
    if not Path(panel_path).is_file():
        parser.error(f"panel does not exist: {panel_path}")
    _set_panel_path(args, panel_path)


def _resolve_call_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if not math.isfinite(args.missing_locus_min_depth) or args.missing_locus_min_depth <= 0:
        parser.error("--missing-locus-min-depth must be finite and positive")
    if not math.isfinite(args.missing_locus_min_fraction) or not 0 < args.missing_locus_min_fraction <= 1:
        parser.error("--missing-locus-min-fraction must be in (0, 1]")
    if not math.isfinite(args.missing_locus_penalty) or args.missing_locus_penalty < 0:
        parser.error("--missing-locus-penalty must be finite and non-negative")
    if not math.isfinite(args.classification_repeat_scale) or args.classification_repeat_scale <= 0:
        parser.error("--classification-repeat-scale must be finite and positive")
    _resolve_panel_option(parser, args)
    database = Path(args.database).resolve() if args.database else None
    legacy_multi_taxon_build = bool(
        database
        and (database / "reference_pipeline_manifest.json").is_file()
        and not (database / "database" / "reference_panel.tsv").is_file()
    )
    if legacy_multi_taxon_build:
        parser.error(
            "This reference database predates the competitive FASTQ mapping schema. "
            "Rebuild it with the current `mlvamaps build-reference`."
        )
    if database:
        root = database / "database" if (database / "database").is_dir() else database
        manifest = database / "manifest.json"
        if not manifest.is_file():
            manifest = root.parent / "manifest.json"
        required = [
            root / "competitive_mapping" / "candidate_contexts.fasta",
            root / "competitive_mapping" / "candidate_metadata.tsv",
            root / "competitive_mapping" / "short.mmi",
            root / "competitive_mapping" / "long.mmi",
        ]
        if not manifest.is_file() or any(not path.is_file() for path in required):
            parser.error(
                "This reference database predates the competitive FASTQ mapping schema. "
                "Rebuild it with the current `mlvamaps build-reference`."
            )
    explicit_short_read_mode = args.input_path == "sr"
    args.short_read_mode = explicit_short_read_mode or args.short_reads
    if explicit_short_read_mode:
        args.input_path = None
    if args.short_reads and (args.reads1 or args.reads2 or args.manifest):
        parser.error("--short-reads directory mode cannot be combined with --fq1, --fq2, or --manifest")
    if args.short_reads and explicit_short_read_mode:
        parser.error("--short-reads requires -i DIRECTORY, not -i sr")
    if not args.short_read_mode and (args.reads1 or args.reads2):
        parser.error("--fq1/--fq2 require the short-read selector: -i sr")
    if args.taxon_identification is True and not args.database:
        parser.error("--taxon-identification requires --database")
    if not args.loci and not args.primers and args.database:
        database = Path(args.database)
        for candidate in (database / "database" / "reference_panel.tsv", database / "reference_panel.tsv"):
            if candidate.is_file():
                args.loci = str(candidate)
                break
    if not args.loci and not args.primers:
        parser.error("call requires -p PANEL unless --database contains reference_panel.tsv")
    if not args.input_path and not args.reads1 and not args.manifest:
        if args.short_read_mode:
            parser.error("-i sr requires --fq1 FASTQ or --manifest TSV")
        parser.error("call requires -i INPUT")


def build_parser(*, advanced: bool = False) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mlvamaps",
        description="Simple MLVA/VNTR calling from primers plus FASTQ, FASTA, or an input directory",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    call = subparsers.add_parser(
        "call",
        help="Call VNTRs from primers plus FASTQ/FASTA files or a directory",
        epilog=(
            "Use mlvamaps call --advanced for all options.\n\n"
            "Examples:\n"
            "  mlvamaps call -p primers.tsv -i sample.fastq.gz\n"
            "  mlvamaps call -p primers.tsv -i assembly.fasta\n"
            "  mlvamaps call -p primers.tsv -i sequence_files/\n"
            "  mlvamaps call -p primers.tsv -i sr --fq1 R1.fastq.gz --fq2 R2.fastq.gz"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    call.add_argument(
        "-i",
        "--input",
        dest="input_path",
        required=False,
        metavar="INPUT",
        help="FASTQ/FASTA path or directory; use 'sr' with --fq1/--fq2 for short reads",
    )
    call.add_argument("--reads", dest="reads_path", metavar="FASTQ", help="Reads to map for assembly depth support")
    call.add_argument(
        "--fq1",
        dest="reads1",
        metavar="FASTQ",
        help="Mate 1 or single-end Illumina FASTQ (explicit short-read interface)",
    )
    call.add_argument(
        "--fq2",
        dest="reads2",
        metavar="FASTQ",
        help="Mate 2 Illumina FASTQ; records must be in the same order as --fq1",
    )
    call.add_argument(
        "--read-technology",
        choices=("auto", "illumina", "accurate-long", "ont", "ont-hq", "hifi"),
        default="auto",
        help="Read evidence model compatibility override (default: %(default)s; -i sr selects Illumina)",
    )
    call.add_argument(
        "--manifest",
        metavar="TSV",
        help="Batch TSV with sample_id, reads1, and optional reads2/metadata_id columns",
    )
    call.add_argument(
        "--short-reads",
        action="store_true",
        help=(
            "Treat -i as a directory and pair exact PREFIX_1.fastq.gz and "
            "PREFIX_2.fastq.gz filenames"
        ),
    )
    call.add_argument(
        "--sample-metadata",
        metavar="CSV_OR_TSV",
        help="One-row-per-sample metadata joined by sample_id and exported for MYOGA",
    )
    call.add_argument("--force", action="store_true", help="Reprocess successful manifest samples")
    call.add_argument("--keep-intermediates", action="store_true", help="Retain FASTQ candidate-mapping alignments")
    call.add_argument("--bam", "--alignments", dest="alignments_path", metavar="BAM/SAM", help="Assembly-aligned BAM/SAM for assembly depth support")
    call.add_argument(
        "-p",
        "--panel",
        dest="panel_path",
        required=False,
        help="Primer/locus panel (optional when --database is a reference build)",
    )
    call.set_defaults(loci=None, primers=None)
    call.add_argument("--profiles")
    call.add_argument(
        "--database",
        help="Reference-build directory for alignment-based reference classification",
    )
    taxon_toggle = call.add_mutually_exclusive_group()
    taxon_toggle.add_argument(
        "--taxon-identification",
        dest="taxon_identification",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    taxon_toggle.add_argument(
        "--no-taxon-identification",
        dest="taxon_identification",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    call.set_defaults(taxon_identification=None)
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
    call.add_argument("--short-min-read-length", type=_positive_int, default=40)
    call.add_argument("--short-min-mean-quality", type=_nonnegative_float, default=15.0)
    call.add_argument("--short-trim-quality", type=_nonnegative_int, default=0)
    call.add_argument("--short-min-pair-retention", type=_fraction, default=0.5)
    call.add_argument("--short-min-mapq", type=_nonnegative_int, default=0)
    call.add_argument("--short-min-spanning-pairs", type=_nonnegative_int, default=2)
    call.add_argument("--short-confidence-threshold", type=_fraction, default=0.8)
    call.add_argument("--short-max-candidate-repeat-count", type=_positive_int, default=100)
    call.add_argument(
        "--no-short-secondary-alignments", action="store_true",
        help="Ignore secondary candidate alignments (not recommended for homologous contexts)",
    )
    call.add_argument(
        "--short-min-informative-molecules",
        type=_positive_int,
        default=3,
        help="Informative molecules required to avoid Illumina LOW_DEPTH (default: %(default)s)",
    )
    call.add_argument("--min-read-length", type=int, default=50)
    call.add_argument("--max-read-length", type=int, default=100000)
    call.add_argument(
        "--min-qscore",
        type=float,
        default=15.0,
        help=(
            "Minimum mean read Phred score (default: %(default)s, approximately "
            "97%% per-base accuracy)"
        ),
    )
    call.add_argument(
        "--recruitment-database",
        help=(
            "Reference-build directory supplying canonical locus products for "
            "FASTQ recruitment without reference classification"
        ),
    )
    call.add_argument("--max-primer-mismatches", type=int, default=2)
    call.add_argument(
        "--assembly-round-tolerance",
        type=_round_tolerance,
        default=0.25,
        metavar="FRACTION",
        help="Legacy integer-rounding tolerance (default: %(default)s)",
    )
    call.add_argument(
        "--sample-mode",
        choices=("isolate", "metagenome"),
        default="metagenome",
        help="Interpret FASTQ evidence as an isolate or metagenome (default: %(default)s)",
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
        "--minimap2-bin",
        default="minimap2",
        metavar="PATH",
        help="minimap2 executable for short-read recruitment, representative mapping, and assembly support (default: %(default)s)",
    )
    call.add_argument(
        "--classification-repeat-scale", type=_positive_float, default=1.0,
        help="Repeat-count difference scale in mapping likelihoods, in repeat units (default: %(default)s)",
    )
    call.add_argument(
        "--missing-locus-min-depth", type=_positive_float, default=3.0,
        help="Minimum supporting molecules per locus for the missing-locus gate (default: %(default)s)",
    )
    call.add_argument(
        "--missing-locus-min-fraction", type=_positive_float, default=0.8,
        help="Fraction of panel loci meeting missing-locus depth, in (0, 1] (default: %(default)s)",
    )
    call.add_argument(
        "--missing-locus-penalty", type=_nonnegative_float, default=1.0,
        help="Log-likelihood penalty per undetected query locus present in a reference after the coverage gate passes; 0 disables (default: %(default)s)",
    )
    call.add_argument(
        "--taxon-min-loci",
        type=_positive_int,
        default=None,
        help="Minimum observed loci for model-supported mapping classification (default: 2)",
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

    call.add_argument("--advanced", action="help", help="Show all options and exit")
    if not advanced:
        basic_options = {
            "help", "advanced", "input_path", "reads1", "reads2", "short_reads",
            "manifest", "panel_path", "database", "outdir", "sample_id",
            "sample_mode", "threads", "quiet",
        }
        for action in call._actions:
            if action.dest not in basic_options:
                action.help = argparse.SUPPRESS

    export = subparsers.add_parser(
        "export-myoga",
        help="Aggregate completed MLVA results into a MYOGA-ready relatedness dataset",
        description=(
            "Discover completed mlvamaps samples, filter exact repeat-count profiles, "
            "join metadata, calculate shared-locus distances, and write a deterministic "
            "neighbor-joining MLVA profile-similarity tree without rerunning calling."
        ),
    )
    export.add_argument(
        "--results",
        required=True,
        metavar="DIR",
        help="Root containing completed per-sample mlvamaps result directories",
    )
    export.add_argument(
        "--metadata",
        required=True,
        metavar="CSV_OR_TSV",
        help="Sample metadata table to join to recorded mlvamaps sample IDs",
    )
    export.add_argument(
        "--metadata-id",
        default="shared_identifier",
        metavar="COLUMN",
        help="Metadata identifier column matched to sample_id (default: %(default)s)",
    )
    export.add_argument(
        "--latitude",
        default="latitude",
        metavar="COLUMN",
        help="Latitude column (default: %(default)s; standard aliases are recognized)",
    )
    export.add_argument(
        "--longitude",
        default="longitude",
        metavar="COLUMN",
        help="Longitude column (default: %(default)s; standard aliases are recognized)",
    )
    export.add_argument(
        "--min-callable-fraction",
        type=_fraction,
        default=0.0,
        metavar="FRACTION",
        help=(
            "Minimum exact VNTR-call fraction per sample; by default, retain any "
            "sample with at least one finite repeat_count (default: %(default)s)"
        ),
    )
    export.add_argument(
        "--min-callable-loci",
        type=_nonnegative_int,
        default=0,
        metavar="COUNT",
        help="Additional minimum exact-call count; both sample thresholds must pass (default: %(default)s)",
    )
    export.add_argument(
        "--min-pairwise-loci",
        type=_positive_int,
        default=1,
        metavar="COUNT",
        help="Minimum shared exact loci per supported pair (default: %(default)s)",
    )
    export.add_argument(
        "--min-pairwise-fraction",
        type=_fraction,
        default=0.0,
        metavar="FRACTION",
        help="Minimum panel fraction callable in both samples (default: %(default)s)",
    )
    export.add_argument(
        "--distance",
        choices=("repeat", "categorical"),
        default="repeat",
        help="Distance used for the matrix and tree (default: %(default)s)",
    )
    export.add_argument(
        "-o",
        "--output",
        "--outdir",
        dest="outdir",
        required=True,
        metavar="DIR",
        help="Output directory for the export",
    )
    export.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing export in the output directory",
    )

    reference = subparsers.add_parser(
        "build-reference",
        help="Build a complete reference database from local assemblies or NCBI taxids",
        description=(
            "Build the complete mlvamaps database: fetch references when needed, "
            "extract real loci, summarize amplifiability, "
            "generate competitive allele contexts, build short/long minimap2 "
            "indexes, and build broad Deacon recruitment assets."
        ),
    )
    reference_source = reference.add_mutually_exclusive_group(required=True)
    reference_source.add_argument(
        "-i",
        "--input",
        dest="assemblies",
        help="Directory containing reference FASTA assemblies",
    )
    reference_source.add_argument("--taxid", help="Single NCBI taxonomy identifier")
    reference_source.add_argument(
        "--taxids",
        "--taxids-csv",
        dest="taxids_csv",
        help="CSV/TSV with a taxid column and optional name column",
    )
    panel = reference.add_mutually_exclusive_group(required=True)
    panel.add_argument(
        "-p",
        "--panel",
        "--primers",
        dest="panel_path",
        help="Primer list or rich locus panel",
    )
    reference.set_defaults(loci=None, primers=None)
    reference.add_argument(
        "--metadata",
        help="Reference metadata CSV/TSV; required with local assembly input via -i",
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
    reference.add_argument(
        "--download-retries",
        type=_positive_int,
        default=3,
        help="NCBI download attempts after transient failures (default: %(default)s)",
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
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="Overall build parallelism (default: %(default)s; 0 auto-detects CPUs)",
    )
    reference.add_argument("--quiet", action="store_true", help="Suppress live progress updates")
    reference.add_argument("--minimap2-bin", default="minimap2", help=argparse.SUPPRESS)
    reference.add_argument("--deacon-bin", default="deacon", help=argparse.SUPPRESS)

    return parser


def _run_single_input(
    args: argparse.Namespace,
    input_path: Path,
    outdir: Path,
    sample_id: str,
) -> dict[str, Path]:
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
            min_mixture_fraction=args.min_mixture_fraction,
            min_secondary_reads=args.min_secondary_reads,
            minimap2_bin=args.minimap2_bin,
            missing_locus_min_depth=args.missing_locus_min_depth,
            missing_locus_min_fraction=args.missing_locus_min_fraction,
            missing_locus_penalty=args.missing_locus_penalty,
            classification_repeat_scale=args.classification_repeat_scale,
            reference_metadata_path=args.reference_metadata,
            taxon_min_loci=args.taxon_min_loci,
            taxon_identification=args.taxon_identification,
            locus_mapping=not args.no_locus_mapping,
            min_mapping_quality=args.min_mapping_quality,
            min_base_quality=args.min_base_quality,
            min_snp_depth=args.min_snp_depth,
            min_snp_alternate_reads=args.min_snp_alternate_reads,
            min_snp_frequency=args.min_snp_frequency,
            threads=args.threads,
            show_progress=not args.quiet,
            sample_mode=args.sample_mode,
            assembly_round_tolerance=args.assembly_round_tolerance,
            max_confidence_depth=args.max_confidence_depth,
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
            print(f"Wrote reference matches to {result['mapping_reference_matches']}")
            if "mlva_profile_tree" in result:
                print(f"Wrote MLVA profile tree to {result['mlva_profile_tree']}")
            if "taxonomic_identification" in result:
                print(f"Wrote automatic taxonomic identification to {result['taxonomic_identification']}")
        return result

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
        threads=args.threads,
        minimap2_preset=args.minimap2_preset,
        minimap2_bin=args.minimap2_bin,
        missing_locus_min_depth=args.missing_locus_min_depth,
        missing_locus_min_fraction=args.missing_locus_min_fraction,
        missing_locus_penalty=args.missing_locus_penalty,
        classification_repeat_scale=args.classification_repeat_scale,
        reference_metadata_path=args.reference_metadata,
        taxon_min_loci=args.taxon_min_loci,
        taxon_identification=args.taxon_identification,
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
        print(f"Wrote reference matches to {result['mapping_reference_matches']}")
        if "mlva_profile_tree" in result:
            print(f"Wrote MLVA profile tree to {result['mlva_profile_tree']}")
        if "taxonomic_identification" in result:
            print(f"Wrote automatic taxonomic identification to {result['taxonomic_identification']}")
    if args.reads_path or args.alignments_path:
        print(f"Wrote read-depth support to {result['read_support']}")
    return result


def _run_short_input(
    args: argparse.Namespace,
    reads1: Path,
    reads2: Path | None,
    outdir: Path,
    sample_id: str,
    metadata: dict[str, str] | None = None,
) -> dict[str, Path]:
    result = run_short_read_call(
        reads1_path=str(reads1),
        reads2_path=None if reads2 is None else str(reads2),
        loci_path=args.loci,
        primers_path=args.primers,
        profiles_path=args.profiles,
        database_path=args.database or args.recruitment_database,
        outdir=str(outdir),
        sample_id=sample_id,
        sample_metadata=metadata,
        short_min_read_length=args.short_min_read_length,
        short_min_mean_quality=args.short_min_mean_quality,
        short_trim_quality=args.short_trim_quality,
        short_min_pair_retention=args.short_min_pair_retention,
        min_depth=args.short_min_informative_molecules,
        threads=args.threads,
        keep_intermediates=args.keep_intermediates,
        sample_mode=args.sample_mode,
        minimap2_bin=args.minimap2_bin,
        short_min_mapping_quality=args.short_min_mapq,
        short_min_spanning_pairs=args.short_min_spanning_pairs,
        short_confidence_threshold=args.short_confidence_threshold,
        short_max_candidate_repeat_count=args.short_max_candidate_repeat_count,
        short_consider_secondary=not args.no_short_secondary_alignments,
        missing_locus_min_depth=args.missing_locus_min_depth,
        missing_locus_min_fraction=args.missing_locus_min_fraction,
        missing_locus_penalty=args.missing_locus_penalty,
        classification_repeat_scale=args.classification_repeat_scale,
        reference_metadata_path=args.reference_metadata,
        taxon_min_loci=args.taxon_min_loci,
        taxon_identification=args.taxon_identification,
        show_progress=not args.quiet,
    )
    print(f"Wrote conservative Illumina calls to {result['calls']}")
    print(f"Wrote short-read QC to {result['short_read_qc']}")
    print(f"Wrote locus recruitment to {result['short_read_recruitment']}")
    print(f"Wrote minimap2-derived mapping evidence to {result['short_read_mapping']}")
    print(f"Wrote MYOGA metadata to {result['myoga_samples']}")
    if "mapping_reference_matches" in result:
        print(f"Wrote mapping reference matches to {result['mapping_reference_matches']}")
    if "mlva_profile_tree" in result:
        print(f"Wrote MLVA profile tree to {result['mlva_profile_tree']}")
    if "taxonomic_identification" in result:
        print(f"Wrote automatic taxonomic identification to {result['taxonomic_identification']}")
    print(f"Wrote report to {result['report']}")
    return result


def _read_manifest(path: str | Path) -> list[dict[str, str]]:
    manifest_path = Path(path)
    with open_text(manifest_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = {str(field).strip().lower() for field in (reader.fieldnames or [])}
        if not {"sample_id", "reads1"}.issubset(fields):
            raise ValueError("manifest requires sample_id and reads1 columns")
        rows = [
            {
                str(key).strip().lower(): "" if value is None else str(value).strip()
                for key, value in row.items()
                if key is not None
            }
            for row in reader
        ]
    sample_ids = [row.get("sample_id", "") for row in rows]
    if any(not sample_id for sample_id in sample_ids):
        raise ValueError("manifest sample_id values cannot be empty")
    duplicates = sorted(
        sample_id for sample_id in set(sample_ids) if sample_ids.count(sample_id) > 1
    )
    if duplicates:
        raise ValueError("manifest sample_id values must be unique: " + ", ".join(duplicates))
    base = manifest_path.parent
    pending: list[tuple[int, dict[str, str], Path, Path | None, Path, object]] = []
    for row_index, row in enumerate(rows):
        for field in ("reads1", "reads2"):
            value = row.get(field, "")
            if value in ("", "."):
                row[field] = ""
            else:
                candidate = Path(value)
                row[field] = str(candidate if candidate.is_absolute() else base / candidate)
    return rows


def _read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _combine_tables(paths: list[Path], output: Path, delimiter: str = "\t") -> Path:
    rows: list[dict[str, str]] = []
    fields: list[str] = []
    for path in paths:
        if not path.exists():
            continue
        for row in _read_table(path, delimiter):
            rows.append(row)
            for field in row:
                if field not in fields:
                    fields.append(field)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    return output


def _metadata_for_manifest_row(
    row: dict[str, str],
    by_sample: dict[str, dict[str, str]],
    metadata_rows: list[dict[str, str]],
) -> dict[str, str] | None:
    if row["sample_id"] in by_sample:
        return by_sample[row["sample_id"]]
    metadata_id = row.get("metadata_id", "")
    if not metadata_id:
        return None
    return next(
        (
            metadata
            for metadata in metadata_rows
            if metadata_id in {metadata.get("biosample", ""), metadata.get("run_accession", "")}
        ),
        None,
    )


def _run_short_batch(
    args: argparse.Namespace,
    parser: argparse.ArgumentParser,
    rows: list[dict[str, str]],
) -> None:
    try:
        metadata_rows = read_sample_metadata(args.sample_metadata)
    except ValueError as exc:
        parser.error(str(exc))
    by_sample = metadata_by_sample(metadata_rows)
    results: list[dict[str, Path]] = []
    statuses: list[dict[str, str]] = []
    output_root = Path(args.outdir)
    output_root.mkdir(parents=True, exist_ok=True)
    batch_outdir = _batch_output_dir(output_root)
    try:
        _check_batch_sample_ids([row["sample_id"] for row in rows])
    except ValueError as exc:
        parser.error(str(exc))
    batch_outdir.mkdir(parents=True, exist_ok=True)
    pending: list[tuple[int, dict[str, str], Path, Path | None, Path, object]] = []
    for row_index, row in enumerate(rows):
        sample_id = row["sample_id"]
        sample_outdir = _sample_output_dir(output_root, sample_id)
        summary_path = sample_outdir / "sample_summary.tsv"
        classification_ready = not (args.database or args.recruitment_database) or (
            sample_outdir / "classification" / "mapping_reference_matches.tsv"
        ).is_file() and (
            sample_outdir / "classification" / "classification.json"
        ).is_file() and summary_path.is_file() and summary_path.stat().st_mtime_ns >= (
            sample_outdir / "classification" / "classification.json"
        ).stat().st_mtime_ns
        if summary_path.exists() and not args.force and classification_ready:
            summary_rows = _read_table(summary_path)
            if summary_rows and summary_rows[0].get("run_status") == "success":
                statuses.append({"sample_id": sample_id, "status": "skipped_success", "message": "existing successful result; use --force to rerun"})
                results.append({
                    "calls": sample_outdir / "calls.tsv",
                    "repeat_counts": sample_outdir / "locus_repeat_counts.tsv",
                    "fingerprint": sample_outdir / "mlva_fingerprint.tsv",
                    "profile_matches": sample_outdir / "profile_matches.tsv",
                    "profile_match_loci": sample_outdir / "profile_match_loci.tsv",
                    "sample_summary": summary_path,
                    "myoga_samples": sample_outdir / "myoga_samples.csv",
                    "myoga_loci": sample_outdir / "myoga_loci.csv",
                    **{
                        key: sample_outdir / "classification" / filename
                        for key, filename in (
                            ("taxonomic_identification", "taxonomic_identification.tsv"),
                            ("taxonomic_identification_evidence", "taxonomic_identification_evidence.tsv"),
                            ("mapping_reference_matches", "mapping_reference_matches.tsv"),
                        )
                        if (sample_outdir / "classification" / filename).is_file()
                    },
                })
                continue
        reads1 = Path(row["reads1"])
        reads2 = Path(row["reads2"]) if row.get("reads2") else None
        pending.append((row_index, row, reads1, reads2, sample_outdir, _metadata_for_manifest_row(row, by_sample, metadata_rows)))
    workers, sample_threads = _batch_allocation(args.threads, len(pending))
    print(f"Processing {len(pending)} Illumina sample(s) with {workers} worker(s), {sample_threads} thread(s) each")
    completed = {}
    def run_one(item):
        index, row, reads1, reads2, sample_outdir, metadata = item
        if not reads1.is_file():
            raise ValueError(f"reads1 does not exist: {reads1}")
        if reads2 is not None and not reads2.is_file():
            raise ValueError(f"reads2 does not exist: {reads2}")
        return index, _run_short_input(_args_with_threads(args, sample_threads), reads1, reads2, sample_outdir, row["sample_id"], metadata)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_one, item): item for item in pending}
        for future in as_completed(futures):
            item = futures[future]
            index, row = item[0], item[1]
            try:
                _, result = future.result()
                completed[index] = (result, {"sample_id": row["sample_id"], "status": "success", "message": ""})
            except Exception as exc:
                completed[index] = (None, {"sample_id": row["sample_id"], "status": "failed", "message": f"{type(exc).__name__}: {exc}"})
                print(f"Sample {row['sample_id']} failed: {type(exc).__name__}: {exc}")
    for index in sorted(completed):
        result, status = completed[index]
        statuses.append(status)
        if result is not None:
            results.append(result)
    write_tsv(statuses, batch_outdir / "batch_status.tsv", ["sample_id", "status", "message"])
    table_keys = {
        "calls": "calls.tsv",
        "repeat_counts": "locus_repeat_counts.tsv",
        "fingerprint": "mlva_fingerprint.tsv",
        "profile_matches": "profile_matches.tsv",
        "profile_match_loci": "profile_match_loci.tsv",
        "sample_summary": "sample_summary.tsv",
        "taxonomic_identification": "taxonomic_identification.tsv",
        "taxonomic_identification_evidence": "taxonomic_identification_evidence.tsv",
        "mapping_reference_matches": "mapping_reference_matches.tsv",
    }
    for key, filename in table_keys.items():
        _combine_tables(
            [result[key] for result in results if key in result],
            batch_outdir / filename,
        )
    myoga_sample_rows = [row for result in results if "myoga_samples" in result for row in _read_table(result["myoga_samples"], ",")]
    write_csv(myoga_sample_rows, batch_outdir / "myoga_samples.csv", MYOGA_SAMPLE_FIELDS)
    _combine_tables(
        [result["myoga_loci"] for result in results if "myoga_loci" in result],
        batch_outdir / "myoga_loci.csv",
        ",",
    )


def _run_manifest(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    try:
        manifest_rows = _read_manifest(args.manifest)
    except ValueError as exc:
        parser.error(str(exc))
    _run_short_batch(args, parser, manifest_rows)


def _combine_legacy_fingerprints(
    fingerprint_paths: list[Path], output_path: Path
) -> Path:
    """Combine single-sample MLVA_finder tables into one compatible CSV.

    Each assembly run writes only one data row, so a streaming CSV pass is
    faster and substantially lighter than constructing a dataframe for this
    operation. The panel header must be identical across all runs.
    """
    header: list[str] | None = None
    if not fingerprint_paths:
        raise ValueError("No MLVA fingerprints were provided for aggregation")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        index = 0
        for fingerprint_path in fingerprint_paths:
            with fingerprint_path.open(newline="", encoding="utf-8") as source:
                reader = csv.reader(source)
                current_header = next(reader, None)
                if current_header is None:
                    raise ValueError(f"Empty MLVA fingerprint: {fingerprint_path}")
                if header is None:
                    header = current_header
                    writer.writerow(header)
                elif current_header != header:
                    raise ValueError(
                        "Cannot combine MLVA fingerprints with different locus "
                        f"columns: {fingerprint_path}"
                    )
                for row in reader:
                    index += 1
                    writer.writerow([f"{index:03d}", *row[1:]])
    return output_path


def _batch_analysis_path(input_path: str, outdir: str) -> Path:
    directory_name = Path(input_path).resolve().name
    return _batch_output_dir(Path(outdir)) / f"MLVA_analysis_{directory_name}.csv"


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser(advanced="--advanced" in argv)
    args = parser.parse_args(argv)
    if args.command == "export-myoga":
        try:
            result = export_myoga(
                args.results,
                args.metadata,
                args.outdir,
                metadata_id=args.metadata_id,
                latitude=args.latitude,
                longitude=args.longitude,
                min_callable_fraction=args.min_callable_fraction,
                min_callable_loci=args.min_callable_loci,
                min_pairwise_loci=args.min_pairwise_loci,
                min_pairwise_fraction=args.min_pairwise_fraction,
                distance=args.distance,
                force=args.force,
            )
        except ValueError as exc:
            parser.error(str(exc))
        print(f"Wrote MYOGA metadata to {result['metadata']}")
        print(f"Wrote MLVA distance matrix to {result['distance_matrix']}")
        if result["tree"]:
            print(f"Wrote MLVA relatedness tree to {result['tree']}")
        else:
            print("No MLVA relatedness tree was written because no samples passed filtering")
        print(f"Wrote export summary to {result['summary']}")
        return 0
    if args.command == "call":
        _resolve_call_args(parser, args)
        if args.reads2 and not args.reads1:
            parser.error("--fq2 requires --fq1")
        if args.manifest and not args.short_read_mode:
            parser.error("--manifest requires the short-read selector: -i sr")
        explicit_sources = sum(bool(value) for value in (args.input_path, args.reads1, args.manifest))
        if explicit_sources != 1:
            parser.error("choose exactly one sample source with -i INPUT or -i sr")
        if args.manifest and args.sample_id:
            parser.error("--sample-id cannot be used with --manifest; sample IDs come from the manifest")
        if args.manifest and args.read_technology == "accurate-long":
            parser.error("--manifest currently implements the Illumina evidence model only")
        if args.short_read_mode and args.read_technology == "accurate-long":
            parser.error("-i sr selects Illumina and cannot use --read-technology accurate-long")
        if not args.short_read_mode and args.read_technology == "illumina":
            parser.error("Illumina input requires -i sr with --fq1 and optional --fq2")
        if args.manifest:
            _run_manifest(args, parser)
            return 0
        if args.short_reads:
            if args.sample_id:
                parser.error(
                    "--sample-id cannot be used with --short-reads; sample IDs come from filename prefixes"
                )
            if args.reads_path or args.alignments_path:
                parser.error(
                    "--reads and --bam/--alignments cannot be combined with --short-reads"
                )
            try:
                directory_rows = _short_read_directory_rows(args.input_path)
            except ValueError as exc:
                parser.error(str(exc))
            _run_short_batch(args, parser, directory_rows)
            return 0
        try:
            metadata_rows = read_sample_metadata(args.sample_metadata)
        except ValueError as exc:
            parser.error(str(exc))
        sample_metadata = metadata_by_sample(metadata_rows)
        if args.reads_path and args.alignments_path:
            parser.error("call accepts either --reads or --bam/--alignments for assembly depth support, not both")
        if args.reads1:
            reads1 = Path(args.reads1)
            reads2 = Path(args.reads2) if args.reads2 else None
            if not reads1.is_file():
                parser.error(f"--fq1 does not exist: {reads1}")
            if reads2 is not None and not reads2.is_file():
                parser.error(f"--fq2 does not exist: {reads2}")
            technology = "illumina" if args.short_read_mode else args.read_technology
            if reads2 is not None and technology != "illumina":
                parser.error("paired --fq1/--fq2 input requires Illumina mode")
            sample_id = args.sample_id or _sample_id_from_path(str(reads1))
            metadata = sample_metadata.get(sample_id)
            if technology == "illumina":
                _run_short_input(args, reads1, reads2, Path(args.outdir), sample_id, metadata)
            else:
                _run_single_input(args, reads1, Path(args.outdir), sample_id)
            return 0
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
        if batch:
            try:
                _check_batch_sample_ids(sample_ids)
            except ValueError as exc:
                parser.error(str(exc))
            _batch_output_dir(Path(args.outdir)).mkdir(parents=True, exist_ok=True)
        legacy_fingerprints = []
        for input_path, derived_sample_id in zip(input_files, sample_ids):
            sample_id = args.sample_id or derived_sample_id
            outdir = (
                _sample_output_dir(Path(args.outdir), sample_id)
                if batch
                else Path(args.outdir)
            )
            if batch:
                print(f"Processing {input_path} as sample {sample_id}")
            technology = args.read_technology
            if technology == "illumina":
                if _input_kind(str(input_path)) != "fastq":
                    parser.error("--read-technology illumina requires FASTQ input; assembly plus Illumina support must use a separate call")
                result = _run_short_input(
                    args,
                    input_path,
                    None,
                    outdir,
                    sample_id,
                    sample_metadata.get(sample_id),
                )
            else:
                result = _run_single_input(args, input_path, outdir, sample_id)
            if result and "legacy_fingerprint" in result:
                legacy_fingerprints.append(Path(result["legacy_fingerprint"]))
        if batch and legacy_fingerprints:
            analysis_path = _combine_legacy_fingerprints(
                legacy_fingerprints,
                _batch_analysis_path(args.input_path, args.outdir),
            )
            print(f"Wrote combined MLVA_finder analysis to {analysis_path}")
        return 0
    if args.command == "build-reference":
        _resolve_panel_option(parser, args)
        if args.assemblies and not args.metadata:
            parser.error("build-reference with -i ASSEMBLIES also requires --metadata")
        if not args.assemblies and args.metadata:
            parser.error("--metadata is only accepted with -i ASSEMBLIES")
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
                download_retries=args.download_retries,
                multiple_products=args.multiple_products,
                max_primer_mismatches=args.max_primer_mismatches,
                threads=args.threads,
                minimap2_bin=args.minimap2_bin,
                deacon_bin=args.deacon_bin,
                show_progress=not args.quiet,
            )
            if not args.quiet:
                print(f"Wrote taxon summary to {result['taxon_summary']}")
                print(f"Wrote taxon/locus summary to {result['locus_amplifiability']}")
                print(f"Wrote combined taxon-identification database to {result['database']}")
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
            threads=args.threads,
            minimap2_bin=args.minimap2_bin,
            deacon_bin=args.deacon_bin,
            show_progress=not args.quiet,
        )
        print(f"Wrote per-locus reference database to {result['database']}")
        print(f"Wrote reference build QC to {result['manifest']}")
        print(f"Reference build completed with status {result['status']}")
        print(f"Wrote MYOGA metadata to {result['myoga_metadata']}")
        return 0
    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
