from __future__ import annotations

import argparse
from pathlib import Path

from .assembly_call import run_assembly_call
from .concurrency import DEFAULT_THREADS
from .in_silico_pcr import run_amplirust
from .io import open_text
from .pipeline import run_call
from .simulation import simulate_reads


def _sample_id_from_path(path: str) -> str:
    sample = Path(path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fna.gz"):
        if sample.endswith(suffix):
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


def _input_kind(path: str) -> str:
    lower = path.lower()
    if lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"
    if lower.endswith((".fasta", ".fa", ".fna", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz")):
        return "fasta"
    with open(path) as handle:
        first = handle.read(1)
    if first == "@":
        return "fastq"
    if first == ">":
        return "fasta"
    raise ValueError(f"Could not tell whether {path!r} is FASTQ reads or FASTA assembly.")


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
        parser.error("call accepts at most two positional paths: primers.tsv and sample.fastq.gz or assembly.fasta")
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
        parser.error("call requires a primer file, for example: mlvamaps call primers.tsv sample.fastq.gz")
    if not args.input_path:
        parser.error("call requires an input FASTQ or FASTA file")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mlvamaps",
        description="Simple MLVA/VNTR calling from primers plus FASTQ or FASTA",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    call = subparsers.add_parser(
        "call",
        help="Call VNTRs from primers plus FASTQ reads or a FASTA assembly",
        epilog=(
            "Examples:\n"
            "  mlvamaps call primers.tsv sample.fastq.gz\n"
            "  mlvamaps call primers.tsv assembly.fasta\n"
            "  mlvamaps call primers.tsv assembly.fasta --reads sample.fastq.gz\n"
            "  mlvamaps call primers.tsv assembly.fasta --bam assembly_reads.bam"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    call.add_argument("paths", nargs="*", metavar="PATH", help="primers.tsv plus sample.fastq.gz or assembly.fasta")
    call.add_argument("--input", dest="input_path", metavar="PATH", help="FASTQ reads or FASTA assembly")
    call.add_argument("--reads", dest="reads_path", metavar="FASTQ", help="Reads to map for assembly depth support")
    call.add_argument("--bam", "--alignments", dest="alignments_path", metavar="BAM/SAM", help="Assembly-aligned BAM/SAM for assembly depth support")
    call.add_argument("--loci")
    call.add_argument("--primers", help="Primer-pair CSV/TSV/whitespace file with locus, forward, reverse columns")
    call.add_argument("--profiles")
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
    call.add_argument("--min-qscore", type=float, default=0.0)
    call.add_argument("--max-primer-mismatches", type=int, default=3)
    call.add_argument("--min-depth", type=int, default=10)
    call.add_argument("--min-posterior", type=float, default=0.75)
    call.add_argument(
        "--min-cluster-size",
        type=_positive_int,
        default=2,
        help="Minimum read support for a retained VNTR cluster (default: %(default)s)",
    )
    call.add_argument(
        "--cluster-min-identity",
        type=_fraction,
        default=0.97,
        help="Minimum VSEARCH global identity within a locus (default: %(default)s)",
    )
    call.add_argument(
        "--vsearch-bin",
        default="vsearch",
        metavar="PATH",
        help="VSEARCH executable (default: %(default)s)",
    )
    call.add_argument(
        "--amplirust-bin",
        default="amplirust",
        metavar="PATH",
        help="Amplirust executable used for degenerate primer pairing (default: %(default)s)",
    )
    call.add_argument(
        "--minimap2-bin",
        default="minimap2",
        metavar="PATH",
        help="minimap2 executable for representative and assembly-support mapping (default: %(default)s)",
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
        help="Extract expected MLVA amplicons from FASTA/GenBank with amplirust",
    )
    extract.add_argument("--input", required=True, help="FASTA/GenBank input or amplirust-supported glob")
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
        help="amplirust threads (default: %(default)s; 0 lets amplirust auto-detect CPUs)",
    )
    extract.add_argument("--circular", action="store_true")
    extract.add_argument("--no-search-rc", action="store_true")
    extract.add_argument("--trim-primers", action="store_true")
    extract.add_argument("--amplirust-bin", default="amplirust")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "call":
        _resolve_call_args(parser, args)
        sample_id = args.sample_id or _sample_id_from_path(args.input_path)
        if args.reads_path and args.alignments_path:
            parser.error("call accepts either --reads or --bam/--alignments for assembly depth support, not both")
        if _input_kind(args.input_path) == "fastq":
            result = run_call(
                reads_path=args.input_path,
                loci_path=args.loci,
                primers_path=args.primers,
                profiles_path=args.profiles,
                outdir=args.outdir,
                sample_id=sample_id,
                min_read_length=args.min_read_length,
                max_read_length=args.max_read_length,
                min_qscore=args.min_qscore,
                max_primer_mismatches=args.max_primer_mismatches,
                min_depth=args.min_depth,
                min_posterior=args.min_posterior,
                min_cluster_size=args.min_cluster_size,
                cluster_min_identity=args.cluster_min_identity,
                vsearch_bin=args.vsearch_bin,
                amplirust_bin=args.amplirust_bin,
                minimap2_bin=args.minimap2_bin,
                locus_mapping=not args.no_locus_mapping,
                min_mapping_quality=args.min_mapping_quality,
                min_base_quality=args.min_base_quality,
                min_snp_depth=args.min_snp_depth,
                min_snp_alternate_reads=args.min_snp_alternate_reads,
                min_snp_frequency=args.min_snp_frequency,
                threads=args.threads,
                show_progress=not args.quiet,
            )
            print(f"Wrote easy MLVA calls to {result['calls']}")
            print(f"Wrote detailed allele evidence to {result['allele_calls']}")
            print(f"Wrote VNTR variant clusters to {result['asv_table']}")
            print(f"Wrote per-read cluster and indel evidence to {result['asv_memberships']}")
            if not args.no_locus_mapping:
                print(f"Wrote locus mapping summaries to {result['mapping_summary']}")
                print(f"Wrote locus SNP evidence to {result['mapping_snps']}")
            print(f"Wrote report to {result['report']}")
        else:
            result = run_assembly_call(
                assembly_path=args.input_path,
                loci_path=args.loci,
                primers_path=args.primers,
                outdir=args.outdir,
                sample_id=sample_id,
                reads_path=args.reads_path,
                alignments_path=args.alignments_path,
                profiles_path=args.profiles,
                max_primer_mismatches=args.max_primer_mismatches,
                threads=args.threads,
                minimap2_preset=args.minimap2_preset,
                minimap2_bin=args.minimap2_bin,
                amplirust_bin=args.amplirust_bin,
                show_progress=not args.quiet,
            )
            print(f"Wrote easy MLVA calls to {result['calls']}")
            print(f"Wrote assembly amplicons to {result['amplicons']}")
            print(f"Wrote report to {result['report']}")
            if args.reads_path or args.alignments_path:
                print(f"Wrote read-depth support to {result['read_support']}")
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
        result = run_amplirust(
            input_path=args.input,
            loci_path=args.loci,
            primers_path=args.primers,
            outdir=args.outdir,
            max_errors=args.max_errors,
            threads=args.threads,
            circular=args.circular,
            search_rc=not args.no_search_rc,
            trim_primers=args.trim_primers,
            executable=args.amplirust_bin,
        )
        print(f"Wrote amplirust primer CSV to {result['primers']}")
        print(f"Wrote extracted amplicons to {result['products']}")
        print(f"Wrote amplirust stats to {result['stats']}")
        return 0
    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
