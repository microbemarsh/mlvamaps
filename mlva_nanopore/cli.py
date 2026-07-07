from __future__ import annotations

import argparse
from pathlib import Path

from .assembly_call import run_assembly_call
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
        parser.error("call requires a primer file, for example: mlva-seer call primers.tsv sample.fastq.gz")
    if not args.input_path:
        parser.error("call requires an input FASTQ or FASTA file")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mlva-seer",
        description="Simple MLVA/VNTR calling from primers plus FASTQ or FASTA",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    call = subparsers.add_parser(
        "call",
        help="Call VNTRs from primers plus FASTQ reads or a FASTA assembly",
        epilog=(
            "Examples:\n"
            "  mlva-seer call primers.tsv sample.fastq.gz\n"
            "  mlva-seer call primers.tsv assembly.fasta\n"
            "  mlva-seer call primers.tsv assembly.fasta --reads sample.fastq.gz"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    call.add_argument("paths", nargs="*", metavar="PATH", help="primers.tsv plus sample.fastq.gz or assembly.fasta")
    call.add_argument("--input", dest="input_path", metavar="PATH", help="FASTQ reads or FASTA assembly")
    call.add_argument("--reads", dest="reads_path", metavar="FASTQ", help="Reads to map for assembly depth support")
    call.add_argument("--loci")
    call.add_argument("--primers", help="Primer-pair CSV/TSV/whitespace file with locus, forward, reverse columns")
    call.add_argument("--profiles")
    call.add_argument("--outdir", default="results")
    call.add_argument("--sample-id")
    call.add_argument("--min-read-length", type=int, default=50)
    call.add_argument("--max-read-length", type=int, default=100000)
    call.add_argument("--min-qscore", type=float, default=0.0)
    call.add_argument("--max-primer-mismatches", type=int, default=3)
    call.add_argument("--min-depth", type=int, default=10)
    call.add_argument("--min-posterior", type=float, default=0.75)
    call.add_argument("--threads", type=int, default=0, help="Worker threads; 0 uses all available CPUs")

    simulate = subparsers.add_parser("simulate", help="Simulate amplicon reads for a VNTR panel")
    simulate.add_argument("--loci", required=True)
    simulate.add_argument("--profile", dest="profiles")
    simulate.add_argument("--profile-id")
    simulate.add_argument("--sample-id", required=True)
    simulate.add_argument("--depth", type=int, default=200)
    simulate.add_argument("--error-rate", type=float, default=0.03)
    simulate.add_argument("--seed", type=int, default=13)
    simulate.add_argument("--outdir", required=True)

    extract = subparsers.add_parser(
        "extract-amplicons",
        help="Extract expected MLVA amplicons from FASTA/GenBank with amplirust",
    )
    extract.add_argument("--input", required=True, help="FASTA/GenBank input or amplirust-supported glob")
    extract.add_argument("--loci")
    extract.add_argument("--primers", help="Primer-pair CSV/TSV/whitespace file with locus, forward, reverse columns")
    extract.add_argument("--outdir", default="assembly_amplicons")
    extract.add_argument("--max-errors", type=int, default=2)
    extract.add_argument("--threads", type=int, default=0, help="amplirust threads; 0 lets amplirust auto-detect CPUs")
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
                threads=args.threads,
            )
            print(f"Wrote easy MLVA calls to {result['calls']}")
            print(f"Wrote detailed allele evidence to {result['allele_calls']}")
            print(f"Wrote report to {result['report']}")
        else:
            result = run_assembly_call(
                assembly_path=args.input_path,
                loci_path=args.loci,
                primers_path=args.primers,
                outdir=args.outdir,
                sample_id=sample_id,
                reads_path=args.reads_path,
                max_primer_mismatches=args.max_primer_mismatches,
                threads=args.threads,
            )
            print(f"Wrote easy MLVA calls to {result['calls']}")
            print(f"Wrote assembly amplicons to {result['amplicons']}")
            if args.reads_path:
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
