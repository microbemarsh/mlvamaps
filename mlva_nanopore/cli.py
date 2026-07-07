from __future__ import annotations

import argparse

from .in_silico_pcr import run_amplirust
from .pipeline import run_call
from .simulation import simulate_reads


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mlva-nanopore", description="Nanopore amplicon MLVA/VNTR typing")
    subparsers = parser.add_subparsers(dest="command", required=True)

    call = subparsers.add_parser("call", help="Call VNTR alleles from FASTQ reads")
    call.add_argument("--input", "--reads", dest="reads", required=True, metavar="INPUT")
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
        if not args.loci and not args.primers:
            parser.error("call requires either --loci or --primers")
        result = run_call(
            reads_path=args.reads,
            loci_path=args.loci,
            primers_path=args.primers,
            profiles_path=args.profiles,
            outdir=args.outdir,
            sample_id=args.sample_id or "sample",
            min_read_length=args.min_read_length,
            max_read_length=args.max_read_length,
            min_qscore=args.min_qscore,
            max_primer_mismatches=args.max_primer_mismatches,
            min_depth=args.min_depth,
            min_posterior=args.min_posterior,
            threads=args.threads,
        )
        print(f"Wrote MLVA calls to {result['allele_calls']}")
        print(f"Wrote report to {result['report']}")
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
