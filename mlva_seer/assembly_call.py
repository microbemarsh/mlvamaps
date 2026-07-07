from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from .calling import estimate_repeat_count_from_product_length
from .concurrency import resolve_threads
from .io import read_fasta, write_fasta, write_tsv
from .models import Locus
from .pipeline import SIMPLE_CALL_FIELDS
from .primers import read_loci_or_primers
from .sequence import find_best, revcomp


AMPLICON_FIELDS = [
    "sample_id",
    "locus_id",
    "product_id",
    "contig",
    "orientation",
    "product_size_bp",
    "forward_mismatches",
    "reverse_mismatches",
]

READ_SUPPORT_FIELDS = ["product_id", "locus_id", "mapped_reads", "mean_coverage"]


def _format_float(value: float | None, digits: int = 3) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def _find_product(sequence: str, locus: Locus, max_primer_mismatches: int) -> tuple[int, int, int, int] | None:
    f_pos, f_mm = find_best(locus.forward_primer, sequence, max_primer_mismatches)
    if f_pos is None:
        return None
    reverse_site = revcomp(locus.reverse_primer)
    search_start = f_pos + len(locus.forward_primer)
    r_rel, r_mm = find_best(reverse_site, sequence[search_start:], max_primer_mismatches)
    if r_rel is None:
        return None
    r_pos = search_start + r_rel
    end = r_pos + len(reverse_site)
    product_size = end - f_pos
    min_len = locus.expected_amplicon_min_bp or len(locus.forward_primer) + len(locus.reverse_primer)
    max_len = locus.expected_amplicon_max_bp or 100000
    if product_size < min_len or product_size > max_len:
        return None
    return f_pos, end, f_mm or 0, r_mm or 0


def extract_primer_products(
    assembly_path: str | Path,
    loci: list[Locus],
    sample_id: str,
    max_primer_mismatches: int = 3,
) -> list[dict]:
    products = []
    for contig_name, sequence in read_fasta(assembly_path):
        for locus in loci:
            candidates = [
                ("forward", sequence),
                ("reverse", revcomp(sequence)),
            ]
            for orientation, oriented_sequence in candidates:
                found = _find_product(oriented_sequence, locus, max_primer_mismatches)
                if found is None:
                    continue
                start, end, f_mm, r_mm = found
                product_id = f"{locus.locus_id}|{contig_name}|{orientation}|{start + 1}-{end}"
                products.append(
                    {
                        "sample_id": sample_id,
                        "locus_id": locus.locus_id,
                        "product_id": product_id,
                        "contig": contig_name,
                        "orientation": orientation,
                        "product_size_bp": end - start,
                        "forward_mismatches": f_mm,
                        "reverse_mismatches": r_mm,
                        "sequence": oriented_sequence[start:end],
                    }
                )
                break
    return products


def build_minimap2_command(
    reference_fasta: str | Path,
    reads_fastq: str | Path,
    threads: int = 0,
    preset: str | None = None,
) -> list[str]:
    command = [
        "minimap2",
        "-a",
        "-t",
        str(resolve_threads(threads)),
    ]
    if preset:
        command.extend(["-x", preset])
    command.extend([str(reference_fasta), str(reads_fastq)])
    return command


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _reference_bases_from_cigar(cigar: str) -> int:
    total = 0
    for length, op in _CIGAR_RE.findall(cigar):
        if op in {"M", "D", "N", "=", "X"}:
            total += int(length)
    return total


def read_minimap2_depth(sam_path: str | Path, reference_lengths: dict[str, int]) -> dict[str, dict[str, float]]:
    support = {name: {"mapped_reads": 0, "aligned_reference_bases": 0} for name in reference_lengths}
    with Path(sam_path).open() as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            flag = int(fields[1])
            if flag & 4 or flag & 256 or flag & 2048:
                continue
            reference_name = fields[2]
            if reference_name not in support:
                continue
            support[reference_name]["mapped_reads"] += 1
            support[reference_name]["aligned_reference_bases"] += _reference_bases_from_cigar(fields[5])
    return {
        name: {
            "mapped_reads": values["mapped_reads"],
            "mean_coverage": values["aligned_reference_bases"] / max(reference_lengths[name], 1),
        }
        for name, values in support.items()
    }


def run_minimap2_depth(
    reference_fasta: str | Path,
    reads_fastq: str | Path,
    sam_path: str | Path,
    reference_lengths: dict[str, int],
    threads: int = 0,
    minimap2_preset: str | None = None,
) -> dict[str, dict[str, float]]:
    if shutil.which("minimap2") is None:
        raise RuntimeError("Could not find 'minimap2' on PATH. Install minimap2 or rerun without --reads.")
    command = build_minimap2_command(reference_fasta, reads_fastq, threads, minimap2_preset)
    with Path(sam_path).open("w") as handle:
        subprocess.run(command, check=True, stdout=handle)
    return read_minimap2_depth(sam_path, reference_lengths)


def assembly_call_rows(
    loci: list[Locus],
    products: list[dict],
    sample_id: str,
    read_support: dict[str, dict[str, float]] | None = None,
) -> list[dict]:
    read_support = read_support or {}
    products_by_locus: dict[str, list[dict]] = {}
    for product in products:
        products_by_locus.setdefault(product["locus_id"], []).append(product)

    rows = []
    for locus in loci:
        locus_products = products_by_locus.get(locus.locus_id, [])
        if not locus_products:
            rows.append(
                {
                    "sample_id": sample_id,
                    "locus_id": locus.locus_id,
                    "present": "no",
                    "repeat_count": "",
                    "repeat_count_raw": "",
                    "product_size_bp": "",
                    "read_depth": 0,
                    "mean_coverage": "",
                    "status": "NOT_FOUND",
                    "evidence": "primer pair not found in assembly",
                }
            )
            continue

        product = sorted(locus_products, key=lambda item: (item["forward_mismatches"] + item["reverse_mismatches"], item["product_size_bp"]))[0]
        raw_count = estimate_repeat_count_from_product_length(locus, int(product["product_size_bp"]))
        support = read_support.get(product["product_id"], {})
        mapped_reads = int(support.get("mapped_reads", 0))
        mean_coverage = support.get("mean_coverage")
        status = "PASS" if raw_count is not None else "PRESENT_COUNT_UNKNOWN"
        rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "present": "yes",
                "repeat_count": round(raw_count) if raw_count is not None else "",
                "repeat_count_raw": _format_float(raw_count),
                "product_size_bp": product["product_size_bp"],
                "read_depth": mapped_reads,
                "mean_coverage": _format_float(mean_coverage),
                "status": status,
                "evidence": product["product_id"],
            }
        )
    return rows


def run_assembly_call(
    assembly_path: str,
    loci_path: str | None,
    outdir: str,
    sample_id: str,
    primers_path: str | None = None,
    reads_path: str | None = None,
    max_primer_mismatches: int = 3,
    threads: int = 0,
    minimap2_preset: str | None = None,
) -> dict[str, Path]:
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    loci = read_loci_or_primers(loci_path, primers_path)
    products = extract_primer_products(assembly_path, loci, sample_id, max_primer_mismatches)
    write_tsv(products, outdir_path / "assembly_amplicons.tsv", AMPLICON_FIELDS)

    product_fasta = outdir_path / "assembly_amplicons.fasta"
    write_fasta(((product["product_id"], product["sequence"]) for product in products), product_fasta)

    read_support = {}
    support_path = outdir_path / "read_support.tsv"
    if reads_path:
        support_rows = []
        if products:
            reference_lengths = {product["product_id"]: int(product["product_size_bp"]) for product in products}
            read_support = run_minimap2_depth(
                product_fasta,
                reads_path,
                outdir_path / "read_support.sam",
                reference_lengths,
                threads,
                minimap2_preset,
            )
            support_rows = [
                {
                    "product_id": product_id,
                    "locus_id": product_id.split("|", 1)[0],
                    "mapped_reads": int(row["mapped_reads"]),
                    "mean_coverage": _format_float(row["mean_coverage"]),
                }
                for product_id, row in read_support.items()
            ]
        write_tsv(support_rows, support_path, READ_SUPPORT_FIELDS)

    calls_path = outdir_path / "calls.tsv"
    write_tsv(assembly_call_rows(loci, products, sample_id, read_support), calls_path, SIMPLE_CALL_FIELDS)
    return {
        "outdir": outdir_path,
        "calls": calls_path,
        "amplicons": outdir_path / "assembly_amplicons.tsv",
        "amplicon_fasta": product_fasta,
        "read_support": support_path,
    }
