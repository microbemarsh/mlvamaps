from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from .calling import estimate_repeat_count_from_product_length
from .concurrency import resolve_threads
from .io import read_fasta, read_profiles, write_fasta, write_tsv
from .models import Locus
from .pipeline import MATCH_FIELDS, NOVELTY_FIELDS, SIMPLE_CALL_FIELDS
from .progress import ProgressReporter
from .profile_matching import build_fingerprint, match_profiles
from .primers import read_loci_or_primers
from .novelty import score_novelty
from .report import write_assembly_report
from .sequence import find_best, revcomp


AMPLICON_FIELDS = [
    "sample_id",
    "locus_id",
    "product_id",
    "contig",
    "contig_start",
    "contig_end",
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
    progress: ProgressReporter | None = None,
) -> list[dict]:
    progress = progress or ProgressReporter(enabled=False)
    products = []
    for contig_idx, (contig_name, sequence) in enumerate(read_fasta(assembly_path), start=1):
        contig_len = len(sequence)
        progress.step(f"Scanning contig {contig_idx:,}: {contig_name} ({len(sequence):,} bp)")
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
                if orientation == "forward":
                    contig_start = start + 1
                    contig_end = end
                else:
                    contig_start = contig_len - end + 1
                    contig_end = contig_len - start
                product_id = f"{locus.locus_id}|{contig_name}|{orientation}|{contig_start}-{contig_end}"
                products.append(
                    {
                        "sample_id": sample_id,
                        "locus_id": locus.locus_id,
                        "product_id": product_id,
                        "contig": contig_name,
                        "contig_start": contig_start,
                        "contig_end": contig_end,
                        "orientation": orientation,
                        "product_size_bp": end - start,
                        "forward_mismatches": f_mm,
                        "reverse_mismatches": r_mm,
                        "sequence": oriented_sequence[start:end],
                    }
                )
                break
        progress.count("Scanned assembly contigs", contig_idx)
    progress.step(f"Found {len(products):,} primer product(s)")
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
    progress: ProgressReporter | None = None,
) -> dict[str, dict[str, float]]:
    if shutil.which("minimap2") is None:
        raise RuntimeError("Could not find 'minimap2' on PATH. Install minimap2 or rerun without --reads.")
    progress = progress or ProgressReporter(enabled=False)
    command = build_minimap2_command(reference_fasta, reads_fastq, threads, minimap2_preset)
    progress.step("Running minimap2 for read-depth support")
    with Path(sam_path).open("w") as handle:
        subprocess.run(command, check=True, stdout=handle)
    progress.step("Parsing minimap2 alignments")
    return read_minimap2_depth(sam_path, reference_lengths)


def _alignment_blocks_from_cigar(pos_1based: int, cigar: str) -> list[tuple[int, int]]:
    blocks = []
    ref_pos = pos_1based
    block_start = None
    for length_text, op in _CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op in {"M", "=", "X"}:
            if block_start is None:
                block_start = ref_pos
            ref_pos += length
        else:
            if block_start is not None:
                blocks.append((block_start, ref_pos - 1))
                block_start = None
            if op in {"D", "N"}:
                ref_pos += length
    if block_start is not None:
        blocks.append((block_start, ref_pos - 1))
    return blocks


def _overlap_bases(blocks: list[tuple[int, int]], start: int, end: int) -> int:
    total = 0
    for block_start, block_end in blocks:
        overlap_start = max(block_start, start)
        overlap_end = min(block_end, end)
        if overlap_start <= overlap_end:
            total += overlap_end - overlap_start + 1
    return total


def read_alignment_depth(alignment_path: str | Path, products: list[dict]) -> dict[str, dict[str, float]]:
    path = Path(alignment_path)
    if path.suffix.lower() == ".bam":
        return _read_bam_depth(path, products)
    return _read_sam_depth(path, products)


def _read_sam_depth(sam_path: Path, products: list[dict]) -> dict[str, dict[str, float]]:
    support = {
        product["product_id"]: {"mapped_reads": 0, "aligned_reference_bases": 0}
        for product in products
    }
    products_by_contig: dict[str, list[dict]] = {}
    for product in products:
        products_by_contig.setdefault(product["contig"], []).append(product)
    with sam_path.open() as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            flag = int(fields[1])
            if flag & 4 or flag & 256 or flag & 2048:
                continue
            contig = fields[2]
            if contig not in products_by_contig:
                continue
            blocks = _alignment_blocks_from_cigar(int(fields[3]), fields[5])
            for product in products_by_contig[contig]:
                overlap = _overlap_bases(blocks, int(product["contig_start"]), int(product["contig_end"]))
                if overlap:
                    support[product["product_id"]]["mapped_reads"] += 1
                    support[product["product_id"]]["aligned_reference_bases"] += overlap
    return _finalize_product_depth(support, products)


def _read_bam_depth(bam_path: Path, products: list[dict]) -> dict[str, dict[str, float]]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - depends on optional package
        raise RuntimeError("Reading BAM support requires pysam. Install pysam or provide SAM/FASTQ input.") from exc

    support = {
        product["product_id"]: {"mapped_reads": 0, "aligned_reference_bases": 0}
        for product in products
    }
    with pysam.AlignmentFile(bam_path) as bam:
        for product in products:
            product_id = product["product_id"]
            start_0based = int(product["contig_start"]) - 1
            end_0based = int(product["contig_end"])
            try:
                alignments = bam.fetch(product["contig"], start_0based, end_0based)
            except ValueError:
                continue
            for alignment in alignments:
                if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                    continue
                blocks = [(start + 1, end) for start, end in alignment.get_blocks()]
                overlap = _overlap_bases(blocks, int(product["contig_start"]), int(product["contig_end"]))
                if overlap:
                    support[product_id]["mapped_reads"] += 1
                    support[product_id]["aligned_reference_bases"] += overlap
    return _finalize_product_depth(support, products)


def _finalize_product_depth(
    support: dict[str, dict[str, float]],
    products: list[dict],
) -> dict[str, dict[str, float]]:
    lengths = {product["product_id"]: int(product["product_size_bp"]) for product in products}
    return {
        product_id: {
            "mapped_reads": values["mapped_reads"],
            "mean_coverage": values["aligned_reference_bases"] / max(lengths[product_id], 1),
        }
        for product_id, values in support.items()
    }


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


def allele_rows_from_assembly_calls(call_rows: list[dict]) -> list[dict]:
    rows = []
    for row in call_rows:
        present = row.get("present") == "yes"
        status = row.get("status", "")
        rows.append(
            {
                "sample_id": row.get("sample_id", ""),
                "locus_id": row.get("locus_id", ""),
                "called_repeat_count": row.get("repeat_count", "") if present else "",
                "posterior_probability": 1.0 if status == "PASS" else 0.0,
                "second_best_repeat_count": "",
                "second_best_posterior": 0.0,
                "read_depth": row.get("read_depth", 0),
                "num_vntr_asvs": 1 if present else 0,
                "dominant_vntr_asv": row.get("evidence", "") if present else "",
                "call_status": status,
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
    alignments_path: str | None = None,
    profiles_path: str | None = None,
    max_primer_mismatches: int = 3,
    threads: int = 0,
    minimap2_preset: str | None = None,
    show_progress: bool = False,
) -> dict[str, Path]:
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    progress = ProgressReporter(enabled=show_progress)
    thread_count = resolve_threads(threads)
    progress.step(f"Starting assembly call for sample {sample_id!r} with {thread_count} worker(s)")

    progress.step("Loading panel")
    loci = read_loci_or_primers(loci_path, primers_path)
    profiles = read_profiles(profiles_path)
    progress.step(f"Loaded {len(loci):,} loci" + (f" and {len(profiles):,} reference profiles" if profiles else ""))
    products = extract_primer_products(assembly_path, loci, sample_id, max_primer_mismatches, progress)
    write_tsv(products, outdir_path / "assembly_amplicons.tsv", AMPLICON_FIELDS)

    progress.step("Writing assembly amplicon FASTA")
    product_fasta = outdir_path / "assembly_amplicons.fasta"
    write_fasta(((product["product_id"], product["sequence"]) for product in products), product_fasta)

    read_support = {}
    support_path = outdir_path / "read_support.tsv"
    if reads_path or alignments_path:
        support_rows = []
        if products and reads_path:
            reference_lengths = {product["product_id"]: int(product["product_size_bp"]) for product in products}
            read_support = run_minimap2_depth(
                product_fasta,
                reads_path,
                outdir_path / "read_support.sam",
                reference_lengths,
                threads,
                minimap2_preset,
                progress,
            )
        elif products and alignments_path:
            progress.step("Reading alignment depth support")
            read_support = read_alignment_depth(alignments_path, products)
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
    progress.step("Writing calls table")
    call_rows = assembly_call_rows(loci, products, sample_id, read_support)
    write_tsv(call_rows, calls_path, SIMPLE_CALL_FIELDS)
    allele_rows = allele_rows_from_assembly_calls(call_rows)
    fingerprint_rows, probabilistic_rows = build_fingerprint(sample_id, allele_rows, loci)
    fingerprint_fields = ["sample_id"] + [locus.locus_id for locus in loci]
    write_tsv(fingerprint_rows, outdir_path / "mlva_fingerprint.tsv", fingerprint_fields)
    write_tsv(
        probabilistic_rows,
        outdir_path / "mlva_fingerprint_probabilistic.tsv",
        ["sample_id", "locus_id", "repeat_count", "posterior_probability"],
    )
    match_rows = match_profiles(sample_id, fingerprint_rows[0], profiles)
    write_tsv(match_rows, outdir_path / "profile_matches.tsv", MATCH_FIELDS)
    novelty_rows = score_novelty(sample_id, allele_rows, match_rows)
    write_tsv(novelty_rows, outdir_path / "novelty_scores.tsv", NOVELTY_FIELDS)
    progress.step("Writing HTML report")
    write_assembly_report(outdir_path, sample_id, call_rows, products, match_rows, profiles, novelty_rows)
    progress.step(f"Done. Main calls: {calls_path}")
    return {
        "outdir": outdir_path,
        "calls": calls_path,
        "amplicons": outdir_path / "assembly_amplicons.tsv",
        "amplicon_fasta": product_fasta,
        "read_support": support_path,
        "fingerprint": outdir_path / "mlva_fingerprint.tsv",
        "profile_matches": outdir_path / "profile_matches.tsv",
        "report": outdir_path / "report.html",
    }
