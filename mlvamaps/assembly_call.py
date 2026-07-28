from __future__ import annotations

import csv
import math
import re
from pathlib import Path

from .calling import (
    allele_grid,
    assembly_equivalent_product_allele,
    estimate_repeat_count_from_product_length,
    gaussian_allele_probabilities,
    legacy_round_repeat_count,
    repeat_unit_length,
)
from .concurrency import DEFAULT_THREADS, resolve_threads
from .in_silico_pcr import read_pcr_results, run_in_silico_pcr_loci
from .io import read_fasta, read_profiles, write_fasta, write_tsv
from .locus_measurement import measure_locus_product
from .models import Locus
from .mapping import (
    build_minimap2_map_command,
    check_minimap2,
    run_minimap2_command,
)
from .pipeline import (
    ALLELE_DISTRIBUTION_FIELDS,
    MATCH_FIELDS,
    REPEAT_COUNT_FIELDS,
    SIMPLE_CALL_FIELDS,
    allele_distribution_rows,
)
from .progress import ProgressReporter
from .profile_matching import (
    PROFILE_MATCH_LOCUS_FIELDS,
    build_fingerprint,
    match_profiles,
    profile_match_locus_rows,
    sequence_reference_match_rows,
)
from .phylogeny import run_phylogenetic_placement
from .primers import LEGACY_NAME_RE, read_loci_or_primers
from .report import write_assembly_report


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
    "primer_error_round",
    "size_derived_repeat_count",
    "forward_primer_match",
    "reverse_primer_match",
]

READ_SUPPORT_FIELDS = ["product_id", "locus_id", "mapped_reads", "mean_coverage"]

ASSEMBLY_ALGORITHMS = ("legacy", "novel")

LEGACY_DETAIL_FIELDS = [
    "strain",
    "primer",
    "position1",
    "position2",
    "size",
    "allele",
    "sequence",
    "nb_mismatch",
    "primer1",
    "mismatch1",
    "primer2",
    "mismatch2",
    "predicted PCR target",
]


_IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def _primer_match_display(primer: str, observed: str) -> str:
    """Render an observed primer match using the legacy case convention.

    Substitutions and insertions are lower-case and deleted primer bases are
    dots. Amplirust normally reports equal-length matches here; handling the
    length difference also keeps the compatibility output useful for indels.
    """
    primer = primer.upper()
    observed = observed.upper()
    rows = len(primer) + 1
    cols = len(observed) + 1
    distance = [[0] * cols for _ in range(rows)]
    trace = [[""] * cols for _ in range(rows)]
    for index in range(1, rows):
        distance[index][0] = index
        trace[index][0] = "deletion"
    for index in range(1, cols):
        distance[0][index] = index
        trace[0][index] = "insertion"
    for primer_index in range(1, rows):
        for observed_index in range(1, cols):
            compatible = observed[observed_index - 1] in _IUPAC_BASES.get(
                primer[primer_index - 1], frozenset(primer[primer_index - 1])
            )
            choices = (
                (distance[primer_index - 1][observed_index - 1] + (not compatible), "diagonal"),
                (distance[primer_index - 1][observed_index] + 1, "deletion"),
                (distance[primer_index][observed_index - 1] + 1, "insertion"),
            )
            distance[primer_index][observed_index], trace[primer_index][observed_index] = min(
                choices, key=lambda choice: choice[0]
            )

    rendered = []
    primer_index = len(primer)
    observed_index = len(observed)
    while primer_index or observed_index:
        operation = trace[primer_index][observed_index]
        if operation == "diagonal":
            primer_base = primer[primer_index - 1]
            observed_base = observed[observed_index - 1]
            compatible = observed_base in _IUPAC_BASES.get(
                primer_base, frozenset(primer_base)
            )
            rendered.append(observed_base if compatible else observed_base.lower())
            primer_index -= 1
            observed_index -= 1
        elif operation == "deletion":
            rendered.append(".")
            primer_index -= 1
        else:
            rendered.append(observed[observed_index - 1].lower())
            observed_index -= 1
    return "".join(reversed(rendered))


def _format_float(value: float | None, digits: int = 3) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def legacy_amplicon_bounds(loci: list[Locus]) -> tuple[int, int]:
    """Return global Amplirust bounds covering every MLVA_finder-valid hit."""
    bounds = []
    for locus in loci:
        repeat_bp = repeat_unit_length(locus)
        if repeat_bp and locus.expected_product_size_bp and locus.nominal_repeat_units:
            lower_exclusive = locus.expected_product_size_bp - repeat_bp * (
                locus.nominal_repeat_units + 100
            )
            upper_exclusive = locus.expected_product_size_bp + repeat_bp * (
                100 - locus.nominal_repeat_units
            )
            lower = max(1, math.floor(lower_exclusive) + 1)
            upper = max(lower, math.ceil(upper_exclusive) - 1)
        else:
            lower = max(1, locus.expected_amplicon_min_bp)
            upper = max(lower, locus.expected_amplicon_max_bp or 100000)
        bounds.append((lower, upper))
    if not bounds:
        return (1, 100000)
    return min(lower for lower, _ in bounds), max(upper for _, upper in bounds)


def _product_in_locus_bounds(product: dict, locus: Locus) -> bool:
    product_size = int(product["product_size_bp"])
    minimum = locus.expected_amplicon_min_bp or len(locus.forward_primer) + len(
        locus.reverse_primer
    )
    maximum = locus.expected_amplicon_max_bp or 100000
    return minimum <= product_size <= maximum


def pcr_rows_to_products(
    rows: list[dict[str, str | int]],
    loci: list[Locus],
    sample_id: str,
    enforce_locus_bounds: bool = True,
    reference_order: dict[str, int] | None = None,
    measure_products: bool = True,
) -> list[dict]:
    """Convert native PCR output into MLVAMaps assembly-product records.

    ``MLVA_finder`` reports product size using the configured primer lengths,
    even when a fuzzy primer match contains an insertion or deletion.  That is
    subtly different from ``full_len`` (the observed sequence
    span), so retain the latter internally and expose the legacy-compatible
    size to the assembly callers.
    """
    locus_by_id = {locus.locus_id: locus for locus in loci}
    reference_order = reference_order or {}
    products = []
    for row in rows:
        locus_id = str(row.get("primer_name", ""))
        locus = locus_by_id.get(locus_id)
        if locus is None or str(row.get("is_circular_wrap", "false")).lower() == "true":
            continue
        amplicon_span_size = int(row["full_len"])
        fwd_start = int(row["fwd_start"])
        fwd_end = int(row["fwd_end"])
        rev_start = int(row["rev_start"])
        rev_end = int(row["rev_end"])
        forward_match_length = max(0, fwd_end - fwd_start)
        reverse_match_length = max(0, rev_end - rev_start)
        # This is algebraically equivalent to both orientation-specific size
        # formulas in MLVA_finder.py.  In particular, the configured reverse
        # primer length is used rather than the observed reverse-match length.
        product_size = (
            rev_start
            + len(locus.reverse_primer)
            - fwd_start
            - forward_match_length
            + len(locus.forward_primer)
        )
        if product_size <= 0:
            continue
        min_len = locus.expected_amplicon_min_bp or len(locus.forward_primer) + len(locus.reverse_primer)
        max_len = locus.expected_amplicon_max_bp or 100000
        if enforce_locus_bounds and not min_len <= product_size <= max_len:
            continue
        contig = str(row["reference_id"]).split()[0]
        contig_start = int(row["original_start"]) + 1
        contig_end = int(row["original_end"])
        if contig_end < contig_start:
            continue
        orientation = "reverse" if row.get("strand") == "-" else "forward"
        product_id = f"{locus_id}|{contig}|{orientation}|{contig_start}-{contig_end}"
        sequence = str(row["product_seq"]).upper()
        original_start = int(row["original_start"])
        original_end = int(row["original_end"])
        if orientation == "forward":
            forward_start_0based = original_start
            reverse_start_0based = original_end - reverse_match_length
        else:
            forward_start_0based = original_end - forward_match_length
            reverse_start_0based = original_start
        measurement = (
            measure_locus_product(
                sequence,
                locus,
                source="assembly",
                sequence_id=product_id,
                calibrated_product_size_bp=product_size,
            )
            if measure_products
            else None
        )
        raw_count = (
            measurement.raw_repeat_count if measurement is not None else None
        )
        if raw_count is None:
            raw_count = estimate_repeat_count_from_product_length(
                locus, product_size
            )
        products.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "product_id": product_id,
                "contig": contig,
                "contig_start": contig_start,
                "contig_end": contig_end,
                "orientation": orientation,
                "product_size_bp": product_size,
                "amplicon_span_size_bp": amplicon_span_size,
                "forward_mismatches": int(row["fwd_mismatches"]),
                "reverse_mismatches": int(row["rev_mismatches"]),
                "primer_error_round": max(int(row["fwd_mismatches"]), int(row["rev_mismatches"])),
                "size_derived_repeat_count": _format_float(raw_count),
                "forward_primer_match": sequence[:forward_match_length],
                "reverse_primer_match": _reverse_complement(sequence[-reverse_match_length:])
                if reverse_match_length
                else "",
                "forward_start_0based": forward_start_0based,
                "reverse_start_0based": reverse_start_0based,
                "forward_match_length": forward_match_length,
                "reverse_match_length": reverse_match_length,
                "contig_index": reference_order.get(contig, 0),
                "sequence": sequence,
                "measurement_status": measurement.status if measurement else "",
                "measured_repeat_length_bp": (
                    measurement.repeat_length_bp if measurement else ""
                ),
                "measured_raw_repeat_count": (
                    measurement.raw_repeat_count if measurement else ""
                ),
                "measured_allele": measurement.called_allele if measurement else "",
            }
        )
    return products


# Source compatibility for the MLVAMaps 0.1 API.
amplirust_rows_to_products = pcr_rows_to_products


def build_minimap2_command(
    reference_fasta: str | Path,
    reads_fastq: str | Path,
    threads: int = DEFAULT_THREADS,
    preset: str | None = None,
    executable: str = "minimap2",
) -> list[str]:
    return build_minimap2_map_command(
        reference_fasta,
        reads_fastq,
        resolve_threads(threads),
        executable=executable,
        preset=preset,
    )


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _reference_bases_from_cigar(cigar: str) -> int:
    total = 0
    for length, op in _CIGAR_RE.findall(cigar):
        if op in {"M", "D", "N", "=", "X"}:
            total += int(length)
    return total


def read_minimap2_depth(
    sam_path: str | Path, reference_lengths: dict[str, int]
) -> dict[str, dict[str, float]]:
    support = {
        name: {"mapped_reads": 0, "aligned_reference_bases": 0}
        for name in reference_lengths
    }
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
    threads: int = DEFAULT_THREADS,
    minimap2_preset: str | None = None,
    executable: str = "minimap2",
    progress: ProgressReporter | None = None,
) -> dict[str, dict[str, float]]:
    progress = progress or ProgressReporter(enabled=False)
    executable_path = check_minimap2(executable)
    command = build_minimap2_command(
        reference_fasta,
        reads_fastq,
        threads,
        minimap2_preset,
        executable_path,
    )
    progress.step("Running minimap2 for read-depth support")
    run_minimap2_command(
        command,
        "assembly read-support mapping",
        stdout_path=Path(sam_path),
    )
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


def _not_found_row(sample_id: str, locus_id: str) -> dict:
    return {
        "sample_id": sample_id,
        "locus_id": locus_id,
        "present": "no",
        "repeat_count": "",
        "repeat_count_raw": "",
        "product_size_bp": "",
        "read_depth": 0,
        "mean_coverage": "",
        "status": "NOT_FOUND",
        "evidence": "primer pair not found in assembly",
    }


def legacy_assembly_call_rows(
    loci: list[Locus],
    products: list[dict],
    sample_id: str,
    round_tolerance: float = 0.25,
    read_support: dict[str, dict[str, float]] | None = None,
) -> list[dict]:
    """Call assembly alleles with MLVA_finder's historical decision rules.

    The original script searched successively larger per-primer mismatch
    rounds. Within the first successful round, each FASTA record replaced the
    preceding record's result, and the smallest unrounded allele on that final
    matching record was retained. It also preferred a forward-orientation
    first-primer match over reverse orientation and equal-length fuzzy primer
    matches over indel matches. Values of 100 or greater were rejected before
    rounding.
    """
    read_support = read_support or {}
    products_by_locus: dict[str, list[dict]] = {}
    for product in products:
        products_by_locus.setdefault(product["locus_id"], []).append(product)

    rows = []
    for locus in loci:
        locus_products = products_by_locus.get(locus.locus_id, [])
        eligible = []
        for product in locus_products:
            raw_count, _called_count = assembly_equivalent_product_allele(
                locus,
                int(product["product_size_bp"]),
                round_tolerance,
            )
            if raw_count is not None and raw_count < 100:
                eligible.append((product, raw_count))
        if not eligible:
            rows.append(_not_found_row(sample_id, locus.locus_id))
            continue

        def compatible_with_legacy_search(product: dict, error_round: int) -> bool:
            same_search = [
                other
                for other in locus_products
                if other.get("contig") == product.get("contig")
                and other.get("orientation") == product.get("orientation")
            ]
            # MLVA_finder searches the forward primer on the input strand
            # first. Any such match prevents its reverse-complement fallback.
            if product.get("orientation") == "reverse" and any(
                other.get("orientation") == "forward"
                and other.get("contig") == product.get("contig")
                and int(other.get("forward_mismatches", 0)) <= error_round
                for other in locus_products
            ):
                return False
            # Its fuzzy matcher discards indel-length matches whenever an
            # equal-length match exists at the same mismatch threshold.
            for length_key, error_key, configured_length in (
                ("forward_match_length", "forward_mismatches", len(locus.forward_primer)),
                ("reverse_match_length", "reverse_mismatches", len(locus.reverse_primer)),
            ):
                if int(product.get(length_key, configured_length)) == configured_length:
                    continue
                if any(
                    int(other.get(error_key, 0)) <= error_round
                    and int(other.get(length_key, configured_length)) == configured_length
                    for other in same_search
                ):
                    return False
            return True

        round_products = []
        max_error_round = max(int(product["primer_error_round"]) for product, _ in eligible)
        for first_error_round in range(max_error_round + 1):
            round_products = [
                item
                for item in eligible
                if int(item[0]["primer_error_round"]) <= first_error_round
                and compatible_with_legacy_search(item[0], first_error_round)
            ]
            if round_products:
                break
        if not round_products:
            rows.append(_not_found_row(sample_id, locus.locus_id))
            continue
        final_contig_index = max(
            int(product.get("contig_index", 0)) for product, _ in round_products
        )
        product, raw_count = min(
            (
                (product, raw_count)
                for product, raw_count in round_products
                if int(product.get("contig_index", 0)) == final_contig_index
            ),
            key=lambda item: item[1],
        )
        _raw_count, called_count = assembly_equivalent_product_allele(
            locus,
            int(product["product_size_bp"]),
            round_tolerance,
        )
        support_probability = 1.0
        support = read_support.get(product["product_id"], {})
        rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "present": "yes",
                "repeat_count": called_count,
                "repeat_count_raw": _format_float(raw_count),
                "product_size_bp": product["product_size_bp"],
                "read_depth": int(support.get("mapped_reads", 0)),
                "mean_coverage": _format_float(support.get("mean_coverage")),
                "allele_confidence": support_probability,
                "second_best_repeat_count": "",
                "second_best_probability": 0.0,
                "inference_method": "legacy_minimum_allele",
                "allele_distribution": f"{called_count}:{support_probability:.6f}",
                "status": "PASS",
                "evidence": product["product_id"],
            }
        )
    return rows


def novel_assembly_call_rows(
    loci: list[Locus],
    products: list[dict],
    sample_id: str,
    read_support: dict[str, dict[str, float]] | None = None,
    min_posterior: float = 0.75,
) -> list[dict]:
    if not 0 <= min_posterior <= 1:
        raise ValueError("minimum allele posterior must be between 0 and 1")
    read_support = read_support or {}
    products_by_locus: dict[str, list[dict]] = {}
    for product in products:
        products_by_locus.setdefault(product["locus_id"], []).append(product)

    rows = []
    for locus in loci:
        locus_products = [
            product
            for product in products_by_locus.get(locus.locus_id, [])
            if _product_in_locus_bounds(product, locus)
        ]
        if not locus_products:
            rows.append(_not_found_row(sample_id, locus.locus_id))
            continue

        legacy_ranked = sorted(
            locus_products,
            key=lambda item: (
                item["primer_error_round"],
                float(item["size_derived_repeat_count"] or "inf"),
                item["forward_mismatches"] + item["reverse_mismatches"],
                item["product_size_bp"],
                item["product_id"],
            ),
        )
        product = legacy_ranked[0]
        candidates = allele_grid(locus, step=0.5)
        raw_by_product = {
            item["product_id"]: estimate_repeat_count_from_product_length(
                locus, int(item["product_size_bp"])
            )
            for item in locus_products
        }
        count_known_products = [
            item
            for item in locus_products
            if raw_by_product[item["product_id"]] is not None
        ]
        has_depth_support = any(
            int(read_support.get(item["product_id"], {}).get("mapped_reads", 0)) > 0
            for item in count_known_products
        )
        observations = count_known_products if has_depth_support else (
            [product] if raw_by_product[product["product_id"]] is not None else []
        )
        weights = {
            item["product_id"]: (
                float(read_support.get(item["product_id"], {}).get("mapped_reads", 0))
                + 0.5
                * math.exp(
                    -(item["forward_mismatches"] + item["reverse_mismatches"])
                )
            )
            if has_depth_support
            else 1.0
            for item in observations
        }
        allele_weights = {candidate: 0.0 for candidate in candidates}
        sigma = max(0.08, 0.5 / max(repeat_unit_length(locus), 1))
        for item in observations:
            raw = raw_by_product[item["product_id"]]
            probabilities = gaussian_allele_probabilities(raw, candidates, sigma)
            for candidate, probability in zip(candidates, probabilities):
                allele_weights[candidate] += probability * weights[item["product_id"]]
        allele_total = sum(allele_weights.values())
        allele_ranking = sorted(
            (
                (candidate, weight / allele_total)
                for candidate, weight in allele_weights.items()
            ),
            key=lambda item: (-item[1], float(item[0])),
        ) if allele_total else []
        if allele_ranking:
            called_count, allele_confidence = allele_ranking[0]
            second_count, second_probability = (
                allele_ranking[1] if len(allele_ranking) > 1 else ("", 0.0)
            )
            product = min(
                observations,
                key=lambda item: (
                    -weights[item["product_id"]],
                    abs(float(raw_by_product[item["product_id"]]) - float(called_count)),
                    item["primer_error_round"],
                    item["product_id"],
                ),
            )
        else:
            called_count = ""
            allele_confidence = 0.0
            second_count = ""
            second_probability = 0.0
        raw_count = estimate_repeat_count_from_product_length(locus, int(product["product_size_bp"]))
        support = read_support.get(product["product_id"], {})
        mapped_reads = int(support.get("mapped_reads", 0))
        mean_coverage = support.get("mean_coverage")
        status = "PASS" if raw_count is not None else "PRESENT_COUNT_UNKNOWN"
        if raw_count is not None and (
            allele_confidence < min_posterior
            or allele_confidence - second_probability < 0.2
        ):
            status = "AMBIGUOUS"
        inference_method = (
            "depth_weighted_product_distribution"
            if has_depth_support
            else "assembly_product_length"
        )
        rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "present": "yes",
                "repeat_count": called_count,
                "repeat_count_raw": _format_float(raw_count),
                "product_size_bp": product["product_size_bp"],
                "read_depth": mapped_reads,
                "mean_coverage": _format_float(mean_coverage),
                "allele_confidence": round(allele_confidence, 6),
                "second_best_repeat_count": second_count,
                "second_best_probability": round(second_probability, 6),
                "inference_method": inference_method,
                "allele_distribution": ";".join(
                    f"{candidate}:{probability:.6f}"
                    for candidate, probability in allele_ranking
                ),
                "status": status,
                "evidence": product["product_id"],
            }
        )
    return rows


def assembly_call_rows(
    loci: list[Locus],
    products: list[dict],
    sample_id: str,
    read_support: dict[str, dict[str, float]] | None = None,
    min_posterior: float = 0.75,
    algorithm: str = "legacy",
    round_tolerance: float = 0.25,
) -> list[dict]:
    """Dispatch to the selected assembly allele-calling algorithm."""
    if algorithm == "legacy":
        return legacy_assembly_call_rows(
            loci,
            products,
            sample_id,
            round_tolerance=round_tolerance,
            read_support=read_support,
        )
    if algorithm == "novel":
        return novel_assembly_call_rows(
            loci, products, sample_id, read_support, min_posterior
        )
    choices = ", ".join(ASSEMBLY_ALGORITHMS)
    raise ValueError(f"unknown assembly algorithm {algorithm!r}; choose one of: {choices}")


def legacy_detail_rows(
    loci: list[Locus],
    products: list[dict],
    call_rows: list[dict],
    sample_id: str,
    round_tolerance: float = 0.25,
) -> list[dict]:
    """Build the row-oriented table emitted by the historical caller."""
    products_by_id = {product["product_id"]: product for product in products}
    loci_by_id = {locus.locus_id: locus for locus in loci}
    rows = []
    for call in call_rows:
        locus = loci_by_id[call["locus_id"]]
        product = products_by_id.get(call.get("evidence", ""))
        if product is None:
            rows.append(
                {
                    "strain": sample_id,
                    "primer": locus.locus_id,
                    "allele": "",
                    "nb_mismatch": "ND",
                    "primer1": locus.forward_primer,
                    "primer2": locus.reverse_primer,
                }
            )
            continue
        forward_errors = int(product["forward_mismatches"])
        reverse_errors = int(product["reverse_mismatches"])
        forward_match = product["forward_primer_match"]
        reverse_match = product["reverse_primer_match"]
        raw_count = estimate_repeat_count_from_product_length(
            locus, int(product["product_size_bp"])
        )
        rows.append(
            {
                "strain": sample_id,
                "primer": locus.locus_id,
                "position1": product["forward_start_0based"],
                "position2": product["reverse_start_0based"],
                "size": product["product_size_bp"],
                "allele": legacy_round_repeat_count(raw_count, round_tolerance)
                if raw_count is not None
                else "",
                "sequence": product["contig"],
                "nb_mismatch": product["primer_error_round"],
                "primer1": locus.forward_primer,
                "mismatch1": _primer_match_display(locus.forward_primer, forward_match)
                if forward_errors
                else "",
                "primer2": locus.reverse_primer,
                "mismatch2": _primer_match_display(locus.reverse_primer, reverse_match)
                if reverse_errors
                else "",
                "predicted PCR target": product["sequence"],
            }
        )
    return rows


def write_legacy_assembly_outputs(
    outdir: str | Path,
    loci: list[Locus],
    products: list[dict],
    call_rows: list[dict],
    sample_id: str,
    round_tolerance: float = 0.25,
) -> dict[str, Path]:
    """Write stable CSV equivalents of all four useful legacy artifacts."""
    outdir = Path(outdir)
    details_path = outdir / "legacy_output.csv"
    fingerprint_path = outdir / "legacy_mlva_analysis.csv"
    sizes_path = outdir / "legacy_predicted_pcr_sizes.csv"
    mismatches_path = outdir / "legacy_primer_mismatches.txt"
    details = legacy_detail_rows(
        loci, products, call_rows, sample_id, round_tolerance
    )
    with details_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEGACY_DETAIL_FIELDS)
        writer.writeheader()
        writer.writerows(details)

    # MLVA_finder shortened encoded names (for example
    # ``Lp03_96bp_941bp_8U``) to their first underscore-delimited field in
    # the wide analysis tables, while retaining the full name in output.csv.
    locus_ids = [
        locus.locus_id.split("_", 1)[0]
        if LEGACY_NAME_RE.match(locus.locus_id)
        else locus.locus_id
        for locus in loci
    ]
    with fingerprint_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "Access_number", *locus_ids])
        writer.writerow(["001", sample_id, *[row["allele"] for row in details]])
    with sizes_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key", "Access_number", *locus_ids])
        writer.writerow(["001", sample_id, *[row.get("size", "") for row in details]])

    mismatch_sections = []
    for row in details:
        for suffix, primer_key, mismatch_key in (
            ("FOR", "primer1", "mismatch1"),
            ("REV", "primer2", "mismatch2"),
        ):
            if row.get(mismatch_key):
                mismatch_sections.append(
                    f"{row['primer']}_{suffix}\n{row[primer_key]}\n{row[mismatch_key]}"
                )
    mismatches_path.write_text("\n\n".join(mismatch_sections) + ("\n" if mismatch_sections else ""))
    return {
        "legacy_details": details_path,
        "legacy_fingerprint": fingerprint_path,
        "legacy_product_sizes": sizes_path,
        "legacy_mismatches": mismatches_path,
    }


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
                "posterior_probability": row.get("allele_confidence", 0.0),
                "second_best_repeat_count": row.get("second_best_repeat_count", ""),
                "second_best_posterior": row.get("second_best_probability", 0.0),
                "read_depth": row.get("read_depth", 0),
                "num_vntr_asvs": 1 if present else 0,
                "dominant_vntr_asv": row.get("evidence", "") if present else "",
                "allele_distribution": row.get("allele_distribution", ""),
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
    database_path: str | None = None,
    max_primer_mismatches: int = 2,
    assembly_round_tolerance: float = 0.25,
    algorithm: str = "legacy",
    min_posterior: float = 0.75,
    threads: int = DEFAULT_THREADS,
    minimap2_preset: str | None = None,
    minimap2_bin: str = "minimap2",
    amplirust_bin: str = "amplirust",
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    epa_ng_bin: str = "epa-ng",
    dnadiff_bin: str = "dnadiff",
    raxml_model: str = "DNA",
    phylogeny_snp_weight: float = 1.0,
    phylogeny_repeat_weight: float = 1.0,
    reference_metadata_path: str | None = None,
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
    progress.step("Finding MLVA_finder-compatible primer products with Sassy")
    pcr_paths = run_in_silico_pcr_loci(
        assembly_path,
        loci,
        outdir_path / "in_silico_pcr",
        max_errors=max_primer_mismatches,
        threads=threads,
        # MLVA_finder does not reject an otherwise valid product because its
        # interior crosses an assembly gap represented by N bases.
        max_n_fraction=1.0,
        amplicon_bounds=legacy_amplicon_bounds(loci),
    )
    reference_order = {
        reference_id: index
        for index, (reference_id, _sequence) in enumerate(read_fasta(assembly_path))
    }
    products = pcr_rows_to_products(
        read_pcr_results(pcr_paths["stats"], pcr_paths["products"]),
        loci,
        sample_id,
        enforce_locus_bounds=False,
        reference_order=reference_order,
        # The historical caller derives its allele exclusively from the
        # calibrated product length. Avoid a second anchor search whose result
        # cannot affect the legacy call.
        measure_products=algorithm != "legacy",
    )
    progress.step(f"Found {len(products):,} primer product(s)")
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
                threads=threads,
                minimap2_preset=minimap2_preset,
                executable=minimap2_bin,
                progress=progress,
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
    call_rows = assembly_call_rows(
        loci,
        products,
        sample_id,
        read_support,
        min_posterior,
        algorithm,
        assembly_round_tolerance,
    )
    write_tsv(call_rows, calls_path, SIMPLE_CALL_FIELDS)
    write_tsv(call_rows, outdir_path / "locus_repeat_counts.tsv", REPEAT_COUNT_FIELDS)
    allele_distribution_path = outdir_path / "allele_probability_distribution.tsv"
    write_tsv(
        allele_distribution_rows(sample_id, call_rows),
        allele_distribution_path,
        ALLELE_DISTRIBUTION_FIELDS,
    )
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
    profile_matches_path = outdir_path / "profile_matches.tsv"
    profile_match_loci_path = outdir_path / "profile_match_loci.tsv"
    write_tsv(
        profile_match_locus_rows(
            sample_id,
            fingerprint_rows[0],
            profiles,
            match_rows,
            allele_rows,
        ),
        profile_match_loci_path,
        PROFILE_MATCH_LOCUS_FIELDS,
    )
    # Re-run product selection without depth so the comparison files retain
    # historical mismatch-round behavior even when modern calls use support.
    legacy_call_rows = legacy_assembly_call_rows(
        loci,
        products,
        sample_id,
        round_tolerance=assembly_round_tolerance,
    )
    legacy_paths = write_legacy_assembly_outputs(
        outdir_path,
        loci,
        products,
        legacy_call_rows,
        sample_id,
        assembly_round_tolerance,
    )
    phylogeny_paths: dict[str, Path] = {}
    phylogenetic_rows: list[dict] = []
    closest_reference_bands: list[dict] = []
    if database_path:
        product_by_id = {product["product_id"]: product for product in products}
        query_sequences = {
            row["locus_id"]: product_by_id[row["evidence"]]["sequence"]
            for row in call_rows
            if row.get("evidence") in product_by_id
        }
        progress.step(
            "Placing assembly loci with EPA-ng using reusable reference trees when available"
        )
        phylogeny_paths = run_phylogenetic_placement(
            query_sequences,
            database_path,
            outdir_path,
            sample_id,
            loci,
            thread_count,
            mafft_bin=mafft_bin,
            raxml_ng_bin=raxml_ng_bin,
            epa_ng_bin=epa_ng_bin,
            raxml_model=raxml_model,
            snp_weight=phylogeny_snp_weight,
            repeat_weight=phylogeny_repeat_weight,
            reference_metadata_path=reference_metadata_path,
            progress=progress,
            query_assembly_path=assembly_path,
            dnadiff_bin=dnadiff_bin,
        )
        phylogenetic_rows = read_profiles(phylogeny_paths["combined_marker_matches"])
        closest_reference_bands = read_profiles(
            phylogeny_paths["closest_reference_bands"]
        )
    output_match_rows = match_rows + sequence_reference_match_rows(
        phylogenetic_rows
    )
    write_tsv(output_match_rows, profile_matches_path, MATCH_FIELDS)
    progress.step("Writing HTML report")
    write_assembly_report(
        outdir_path,
        sample_id,
        call_rows,
        products,
        match_rows,
        profiles,
        loci,
        phylogenetic_rows,
        closest_reference_bands,
    )
    progress.step(f"Done. Main calls: {calls_path}")
    return {
        "outdir": outdir_path,
        "calls": calls_path,
        "repeat_counts": outdir_path / "locus_repeat_counts.tsv",
        "allele_distribution": allele_distribution_path,
        "amplicons": outdir_path / "assembly_amplicons.tsv",
        "amplicon_fasta": product_fasta,
        "in_silico_pcr": outdir_path / "in_silico_pcr",
        "read_support": support_path,
        "fingerprint": outdir_path / "mlva_fingerprint.tsv",
        "profile_matches": profile_matches_path,
        "profile_match_loci": profile_match_loci_path,
        "report": outdir_path / "report.html",
        **legacy_paths,
        **phylogeny_paths,
    }
