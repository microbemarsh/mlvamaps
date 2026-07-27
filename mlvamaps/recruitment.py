from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from .calling import (
    estimate_repeat_count_from_product_length,
    legacy_round_repeat_count,
    normalize_allele,
)
from .io import read_fasta, write_fasta, write_fastq, write_tsv
from .mapping import (
    build_minimap2_map_command,
    check_minimap2,
    run_minimap2_command,
)
from .models import Assignment, Locus, ReadPrediction, ReadRecord
from .sequence import majority_consensus, revcomp


RECRUITMENT_READ_FIELDS = [
    "read_id",
    "locus_id",
    "reference_name",
    "reference_source",
    "candidate_allele",
    "mapping_quality",
    "alignment_identity",
    "locus_score_margin",
    "aligned_query_bp",
    "reference_coverage",
    "full_product",
    "genotype_informative",
    "evidence_class",
]

RECRUITMENT_SUMMARY_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_source",
    "mapped_reads",
    "full_product_reads",
    "genotype_informative_reads",
    "candidate_alleles",
    "presence_status",
]

_SAFE_NAME = re.compile(r"[^A-Za-z0-9_.:-]+")


def _database_root(database_path: str | Path | None, loci: list[Locus]) -> Path | None:
    if not database_path:
        return None
    path = Path(database_path)
    candidates = [path, path / "database"]
    for candidate in candidates:
        if any((candidate / f"{locus.locus_id}.fasta").exists() for locus in loci):
            return candidate
    return None


def _canonical_product(
    locus: Locus,
    database_root: Path | None,
) -> tuple[str, str]:
    if database_root is not None:
        path = database_root / f"{locus.locus_id}.fasta"
        if path.exists():
            records = list(read_fasta(path))
            if records:
                # Prefer the modal product length so one unusual reference does
                # not define recruitment for the whole panel.
                lengths = Counter(len(sequence) for _name, sequence in records)
                modal_length = min(
                    lengths,
                    key=lambda length: (-lengths[length], length),
                )
                sequence = next(
                    sequence
                    for _name, sequence in records
                    if len(sequence) == modal_length
                )
                return sequence, "database_product"

    repeat_length = locus.repeat_unit_length_bp or len(locus.repeat_motif)
    if (
        locus.forward_primer
        and locus.reverse_primer
        and locus.repeat_motif
        and repeat_length
    ):
        nominal = locus.nominal_repeat_units
        if not nominal:
            nominal = max(
                locus.expected_min_repeats,
                min(
                    locus.expected_max_repeats,
                    round(
                        (locus.expected_min_repeats + locus.expected_max_repeats)
                        / 2
                    ),
                ),
            )
        repeat_bp = round(float(nominal) * repeat_length)
        repeated = (locus.repeat_motif * math.ceil(repeat_bp / len(locus.repeat_motif)))[
            :repeat_bp
        ]
        return (
            locus.forward_primer
            + locus.left_flank_sequence
            + repeated
            + locus.right_flank_sequence
            + revcomp(locus.reverse_primer),
            "synthetic_panel_product",
        )
    return "", "unavailable"


def _repeat_bounds(locus: Locus, sequence: str) -> tuple[int, int] | None:
    if not locus.left_flank_sequence or not locus.right_flank_sequence:
        return None
    left = sequence.find(locus.left_flank_sequence)
    if left < 0:
        return None
    start = left + len(locus.left_flank_sequence)
    end = sequence.find(locus.right_flank_sequence, start)
    if end < start:
        return None
    return start, end


def build_recruitment_references(
    loci: list[Locus],
    database_path: str | Path | None = None,
) -> list[dict]:
    """Build competitive locus/allele products for long-read recruitment."""
    database_root = _database_root(database_path, loci)
    references = []
    reference_index = 0
    for locus in loci:
        canonical, source = _canonical_product(locus, database_root)
        if not canonical:
            continue
        bounds = _repeat_bounds(locus, canonical)
        products: list[tuple[int | float | str, str, tuple[int, int] | None]] = []
        if bounds and locus.repeat_motif:
            start, end = bounds
            repeat_length = locus.repeat_unit_length_bp or len(locus.repeat_motif)
            minimum_bp = max(
                0, math.floor(float(locus.expected_min_repeats) * repeat_length)
            )
            maximum_bp = max(
                minimum_bp,
                math.ceil(float(locus.expected_max_repeats) * repeat_length),
            )
            if maximum_bp - minimum_bp <= 500:
                repeat_sizes = range(minimum_bp, maximum_bp + 1)
            else:
                # Very broad default panels would otherwise create thousands
                # of nearly identical mapping targets. Retain motif-resolution
                # coverage until the panel supplies tighter bounds.
                repeat_sizes = range(
                    minimum_bp,
                    maximum_bp + 1,
                    max(repeat_length, 1),
                )
            for repeat_bp in repeat_sizes:
                repeated = (
                    locus.repeat_motif
                    * math.ceil(max(repeat_bp, 1) / len(locus.repeat_motif))
                )[:repeat_bp]
                product = canonical[:start] + repeated + canonical[end:]
                raw = estimate_repeat_count_from_product_length(
                    locus, len(product)
                )
                if raw is None:
                    raw = repeat_bp / repeat_length
                allele = legacy_round_repeat_count(raw)
                products.append(
                    (
                        allele,
                        product,
                        (start, start + len(repeated)),
                    )
                )
        else:
            raw = estimate_repeat_count_from_product_length(locus, len(canonical))
            products.append(("" if raw is None else raw, canonical, bounds))

        seen_sequences: set[str] = set()
        for allele, product, product_bounds in products:
            if product in seen_sequences:
                continue
            seen_sequences.add(product)
            safe_locus = _SAFE_NAME.sub("_", locus.locus_id)
            reference_name = f"recruit_{reference_index:06d}_{safe_locus}"
            reference_index += 1
            references.append(
                {
                    "reference_name": reference_name,
                    "locus_id": locus.locus_id,
                    "reference_source": source,
                    "candidate_allele": allele,
                    "sequence": product,
                    "repeat_start": (
                        "" if product_bounds is None else product_bounds[0]
                    ),
                    "repeat_end": (
                        "" if product_bounds is None else product_bounds[1]
                    ),
                }
            )
    return references


def _quality_string(values: list[int] | None, length: int) -> str:
    if values is None or len(values) != length:
        return "I" * length
    return "".join(chr(min(max(value, 0), 93) + 33) for value in values)


def parse_recruitment_sam(
    sam_path: str | Path,
    references: list[dict],
    reads: list[ReadRecord],
    loci: list[Locus],
    min_mapping_quality: int = 0,
    min_alignment_identity: float = 0.9,
    min_aligned_bp: int = 100,
    min_locus_score_margin: int = 10,
    full_product_edge_tolerance: int = 0,
) -> tuple[list[dict], list[Assignment]]:
    """Parse competitive mappings into presence evidence and full products."""
    if min_mapping_quality < 0 or min_aligned_bp < 1:
        raise ValueError("recruitment mapping thresholds are invalid")
    if not 0 <= min_alignment_identity <= 1:
        raise ValueError("minimum recruitment identity must be between 0 and 1")
    if min_locus_score_margin < 0:
        raise ValueError("minimum locus score margin must be non-negative")
    reference_by_name = {row["reference_name"]: row for row in references}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    rows = []
    assignments = []
    candidates_by_read: dict[str, list[dict]] = defaultdict(list)
    with pysam.AlignmentFile(str(sam_path), "r", check_sq=False) as alignments:
        for alignment in alignments.fetch(until_eof=True):
            if (
                alignment.is_unmapped
                or alignment.is_supplementary
            ):
                continue
            reference = reference_by_name.get(alignment.reference_name)
            if reference is None:
                continue
            aligned_bp = int(alignment.query_alignment_length or 0)
            try:
                edits = int(alignment.get_tag("NM"))
            except KeyError:
                edits = 0
            identity = max(0.0, 1 - edits / max(aligned_bp, 1))
            reference_length = len(reference["sequence"])
            required_bp = min(min_aligned_bp, max(1, reference_length // 2))
            if (
                alignment.mapping_quality < min_mapping_quality
                or identity < min_alignment_identity
                or aligned_bp < required_bp
            ):
                continue
            try:
                alignment_score = int(alignment.get_tag("AS"))
            except KeyError:
                alignment_score = aligned_bp - edits
            candidates_by_read[alignment.query_name].append(
                {
                    "alignment": alignment,
                    "reference": reference,
                    "aligned_bp": aligned_bp,
                    "identity": identity,
                    "score": alignment_score,
                }
            )

    for query_name, candidates in candidates_by_read.items():
        best_by_locus: dict[str, dict] = {}
        for candidate in candidates:
            locus_id = str(candidate["reference"]["locus_id"])
            current = best_by_locus.get(locus_id)
            key = (
                int(candidate["score"]),
                int(not candidate["alignment"].is_secondary),
                int(candidate["aligned_bp"]),
            )
            if current is None or key > (
                int(current["score"]),
                int(not current["alignment"].is_secondary),
                int(current["aligned_bp"]),
            ):
                best_by_locus[locus_id] = candidate
        ranked_loci = sorted(
            best_by_locus.values(),
            key=lambda candidate: (
                -int(candidate["score"]),
                str(candidate["reference"]["locus_id"]),
            ),
        )
        selected = ranked_loci[0]
        next_score = (
            int(ranked_loci[1]["score"]) if len(ranked_loci) > 1 else None
        )
        locus_margin = (
            int(selected["score"]) - next_score
            if next_score is not None
            else ""
        )
        if next_score is not None and int(locus_margin) < min_locus_score_margin:
            continue
        alignment = selected["alignment"]
        reference = selected["reference"]
        aligned_bp = int(selected["aligned_bp"])
        identity = float(selected["identity"])
        reference_length = len(reference["sequence"])
        reference_start = int(alignment.reference_start)
        reference_end = int(alignment.reference_end or 0)
        repeat_start = reference.get("repeat_start")
        repeat_end = reference.get("repeat_end")
        aligned_reference_positions = {
            int(reference_position)
            for query_position, reference_position in alignment.get_aligned_pairs(
                matches_only=False
            )
            if query_position is not None and reference_position is not None
        }

        def anchor_supported(window_start: int, window_end: int) -> bool:
            window_length = max(0, window_end - window_start)
            required = min(4, window_length)
            return (
                sum(
                    position in aligned_reference_positions
                    for position in range(window_start, window_end)
                )
                >= required
            )

        left_anchor = (
            repeat_start not in ("", None)
            and anchor_supported(max(0, int(repeat_start) - 8), int(repeat_start))
        )
        right_anchor = (
            repeat_end not in ("", None)
            and anchor_supported(
                int(repeat_end),
                min(reference_length, int(repeat_end) + 8),
            )
        )
        genotype_informative = (
            repeat_start not in ("", None)
            and repeat_end not in ("", None)
            and reference_start <= int(repeat_start)
            and reference_end >= int(repeat_end)
            and left_anchor
            and right_anchor
        )
        full_product = (
            reference_start <= full_product_edge_tolerance
            and reference_end >= reference_length - full_product_edge_tolerance
        )
        evidence_class = (
            "FULL_PRODUCT"
            if full_product
            else "REPEAT_INFORMATIVE"
            if genotype_informative
            else "PRESENCE_ONLY"
        )
        rows.append(
            {
                "read_id": query_name,
                "locus_id": reference["locus_id"],
                "reference_name": reference["reference_name"],
                "reference_source": reference["reference_source"],
                "candidate_allele": reference["candidate_allele"],
                "mapping_quality": alignment.mapping_quality,
                "alignment_identity": round(identity, 6),
                "locus_score_margin": locus_margin,
                "aligned_query_bp": aligned_bp,
                "reference_coverage": round(
                    (reference_end - reference_start)
                    / max(reference_length, 1),
                    6,
                ),
                "full_product": "yes" if full_product else "no",
                "genotype_informative": "yes" if genotype_informative else "no",
                "evidence_class": evidence_class,
            }
        )
        if not full_product:
            continue
        locus = locus_by_id[reference["locus_id"]]
        product_sequence = (alignment.query_alignment_sequence or "").upper()
        if not product_sequence:
            continue
        quality_values = (
            list(alignment.query_alignment_qualities)
            if alignment.query_alignment_qualities is not None
            else None
        )
        quality = _quality_string(quality_values, len(product_sequence))
        assignments.append(
            Assignment(
                read_id=query_name,
                sample_id="",
                assigned_locus=locus.locus_id,
                assignment_score=round(identity, 4),
                orientation="reverse" if alignment.is_reverse else "forward",
                primer_forward_detected=True,
                primer_reverse_detected=True,
                passes_assignment_qc=True,
                oriented_sequence=product_sequence,
                oriented_quality=quality,
                forward_start=0,
                forward_end=min(len(locus.forward_primer), len(product_sequence)),
                reverse_start=max(
                    0, len(product_sequence) - len(locus.reverse_primer)
                ),
                reverse_end=len(product_sequence),
                forward_mismatches=0,
                reverse_mismatches=0,
                product_size_bp=len(product_sequence),
            )
        )
    return rows, assignments


def recruitment_summary_rows(
    sample_id: str,
    loci: list[Locus],
    rows: list[dict],
) -> list[dict]:
    by_locus: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_locus[str(row["locus_id"])].append(row)
    output = []
    for locus in loci:
        locus_rows = by_locus.get(locus.locus_id, [])
        full = sum(row["full_product"] == "yes" for row in locus_rows)
        informative = sum(
            row["genotype_informative"] == "yes" for row in locus_rows
        )
        if full > 1:
            status = "PRESENT_GENOTYPED"
        elif full == 1 or informative:
            status = "PRESENT_PROVISIONAL"
        elif locus_rows:
            status = "PRESENT_UNTYPED"
        else:
            status = "NO_EVIDENCE"
        sources = sorted({str(row["reference_source"]) for row in locus_rows})
        alleles = sorted(
            {
                str(row["candidate_allele"])
                for row in locus_rows
                if row["genotype_informative"] == "yes"
            }
        )
        output.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "reference_source": ",".join(sources),
                "mapped_reads": len(locus_rows),
                "full_product_reads": full,
                "genotype_informative_reads": informative,
                "candidate_alleles": ",".join(alleles),
                "presence_status": status,
            }
        )
    return output


def recruitment_fallback_evidence(
    rows: list[dict],
    loci: list[Locus],
    feature_loci: set[str],
    sample_id: str,
) -> tuple[list[dict], list[ReadPrediction]]:
    """Create provisional allele evidence when no complete product was found."""
    locus_by_id = {locus.locus_id: locus for locus in loci}
    grouped: dict[tuple[str, int | float], list[dict]] = defaultdict(list)
    for row in rows:
        locus_id = str(row.get("locus_id", ""))
        if locus_id in feature_loci or row.get("genotype_informative") != "yes":
            continue
        try:
            allele = normalize_allele(float(row.get("candidate_allele")))
        except (TypeError, ValueError):
            continue
        grouped[(locus_id, allele)].append(row)

    totals = Counter(
        locus_id
        for (locus_id, _allele), locus_rows in grouped.items()
        for _row in locus_rows
    )
    asv_rows = []
    predictions = []
    rank_by_locus: Counter[str] = Counter()
    for (locus_id, allele), locus_rows in sorted(
        grouped.items(),
        key=lambda item: (
            item[0][0],
            -len(item[1]),
            float(item[0][1]),
        ),
    ):
        rank_by_locus[locus_id] += 1
        variant_id = f"{locus_id}_RECRUIT{rank_by_locus[locus_id]}"
        locus = locus_by_id[locus_id]
        motif = locus.repeat_motif or "N"
        repeat_length = locus.repeat_unit_length_bp or len(motif)
        repeat_bp = round(float(allele) * repeat_length)
        repeat_sequence = (motif * math.ceil(max(repeat_bp, 1) / len(motif)))[
            :repeat_bp
        ]
        support = len(locus_rows)
        asv_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "variant_id": variant_id,
                "repeat_count": allele,
                "support_reads": support,
                "unique_sequences": 0,
                "frequency": round(support / max(totals[locus_id], 1), 6),
                "representative_read_id": str(locus_rows[0]["read_id"]),
                "representative_pattern": "recruitment_candidate",
                "representative_sequence": repeat_sequence,
                "representative_length_bp": len(repeat_sequence),
                "reads_with_indels": 0,
                "total_insertions": 0,
                "total_deletions": 0,
                "total_substitutions": 0,
                "mean_edit_distance_to_representative": 0,
                "max_edit_distance_to_representative": 0,
            }
        )
        for row in locus_rows:
            identity = float(row.get("alignment_identity") or 0)
            predictions.append(
                ReadPrediction(
                    read_id=str(row["read_id"]),
                    locus_id=locus_id,
                    predicted_repeat_count=allele,
                    probability=round(max(0.0, min(identity, 1.0)), 6),
                    top_alt_repeat_count=None,
                    top_alt_probability=0.0,
                    variant_id=variant_id,
                    insertions_vs_representative=0,
                    deletions_vs_representative=0,
                    substitutions_vs_representative=0,
                    evidence_weight=round(max(0.0, min(identity, 1.0)), 6),
                    raw_repeat_count_estimate=float(allele),
                    measurement_sigma=0.2,
                    measurement_repeat_count_estimate=float(allele),
                )
            )
    return asv_rows, predictions


def local_product_records(
    assignments: list[Assignment],
) -> list[tuple[str, str]]:
    by_locus: dict[str, list[Assignment]] = defaultdict(list)
    for assignment in assignments:
        by_locus[assignment.assigned_locus].append(assignment)
    products = []
    for locus_id, locus_assignments in sorted(by_locus.items()):
        product_sequences = [
            assignment.oriented_sequence[
                assignment.forward_start or 0 : (
                    assignment.reverse_end
                    if assignment.reverse_end is not None
                    else len(assignment.oriented_sequence)
                )
            ]
            for assignment in locus_assignments
        ]
        lengths = Counter(
            len(sequence) for sequence in product_sequences
        )
        modal_length = min(lengths, key=lambda length: (-lengths[length], length))
        sequences = [
            sequence
            for sequence in product_sequences
            if len(sequence) == modal_length
        ]
        products.append((f"{locus_id}_local_primary", majority_consensus(sequences)))
    return products


def run_read_recruitment(
    reads: list[ReadRecord],
    loci: list[Locus],
    outdir: str | Path,
    sample_id: str,
    database_path: str | Path | None,
    threads: int,
    executable: str = "minimap2",
    preset: str | None = None,
    min_mapping_quality: int = 0,
    min_alignment_identity: float = 0.9,
    min_aligned_bp: int = 100,
    min_locus_score_margin: int = 10,
) -> tuple[list[dict], list[dict], list[Assignment], dict[str, Path]]:
    outdir = Path(outdir)
    root = outdir / "recruitment"
    root.mkdir(parents=True, exist_ok=True)
    references_path = root / "locus_recruitment_references.fasta"
    reads_path = root / "filtered_reads.fastq"
    sam_path = root / "read_recruitment.sam"
    read_rows_path = outdir / "locus_recruited_reads.tsv"
    summary_path = outdir / "locus_presence.tsv"
    local_products_path = outdir / "local_locus_products.fasta"
    references = build_recruitment_references(loci, database_path)
    write_fasta(
        (
            (
                f"{row['reference_name']} locus={row['locus_id']} "
                f"allele={row['candidate_allele']} source={row['reference_source']}",
                row["sequence"],
            )
            for row in references
        ),
        references_path,
    )
    write_fastq(reads, reads_path)
    read_rows: list[dict] = []
    assignments: list[Assignment] = []
    if references and reads:
        executable_path = check_minimap2(executable)
        run_minimap2_command(
            build_minimap2_map_command(
                references_path,
                reads_path,
                threads,
                executable=executable_path,
                preset=preset,
            ),
            "competitive locus recruitment",
            stdout_path=sam_path,
        )
        read_rows, assignments = parse_recruitment_sam(
            sam_path,
            references,
            reads,
            loci,
            min_mapping_quality=min_mapping_quality,
            min_alignment_identity=min_alignment_identity,
            min_aligned_bp=min_aligned_bp,
            min_locus_score_margin=min_locus_score_margin,
        )
    else:
        sam_path.write_text("")
    for assignment_index, assignment in enumerate(assignments):
        assignments[assignment_index] = Assignment(
            **{**assignment.__dict__, "sample_id": sample_id}
        )
    summaries = recruitment_summary_rows(sample_id, loci, read_rows)
    write_tsv(read_rows, read_rows_path, RECRUITMENT_READ_FIELDS)
    write_tsv(summaries, summary_path, RECRUITMENT_SUMMARY_FIELDS)
    write_fasta(local_product_records(assignments), local_products_path)
    return read_rows, summaries, assignments, {
        "recruitment_references": references_path,
        "recruitment_alignments": sam_path,
        "recruited_reads": read_rows_path,
        "locus_presence": summary_path,
        "local_products": local_products_path,
    }
