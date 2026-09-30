from __future__ import annotations

import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import pysam

from .models import RepeatFeature


_SAM_NAME = re.compile(r"[^A-Za-z0-9_.:-]+")
_DNA_BASES = "ACGT"


MAPPING_SUMMARY_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_variant_id",
    "reference_read_id",
    "reference_length_bp",
    "total_reads",
    "mapped_reads",
    "mapping_rate",
    "mean_mapping_quality",
    "mean_depth",
    "covered_bases",
    "coverage_percent",
    "snp_count",
]

SNP_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_variant_id",
    "position",
    "reference_base",
    "alternate_base",
    "depth",
    "alternate_depth",
    "alternate_frequency",
    "mean_alternate_base_quality",
]


def check_minimap2(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"minimap2 executable {executable!r} was not found. Install minimap2 "
            "from Bioconda, pass --minimap2-bin, or disable the requested "
            "read-mapping stage."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    detail = f"{result.stdout}\n{result.stderr}".lower()
    if result.returncode or not detail.strip():
        raise RuntimeError(f"Could not run minimap2 at {path}")
    return path


def build_minimap2_map_command(
    reference_path: str | Path,
    reads_path: str | Path,
    threads: int,
    executable: str = "minimap2",
    preset: str | None = None,
) -> list[str]:
    command = [executable, "-a", "-t", str(threads)]
    if preset:
        command.extend(["-x", preset])
    command.extend([str(reference_path), str(reads_path)])
    return command


def run_minimap2_command(
    command: list[str], stage: str, stdout_path: Path | None = None
) -> None:
    if stdout_path is None:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    else:
        with stdout_path.open("w") as output:
            result = subprocess.run(
                command,
                stdout=output,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(
            f"minimap2 {stage} failed (exit {result.returncode}): {detail}"
        )


def _safe_name(value: str, fallback: str) -> str:
    cleaned = _SAM_NAME.sub("_", value).strip("_")
    return cleaned or fallback


def _dominant_references(
    features: list[RepeatFeature], asv_rows: list[dict]
) -> list[dict]:
    feature_by_key = {(feature.locus_id, feature.read_id): feature for feature in features}
    dominant_by_locus: dict[str, dict] = {}
    for row in asv_rows:
        locus_id = str(row["locus_id"])
        current = dominant_by_locus.get(locus_id)
        if current is None:
            dominant_by_locus[locus_id] = row
            continue
        if int(row.get("support_reads") or 0) > int(
            current.get("support_reads") or 0
        ):
            dominant_by_locus[locus_id] = row

    references = []
    used_names: set[str] = set()
    for index, (locus_id, row) in enumerate(sorted(dominant_by_locus.items())):
        read_id = str(row.get("representative_read_id") or "")
        feature = feature_by_key.get((locus_id, read_id))
        if feature is None:
            continue
        sequence = (feature.amplicon_sequence or feature.repeat_sequence).upper()
        if not sequence:
            continue
        name = _safe_name(str(row.get("variant_id") or locus_id), f"locus_{index:04d}")
        if name in used_names:
            name = f"{name}_{index:04d}"
        used_names.add(name)
        references.append(
            {
                "reference_name": name,
                "locus_id": locus_id,
                "reference_variant_id": str(row.get("variant_id") or ""),
                "reference_read_id": read_id,
                "sequence": sequence,
            }
        )
    return references


def _write_mapping_inputs(
    features: list[RepeatFeature],
    references: list[dict],
    reference_path: Path,
    reads_path: Path,
) -> tuple[dict[str, dict], dict[str, dict]]:
    reference_by_name = {row["reference_name"]: row for row in references}
    locus_to_reference = {row["locus_id"]: row for row in references}
    with reference_path.open("w") as handle:
        for row in references:
            handle.write(
                f">{row['reference_name']} locus={row['locus_id']} "
                f"representative_read={row['reference_read_id']}\n{row['sequence']}\n"
            )

    query_metadata: dict[str, dict] = {}
    with reads_path.open("w") as handle:
        query_index = 0
        for feature in features:
            reference = locus_to_reference.get(feature.locus_id)
            if reference is None:
                continue
            sequence = (feature.amplicon_sequence or feature.repeat_sequence).upper()
            if not sequence:
                continue
            quality = feature.amplicon_quality or ("I" * len(sequence))
            if len(quality) != len(sequence):
                quality = "I" * len(sequence)
            query_name = f"mlvamaps_read_{query_index}"
            query_index += 1
            query_metadata[query_name] = {
                "read_id": feature.read_id,
                "locus_id": feature.locus_id,
                "expected_reference": reference["reference_name"],
            }
            handle.write(f"@{query_name}\n{sequence}\n+\n{quality}\n")
    return reference_by_name, query_metadata


def parse_minimap2_sam(
    sam_path: str | Path,
    references: dict[str, dict],
    queries: dict[str, dict],
    sample_id: str,
    min_mapping_quality: int = 0,
    min_base_quality: int = 20,
    min_snp_depth: int = 3,
    min_snp_alternate_reads: int = 2,
    min_snp_frequency: float = 0.2,
) -> tuple[list[dict], list[dict]]:
    """Summarize locus-relative mappings and call simple quality-filtered SNPs."""
    if min_mapping_quality < 0 or min_base_quality < 0:
        raise ValueError("mapping and base quality thresholds must be non-negative")
    if min_snp_depth < 1 or min_snp_alternate_reads < 1:
        raise ValueError("SNP depth and alternate-read thresholds must be at least 1")
    if not 0.0 <= min_snp_frequency <= 1.0:
        raise ValueError("min_snp_frequency must be between 0 and 1")
    total_reads = Counter(row["locus_id"] for row in queries.values())
    mapped_reads: Counter[str] = Counter()
    mapping_qualities: dict[str, list[int]] = defaultdict(list)
    base_counts: dict[str, dict[int, Counter[str]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    base_quality_sums: dict[str, dict[int, Counter[str]]] = defaultdict(
        lambda: defaultdict(Counter)
    )
    seen_primary: set[str] = set()

    with pysam.AlignmentFile(str(sam_path), "r", check_sq=False) as alignments:
        for alignment in alignments.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.query_name in seen_primary:
                continue
            query = queries.get(alignment.query_name)
            if query is None or alignment.reference_name != query["expected_reference"]:
                continue
            if alignment.mapping_quality < min_mapping_quality:
                continue
            reference = references.get(alignment.reference_name)
            if reference is None:
                continue
            seen_primary.add(alignment.query_name)
            locus_id = reference["locus_id"]
            mapped_reads[locus_id] += 1
            mapping_qualities[locus_id].append(alignment.mapping_quality)
            sequence = alignment.query_sequence or ""
            qualities = alignment.query_qualities
            for query_position, reference_position in alignment.get_aligned_pairs(
                matches_only=False
            ):
                if query_position is None or reference_position is None:
                    continue
                if query_position >= len(sequence):
                    continue
                quality = qualities[query_position] if qualities is not None else 40
                if quality < min_base_quality:
                    continue
                base = sequence[query_position].upper()
                if base not in _DNA_BASES:
                    continue
                base_counts[locus_id][reference_position][base] += 1
                base_quality_sums[locus_id][reference_position][base] += quality

    snp_rows: list[dict] = []
    for reference in references.values():
        locus_id = reference["locus_id"]
        sequence = reference["sequence"]
        for position, counts in sorted(base_counts[locus_id].items()):
            if position < 0 or position >= len(sequence):
                continue
            depth = sum(counts.values())
            if depth < min_snp_depth:
                continue
            reference_base = sequence[position]
            if reference_base not in _DNA_BASES:
                continue
            for alternate_base in _DNA_BASES:
                if alternate_base == reference_base:
                    continue
                alternate_depth = counts[alternate_base]
                frequency = alternate_depth / depth
                if (
                    alternate_depth < min_snp_alternate_reads
                    or frequency < min_snp_frequency
                ):
                    continue
                snp_rows.append(
                    {
                        "sample_id": sample_id,
                        "locus_id": locus_id,
                        "reference_variant_id": reference["reference_variant_id"],
                        "position": position + 1,
                        "reference_base": reference_base,
                        "alternate_base": alternate_base,
                        "depth": depth,
                        "alternate_depth": alternate_depth,
                        "alternate_frequency": round(frequency, 6),
                        "mean_alternate_base_quality": round(
                            base_quality_sums[locus_id][position][alternate_base]
                            / alternate_depth,
                            2,
                        ),
                    }
                )

    snp_counts = Counter(row["locus_id"] for row in snp_rows)
    summary_rows = []
    for reference in sorted(references.values(), key=lambda row: row["locus_id"]):
        locus_id = reference["locus_id"]
        sequence_length = len(reference["sequence"])
        depths = [
            sum(base_counts[locus_id][position].values())
            for position in range(sequence_length)
        ]
        covered_bases = sum(depth > 0 for depth in depths)
        mapped = mapped_reads[locus_id]
        total = total_reads[locus_id]
        mapqs = mapping_qualities[locus_id]
        summary_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "reference_variant_id": reference["reference_variant_id"],
                "reference_read_id": reference["reference_read_id"],
                "reference_length_bp": sequence_length,
                "total_reads": total,
                "mapped_reads": mapped,
                "mapping_rate": round(mapped / total, 6) if total else 0,
                "mean_mapping_quality": round(sum(mapqs) / len(mapqs), 3)
                if mapqs
                else 0,
                "mean_depth": round(sum(depths) / sequence_length, 3)
                if sequence_length
                else 0,
                "covered_bases": covered_bases,
                "coverage_percent": round(100 * covered_bases / sequence_length, 3)
                if sequence_length
                else 0,
                "snp_count": snp_counts[locus_id],
            }
        )
    return summary_rows, snp_rows


def run_locus_mapping(
    features: list[RepeatFeature],
    asv_rows: list[dict],
    outdir: str | Path,
    sample_id: str,
    threads: int,
    executable: str = "minimap2",
    preset: str | None = "sr",
    min_mapping_quality: int = 0,
    min_base_quality: int = 20,
    min_snp_depth: int = 3,
    min_snp_alternate_reads: int = 2,
    min_snp_frequency: float = 0.2,
    primary_product_sequences: dict[str, str] | None = None,
) -> tuple[list[dict], list[dict], dict[str, Path]]:
    """Map locus reads to POA products and report support/SNP evidence."""
    outdir = Path(outdir)
    minimap2_root = outdir / "minimap2"
    minimap2_root.mkdir(parents=True, exist_ok=True)
    work_dir = minimap2_root / "locus_mapping"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    if primary_product_sequences:
        references = [
            {
                "reference_name": _safe_name(
                    f"{locus_id}_POA", f"locus_{index:04d}"
                ),
                "locus_id": locus_id,
                "reference_variant_id": f"{locus_id}_POA",
                "reference_read_id": "",
                "sequence": sequence.upper(),
            }
            for index, (locus_id, sequence) in enumerate(
                sorted(primary_product_sequences.items())
            )
            if sequence
        ]
    else:
        references = _dominant_references(features, asv_rows)
    public_reference = outdir / "locus_mapping_references.fasta.gz"
    internal_reference = work_dir / "locus_mapping_references.fasta.gz"
    reads_path = work_dir / "locus_reads.fastq.gz"
    sam_path = outdir / "locus_read_alignments.sam"

    reference_by_name, query_metadata = _write_mapping_inputs(
        features, references, internal_reference, reads_path
    )
    shutil.copyfile(internal_reference, public_reference)
    paths = {
        "mapping_references": public_reference,
        "mapping_alignments": sam_path,
        "minimap2": minimap2_root,
    }
    if not references or not query_metadata:
        sam_path.write_text("")
        return [], [], paths

    executable_path = check_minimap2(executable)
    run_minimap2_command(
        build_minimap2_map_command(
            internal_reference,
            reads_path,
            threads,
            executable_path,
            preset=preset,
        ),
        "mapping",
        stdout_path=sam_path,
    )
    summary_rows, snp_rows = parse_minimap2_sam(
        sam_path,
        reference_by_name,
        query_metadata,
        sample_id,
        min_mapping_quality=min_mapping_quality,
        min_base_quality=min_base_quality,
        min_snp_depth=min_snp_depth,
        min_snp_alternate_reads=min_snp_alternate_reads,
        min_snp_frequency=min_snp_frequency,
    )
    return summary_rows, snp_rows, paths
