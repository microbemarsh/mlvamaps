from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import shutil
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .calling import estimate_repeat_count_from_product_length, normalize_allele, repeat_unit_length
from .concurrency import resolve_threads
from .io import gzip_output_file
from .models import Locus, RepeatFeature
from .progress import ProgressReporter
from .taxon_assignment import (
    TAXONOMIC_EVIDENCE_FIELDS,
    TAXONOMIC_SUMMARY_FIELDS,
    TAXON_ASSIGNMENT_FIELDS,
    TAXON_CANDIDATE_FIELDS,
    TAXON_LOCUS_FIELDS,
    TaxonCalibration,
    _finite,
    _reference_taxa,
    assign_best_taxon,
    assign_target_taxon,
)


PLACEMENT_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_id",
    "phylogenetic_distance",
    "likelihood_weighted_phylogenetic_distance",
    "placement_edge",
    "like_weight_ratio",
    "pendant_length",
    "distal_length",
    "placement_count",
    "placement_entropy",
    "reference_tree_scale",
]

SUMMARY_FIELDS = [
    "sample_id",
    "reference_id",
    "total_phylogenetic_distance",
    "total_likelihood_weighted_distance",
    "compared_loci",
    "mean_phylogenetic_distance",
    "mean_likelihood_weighted_distance",
    "distance_gap_to_next",
    "relative_distance_gap_to_next",
    "rank",
]

STATUS_FIELDS = ["locus_id", "reference_sequences", "query_sequence", "status"]

MARKER_COMPONENT_FIELDS = [
    "sample_id",
    "locus_id",
    "record_type",
    "reference_id",
    "repeat_count_raw",
    "repeat_count",
    "repeat_sequence",
    "repeat_haplotype",
    "repeat_region_start",
    "repeat_region_end",
    "snp_sequence_length",
    "masking_method",
]

LOCUS_MARKER_DISTANCE_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_id",
    "query_repeat_count",
    "reference_repeat_count",
    "repeat_count_delta",
    "normalized_repeat_distance",
    "likelihood_weighted_snp_distance",
    "reference_tree_scale",
    "placement_normalized_snp_distance",
    "direct_snp_distance",
    "direct_snp_scale",
    "normalized_direct_snp_distance",
    "exact_snp_match",
    "normalized_snp_distance",
    "placement_entropy",
]

COMBINED_MARKER_FIELDS = [
    "sample_id",
    "reference_id",
    "total_likelihood_weighted_snp_distance",
    "total_placement_normalized_snp_distance",
    "total_normalized_direct_snp_distance",
    "total_normalized_snp_distance",
    "total_repeat_count_distance",
    "total_normalized_repeat_distance",
    "snp_weight",
    "repeat_weight",
    "combined_marker_distance",
    "compared_loci",
    "repeat_compared_loci",
    "exact_snp_loci",
    "exact_marker_loci",
    "match_status",
    "ranking_warning",
    "whole_genome_exact_match",
    "whole_genome_snps",
    "whole_genome_indel_bases",
    "whole_genome_align_fraction_ref",
    "whole_genome_align_fraction_query",
    "tie_break_method",
    "tie_break_status",
    "distance_gap_to_next",
    "relative_distance_gap_to_next",
    "collection_date",
    "latitude",
    "longitude",
    "location",
    "source",
    "rank",
]

CLOSEST_REFERENCE_BAND_FIELDS = [
    "reference_id",
    "locus_id",
    "product_size_bp",
    "repeat_count",
]

REFERENCE_SEQUENCE_INDEX_FIELDS = [
    "index_version",
    "panel_sha256",
    "database_signature",
    "locus_id",
    "reference_id",
    "amplicon_sha256",
    "snp_sha256",
    "marker_sha256",
    "product_size_bp",
    "snp_sequence_length",
    "repeat_count_raw",
    "repeat_count",
]

REFERENCE_HAPLOTYPE_FIELDS = [
    "locus_id",
    "haplotype_id",
    "reference_id",
    "snp_sha256",
]

_FASTA_SUFFIXES = (".fa", ".fas", ".fasta", ".fna", ".ffn")


def _is_fasta_path(path: Path) -> bool:
    name = path.name.lower()
    return any(
        name.endswith(suffix) or name.endswith(f"{suffix}.gz")
        for suffix in _FASTA_SUFFIXES
    )


def _fasta_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in _FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem

REFERENCE_ASSEMBLY_FIELDS = ["reference_id", "assembly_file", "assembly_sha256"]
DNADIFF_RESULT_FIELDS = [
    "reference_id",
    "reference_file",
    "query_file",
    "exact_genome_match",
    "snps",
    "indel_bases",
    "align_fraction_ref",
    "align_fraction_query",
    "report_file",
]
_SAFE_FILE = re.compile(r"[^A-Za-z0-9_.-]+")
_IUPAC = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"), "K": frozenset("GT"),
    "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}
_RC = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


@dataclass(frozen=True)
class MarkerComponents:
    oriented_sequence: str
    snp_sequence: str
    repeat_sequence: str
    repeat_count_raw: float | None
    repeat_count: int | float | None
    repeat_region_start: int | None
    repeat_region_end: int | None
    masking_method: str


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_RC)[::-1]


def _iupac_find(pattern: str, sequence: str, start: int = 0) -> int:
    pattern = pattern.upper()
    sequence = sequence.upper()
    if not pattern:
        return -1
    for index in range(start, len(sequence) - len(pattern) + 1):
        if all(base in _IUPAC.get(code, frozenset(code)) for code, base in zip(pattern, sequence[index:])):
            return index
    return -1


def _longest_motif_run(sequence: str, motif: str) -> tuple[int, int] | None:
    motif = motif.upper()
    if not motif or set(motif) == {"N"}:
        return None
    motif_length = len(motif)
    best: tuple[int, int] | None = None
    for offset in range(len(sequence) - motif_length + 1):
        end = offset
        while end + motif_length <= len(sequence):
            chunk = sequence[end : end + motif_length]
            mismatches = sum(
                base not in _IUPAC.get(code, frozenset(code))
                for code, base in zip(motif, chunk)
            )
            if mismatches > motif_length // 8:
                break
            end += motif_length
        if end > offset and (best is None or end - offset > best[1] - best[0]):
            best = (offset, end)
    return best


def decompose_marker_sequence(locus: Locus, sequence: str) -> MarkerComponents:
    """Separate explicit VNTR characters from sequence used for SNP placement."""
    sequence = sequence.upper().replace("-", "")
    reverse_primer_site = _reverse_complement(locus.reverse_primer)
    forward = _iupac_find(locus.forward_primer, sequence)
    reverse = _iupac_find(reverse_primer_site, sequence, max(forward, 0))
    if forward < 0 or reverse < 0 or reverse <= forward:
        reverse_sequence = _reverse_complement(sequence)
        reverse_forward = _iupac_find(locus.forward_primer, reverse_sequence)
        reverse_reverse = _iupac_find(
            reverse_primer_site, reverse_sequence, max(reverse_forward, 0)
        )
        if reverse_forward >= 0 and reverse_reverse > reverse_forward:
            sequence = reverse_sequence
            forward, reverse = reverse_forward, reverse_reverse

    inner_start = forward + len(locus.forward_primer) if forward >= 0 else 0
    inner_end = reverse if reverse > inner_start else len(sequence)
    repeat_start: int | None = None
    repeat_end: int | None = None
    method = "unmasked"
    if locus.left_flank_sequence and locus.right_flank_sequence:
        left = _iupac_find(locus.left_flank_sequence, sequence, inner_start)
        right_start = left + len(locus.left_flank_sequence) if left >= 0 else inner_start
        right = _iupac_find(locus.right_flank_sequence, sequence, right_start)
        if left >= 0 and right >= right_start and right <= inner_end:
            repeat_start = right_start
            repeat_end = right
            method = "flank_bounded"
    if repeat_start is None:
        motif_run = _longest_motif_run(sequence[inner_start:inner_end], locus.repeat_motif)
        if motif_run is not None:
            repeat_start = inner_start + motif_run[0]
            repeat_end = inner_start + motif_run[1]
            method = "motif_run"

    repeat_sequence = (
        sequence[repeat_start:repeat_end]
        if repeat_start is not None and repeat_end is not None
        else ""
    )
    unit_length = repeat_unit_length(locus)
    if repeat_sequence and unit_length:
        raw_count = len(repeat_sequence) / unit_length
    else:
        raw_count = estimate_repeat_count_from_product_length(locus, len(sequence))
    repeat_count = normalize_allele(raw_count) if raw_count is not None else None
    snp_sequence = (
        sequence[:repeat_start] + sequence[repeat_end:]
        if repeat_start is not None and repeat_end is not None
        else sequence
    )
    if not snp_sequence:
        snp_sequence = sequence
        method += "_empty_mask_fallback"
    return MarkerComponents(
        sequence,
        snp_sequence,
        repeat_sequence,
        raw_count,
        repeat_count,
        repeat_start,
        repeat_end,
        method,
    )


def _repeat_haplotype(sequence: str, unit_length: int) -> str:
    if not sequence or unit_length <= 0:
        return ""
    return "|".join(
        sequence[index : index + unit_length]
        for index in range(0, len(sequence), unit_length)
    )


def _marker_component_row(
    sample_id: str,
    locus: Locus,
    record_type: str,
    reference_id: str,
    components: MarkerComponents,
) -> dict:
    return {
        "sample_id": sample_id,
        "locus_id": locus.locus_id,
        "record_type": record_type,
        "reference_id": reference_id,
        "repeat_count_raw": ""
        if components.repeat_count_raw is None
        else f"{components.repeat_count_raw:.6f}",
        "repeat_count": "" if components.repeat_count is None else components.repeat_count,
        "repeat_sequence": components.repeat_sequence,
        "repeat_haplotype": _repeat_haplotype(
            components.repeat_sequence, repeat_unit_length(locus)
        ),
        "repeat_region_start": ""
        if components.repeat_region_start is None
        else components.repeat_region_start,
        "repeat_region_end": ""
        if components.repeat_region_end is None
        else components.repeat_region_end,
        "snp_sequence_length": len(components.snp_sequence),
        "masking_method": components.masking_method,
    }


def _read_fasta(path: str | Path) -> list[tuple[str, str]]:
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0]
                sequence = []
            elif name is None:
                raise ValueError(f"Sequence appeared before a FASTA header in {path}")
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def _write_tsv(rows: list[dict], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _aligned_snp_distance(query: str, reference: str) -> float:
    """Return an ambiguity-aware mismatch fraction from a shared alignment."""
    if len(query) != len(reference):
        raise ValueError("Aligned query and reference sequences have different lengths")
    mismatches = 0
    comparable = 0
    for query_base, reference_base in zip(query.upper(), reference.upper()):
        if query_base == "-" and reference_base == "-":
            continue
        comparable += 1
        if query_base == reference_base:
            continue
        if query_base == "-" or reference_base == "-":
            mismatches += 1
            continue
        query_states = _IUPAC.get(query_base)
        reference_states = _IUPAC.get(reference_base)
        if query_states is None or reference_states is None or query_states.isdisjoint(
            reference_states
        ):
            mismatches += 1
    return mismatches / comparable if comparable else 0.0


def _sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def _marker_digest(components: MarkerComponents) -> str:
    payload = "\0".join(
        (
            components.snp_sequence,
            components.repeat_sequence,
            "" if components.repeat_count_raw is None else f"{components.repeat_count_raw:.12g}",
        )
    )
    return _sequence_digest(payload)


def _panel_digest(loci: list[Locus]) -> str:
    payload = json.dumps(
        [asdict(locus) for locus in sorted(loci, key=lambda item: item.locus_id)],
        sort_keys=True,
        separators=(",", ":"),
    )
    return _sequence_digest(payload)


def _database_stat_digest(database_path: str | Path) -> str:
    path = Path(database_path)
    candidates = (
        sorted(
            item
            for item in path.iterdir()
            if item.is_file() and _is_fasta_path(item)
        )
        if path.is_dir()
        else [path]
    )
    payload = "\n".join(
        f"{item.name}\t{item.stat().st_size}\t{item.stat().st_mtime_ns}"
        for item in candidates
    )
    return _sequence_digest(payload)


def reference_sequence_index_rows(
    references: dict[str, list[tuple[str, str]]],
    loci: list[Locus],
    database_path: str | Path | None = None,
) -> list[dict]:
    """Build compact canonical identity keys for fast reference lookup."""
    locus_by_id = {locus.locus_id: locus for locus in loci}
    panel_sha256 = _panel_digest(loci)
    database_signature = (
        _database_stat_digest(database_path) if database_path is not None else ""
    )
    rows: list[dict] = []
    for locus_id in sorted(references):
        locus = locus_by_id.get(locus_id)
        if locus is None:
            continue
        for reference_id, sequence in references[locus_id]:
            components = decompose_marker_sequence(locus, sequence)
            rows.append(
                {
                    "index_version": "1",
                    "panel_sha256": panel_sha256,
                    "database_signature": database_signature,
                    "locus_id": locus_id,
                    "reference_id": reference_id,
                    "amplicon_sha256": _sequence_digest(components.oriented_sequence),
                    "snp_sha256": _sequence_digest(components.snp_sequence),
                    "marker_sha256": _marker_digest(components),
                    "product_size_bp": len(components.oriented_sequence),
                    "snp_sequence_length": len(components.snp_sequence),
                    "repeat_count_raw": ""
                    if components.repeat_count_raw is None
                    else f"{components.repeat_count_raw:.6f}",
                    "repeat_count": ""
                    if components.repeat_count is None
                    else components.repeat_count,
                }
            )
    return rows


def _read_reference_sequence_index(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(REFERENCE_SEQUENCE_INDEX_FIELDS).issubset(
            reader.fieldnames
        ):
            return []
        return [dict(row) for row in reader]


def _exact_reference_group(
    query_sequences: dict[str, str],
    locus_by_id: dict[str, Locus],
    index_rows: list[dict[str, str]],
) -> tuple[str, list[str], list[str], dict[str, MarkerComponents]]:
    """Return the strongest reference group identical across every callable locus."""
    rows_by_locus: dict[str, list[dict[str, str]]] = {}
    for row in index_rows:
        rows_by_locus.setdefault(str(row["locus_id"]), []).append(row)
    required_loci = sorted(locus_by_id)
    callable_loci = [
        locus_id
        for locus_id in required_loci
        if query_sequences.get(locus_id) and locus_id in rows_by_locus
    ]
    if not required_loci or callable_loci != required_loci:
        return "", [], [], {}
    query_components = {
        locus_id: decompose_marker_sequence(locus_by_id[locus_id], query_sequences[locus_id])
        for locus_id in callable_loci
    }
    for match_type, field, query_key in (
        (
            "EXACT_AMPLICON_MATCH",
            "amplicon_sha256",
            lambda components: _sequence_digest(components.oriented_sequence),
        ),
        ("EXACT_MARKER_MATCH", "marker_sha256", _marker_digest),
    ):
        matching_references: set[str] | None = None
        for locus_id in callable_loci:
            digest = query_key(query_components[locus_id])
            locus_matches = {
                str(row["reference_id"])
                for row in rows_by_locus[locus_id]
                if row.get(field) == digest
            }
            matching_references = (
                locus_matches
                if matching_references is None
                else matching_references & locus_matches
            )
        if matching_references:
            return (
                match_type,
                sorted(matching_references),
                callable_loci,
                query_components,
            )
    return "", [], callable_loci, query_components


@dataclass
class _Node:
    name: str | None
    children: list[tuple["_Node", float]]
    edge_num: int | None = None


def canonical_assembly_digest(path: str | Path) -> str:
    """Hash assembly sequence content independent of headers/order/orientation."""
    canonical_contigs = sorted(
        min(sequence, _reverse_complement(sequence))
        for _name, sequence in _read_fasta(path)
    )
    payload = "".join(f"{len(sequence)}:{sequence};" for sequence in canonical_contigs)
    return _sequence_digest(payload)


def check_dnadiff(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"dnadiff executable {executable!r} was not found. Install mummer4 "
            "from Bioconda or pass --dnadiff-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run dnadiff at {path}: {detail}")
    return path


def _parse_dnadiff_report(path: Path) -> dict[str, float | int]:
    section = ""
    values: dict[str, float | int] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue
        parts = line.split()
        if not parts:
            continue
        if section == "Bases" and parts[0] == "AlignedBases" and len(parts) >= 3:
            ref_match = re.search(r"\(([0-9.]+)%\)", parts[1])
            query_match = re.search(r"\(([0-9.]+)%\)", parts[2])
            if ref_match and query_match:
                values["align_fraction_ref"] = float(ref_match.group(1))
                values["align_fraction_query"] = float(query_match.group(1))
        elif section == "SNPs" and parts[0] == "TotalSNPs" and len(parts) >= 2:
            values["snps"] = int(parts[1])
        elif section == "SNPs" and parts[0] == "TotalIndels" and len(parts) >= 2:
            values["indel_bases"] = int(parts[1])
    required = {
        "snps",
        "indel_bases",
        "align_fraction_ref",
        "align_fraction_query",
    }
    missing = sorted(required - set(values))
    if missing:
        raise RuntimeError(f"Could not parse {', '.join(missing)} from {path}")
    return values


def _run_dnadiff_comparison(
    dnadiff: str,
    reference_id: str,
    reference_path: Path,
    query_path: Path,
    output_dir: Path,
) -> tuple[str, dict[str, float | int | bool], Path]:
    safe_reference = _SAFE_FILE.sub("_", reference_id).strip("_") or "reference"
    comparison_dir = output_dir / f"{safe_reference}__{_sequence_digest(reference_id)[:12]}"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    prefix = comparison_dir / "dnadiff"
    result = subprocess.run(
        [dnadiff, "-p", str(prefix), str(reference_path), str(query_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(
            f"dnadiff comparison for {reference_id!r} failed "
            f"(exit {result.returncode}): {detail}"
        )
    report_path = Path(f"{prefix}.report")
    if not report_path.exists():
        raise RuntimeError(
            f"dnadiff comparison for {reference_id!r} completed without {report_path}"
        )
    metrics: dict[str, float | int | bool] = {
        **_parse_dnadiff_report(report_path),
        "exact_genome_match": False,
    }
    return reference_id, metrics, report_path


def _dnadiff_tie_break(
    query_assembly: str | Path,
    reference_ids: list[str],
    sequence_database_path: Path,
    output_dir: Path,
    threads: int,
    executable: str,
) -> tuple[dict[str, dict[str, float | int | bool]], Path]:
    """Compare tied references by exact whole-genome 1-to-1 MUMmer alignments."""
    mapping_path = sequence_database_path / "reference_assemblies.tsv"
    interpreted_path = output_dir / "whole_genome_dnadiff.tsv"
    if not mapping_path.exists():
        return {}, interpreted_path

    with mapping_path.open(newline="") as handle:
        mapping_rows = list(csv.DictReader(handle, delimiter="\t"))
    wanted = set(reference_ids)
    references: dict[str, tuple[Path, str]] = {}
    for row in mapping_rows:
        reference_id = str(row.get("reference_id", ""))
        assembly_file = str(row.get("assembly_file", ""))
        if reference_id not in wanted or not assembly_file:
            continue
        reference_path = Path(assembly_file)
        references[reference_id] = (
            reference_path,
            str(row.get("assembly_sha256", "")),
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    query_path = Path(query_assembly).resolve()
    query_digest = canonical_assembly_digest(query_path)
    metrics: dict[str, dict[str, float | int | bool]] = {}
    report_paths: dict[str, Path] = {}
    comparisons: list[tuple[str, Path]] = []
    for reference_id, (reference_path, stored_digest) in references.items():
        reference_digest = stored_digest or (
            canonical_assembly_digest(reference_path) if reference_path.exists() else ""
        )
        if reference_digest == query_digest:
            metrics[reference_id] = {
                "exact_genome_match": True,
                "snps": 0,
                "indel_bases": 0,
                "align_fraction_ref": 100.0,
                "align_fraction_query": 100.0,
            }
        elif reference_path.exists():
            comparisons.append((reference_id, reference_path))

    if comparisons:
        dnadiff = check_dnadiff(executable)
        worker_count = min(max(1, threads), len(comparisons))
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(
                    _run_dnadiff_comparison,
                    dnadiff,
                    reference_id,
                    reference_path,
                    query_path,
                    output_dir / "dnadiff",
                ): reference_id
                for reference_id, reference_path in comparisons
            }
            for future in as_completed(futures):
                reference_id, values, report_path = future.result()
                metrics[reference_id] = values
                report_paths[reference_id] = report_path

    result_rows: list[dict] = []
    for reference_id in sorted(metrics):
        values = metrics[reference_id]
        reference_path = references[reference_id][0]
        result_rows.append(
            {
                "reference_id": reference_id,
                "reference_file": str(reference_path),
                "query_file": str(query_path),
                "exact_genome_match": "yes"
                if values["exact_genome_match"]
                else "no",
                "snps": values["snps"],
                "indel_bases": values["indel_bases"],
                "align_fraction_ref": f"{float(values['align_fraction_ref']):.8f}",
                "align_fraction_query": f"{float(values['align_fraction_query']):.8f}",
                "report_file": str(report_paths.get(reference_id, "")),
            }
        )
    _write_tsv(result_rows, interpreted_path, DNADIFF_RESULT_FIELDS)
    return metrics, interpreted_path


def check_mafft(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"MAFFT executable {executable!r} was not found. Install mafft from "
            "Bioconda or pass --mafft-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run MAFFT at {path}: {detail}")
    return path


def build_mafft_reference_command(
    reference_fasta: str | Path, threads: int, executable: str = "mafft"
) -> list[str]:
    return [executable, "--auto", "--thread", str(threads), str(reference_fasta)]


def build_mafft_add_command(
    query_fasta: str | Path,
    reference_alignment: str | Path,
    threads: int,
    executable: str = "mafft",
) -> list[str]:
    return [
        executable,
        "--add",
        str(query_fasta),
        "--keeplength",
        "--thread",
        str(threads),
        str(reference_alignment),
    ]


def _run_mafft(command: list[str], output_path: Path, stage: str) -> None:
    with output_path.open("w") as output:
        result = subprocess.run(
            command,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if result.returncode:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"MAFFT {stage} failed (exit {result.returncode}): "
            f"{(result.stderr or '').strip()}"
        )


def check_raxml_ng(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"RAxML-NG executable {executable!r} was not found. Install raxml-ng "
            "from Bioconda or pass --raxml-ng-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run RAxML-NG at {path}: {detail}")
    return path


def build_raxml_ng_command(
    alignment_path: str | Path,
    prefix: str | Path,
    threads: int,
    executable: str = "raxml-ng",
    model: str = "DNA",
) -> list[str]:
    return [
        executable,
        "--search",
        "--msa",
        str(alignment_path),
        "--model",
        model,
        "--prefix",
        str(prefix),
        "--seed",
        "12345",
        "--threads",
        str(threads),
        "--redo",
    ]


def _run_raxml_ng(
    command: list[str],
    prefix: Path,
    output_tree: Path,
    stage: str,
    progress: ProgressReporter | None = None,
) -> None:
    """Run RAxML-NG, reducing only its threads for low-pattern alignments."""
    attempted_command = list(command)
    thread_index = attempted_command.index("--threads") + 1
    raxml_threads = int(attempted_command[thread_index])
    log_path = Path(f"{prefix}.mlvamaps.raxml.log")
    attempts: list[str] = []

    while True:
        attempted_command[thread_index] = str(raxml_threads)
        result = subprocess.run(
            attempted_command, capture_output=True, text=True, check=False
        )
        detail = "\n".join(
            part.strip() for part in (result.stdout, result.stderr) if part.strip()
        )
        attempts.append(
            f"$ {' '.join(attempted_command)}\n{detail}".rstrip()
        )
        log_path.write_text("\n\n".join(attempts) + "\n")
        if not result.returncode:
            break
        if "Too few patterns per thread" not in detail or raxml_threads == 1:
            detail_tail = "\n".join(detail.splitlines()[-40:])
            raise RuntimeError(
                f"RAxML-NG {stage} failed (exit {result.returncode}). Full output: "
                f"{log_path}\n{detail_tail}"
            )
        next_threads = max(1, raxml_threads // 2)
        if progress is not None:
            progress.step(
                f"RAxML-NG found too few alignment patterns for {raxml_threads} "
                f"threads; retrying {stage} with {next_threads}"
            )
        raxml_threads = next_threads

    best_tree = Path(f"{prefix}.raxml.bestTree")
    if not best_tree.exists():
        raise RuntimeError(
            f"RAxML-NG {stage} completed without producing {best_tree}"
        )
    shutil.copyfile(best_tree, output_tree)


def check_epa_ng(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"EPA-ng executable {executable!r} was not found. Install epa-ng "
            "from Bioconda or pass --epa-ng-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run EPA-ng at {path}: {detail}")
    return path


def build_epa_ng_command(
    reference_alignment: str | Path,
    reference_tree: str | Path,
    query_alignment: str | Path,
    model_path: str | Path,
    outdir: str | Path,
    threads: int,
    executable: str = "epa-ng",
) -> list[str]:
    return [
        executable,
        "--ref-msa",
        str(reference_alignment),
        "--tree",
        str(reference_tree),
        "--query",
        str(query_alignment),
        "--model",
        str(model_path),
        "--outdir",
        str(outdir),
        "--threads",
        str(threads),
    ]


def _run_epa_ng(command: list[str], outdir: Path, stage: str) -> Path:
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"EPA-ng {stage} failed (exit {result.returncode}): {detail}")
    jplace = outdir / "epa_result.jplace"
    if not jplace.exists():
        raise RuntimeError(f"EPA-ng {stage} completed without producing {jplace}")
    return jplace


@dataclass(frozen=True)
class _ReferenceTreeJob:
    locus_id: str
    reference_fasta: Path
    reference_alignment: Path
    reference_prefix: Path
    reference_tree: Path
    reference_model: Path


def _run_reference_tree_job(
    job: _ReferenceTreeJob,
    mafft: str,
    raxml_ng: str,
    raxml_model: str,
    native_threads: int,
) -> None:
    """Build one locus tree; callers parallelize independent loci."""
    _run_mafft(
        build_mafft_reference_command(
            job.reference_fasta, native_threads, mafft
        ),
        job.reference_alignment,
        f"reference alignment for {job.locus_id}",
    )
    _run_raxml_ng(
        build_raxml_ng_command(
            job.reference_alignment,
            job.reference_prefix,
            native_threads,
            raxml_ng,
            raxml_model,
        ),
        job.reference_prefix,
        job.reference_tree,
        f"reference tree search for {job.locus_id}",
    )
    if not job.reference_model.exists():
        raise RuntimeError(
            "RAxML-NG reference search did not produce model file "
            f"{job.reference_model}"
        )


@dataclass(frozen=True)
class _PlacementJob:
    locus_id: str
    query_name: str
    query_fasta: Path
    reference_alignment: Path
    reference_tree: Path
    reference_model: Path
    placed_alignment: Path
    query_alignment: Path
    epa_outdir: Path


def _run_placement_job(
    job: _PlacementJob,
    mafft: str,
    epa_ng: str,
    native_threads: int = 1,
) -> Path:
    """Add and place one query; jobs are parallelized across independent loci."""
    _run_mafft(
        build_mafft_add_command(
            job.query_fasta, job.reference_alignment, native_threads, mafft
        ),
        job.placed_alignment,
        f"--add placement for {job.locus_id}",
    )
    aligned_query_records = [
        (name, sequence)
        for name, sequence in _read_fasta(job.placed_alignment)
        if name == job.query_name
    ]
    if len(aligned_query_records) != 1:
        raise RuntimeError(
            f"MAFFT placement alignment for {job.locus_id} did not contain exactly "
            f"one {job.query_name!r}"
        )
    _write_fasta(aligned_query_records, job.query_alignment)
    return _run_epa_ng(
        build_epa_ng_command(
            job.reference_alignment,
            job.reference_tree,
            job.query_alignment,
            job.reference_model,
            job.epa_outdir,
            native_threads,
            epa_ng,
        ),
        job.epa_outdir,
        f"placement for {job.locus_id}",
    )


def _read_database_fasta(path: Path, locus_ids: set[str]) -> dict[str, list[tuple[str, str]]]:
    records = _read_fasta(path)
    if path.parent != path and _fasta_stem(path) in locus_ids:
        return {_fasta_stem(path): records}
    by_locus: dict[str, list[tuple[str, str]]] = {}
    for header, sequence in records:
        parts = header.split("|")
        matching = [part for part in parts if part in locus_ids]
        if len(matching) != 1:
            raise ValueError(
                f"Could not identify one panel locus in FASTA header {header!r}. "
                "Use LOCUS.fasta files or headers such as reference_id|LOCUS."
            )
        locus_id = matching[0]
        reference_id = next((part for part in parts if part != locus_id), "")
        if not reference_id:
            raise ValueError(f"No reference id was present in FASTA header {header!r}")
        by_locus.setdefault(locus_id, []).append((reference_id, sequence))
    return by_locus


def read_sequence_database(
    database_path: str | Path, locus_ids: set[str]
) -> dict[str, list[tuple[str, str]]]:
    """Read per-locus references from a directory, FASTA, or long-form TSV."""
    path = Path(database_path)
    if not path.exists():
        raise ValueError(f"Sequence database path does not exist: {path}")
    by_locus: dict[str, list[tuple[str, str]]] = {}
    if path.is_dir():
        fasta_paths = sorted(
            item for item in path.iterdir() if item.is_file() and _is_fasta_path(item)
        )
        if not fasta_paths:
            raise ValueError(f"Sequence database directory contains no FASTA files: {path}")
        for fasta_path in fasta_paths:
            locus_name = _fasta_stem(fasta_path)
            if locus_name not in locus_ids:
                continue
            by_locus[locus_name] = _read_fasta(fasta_path)
    elif _is_fasta_path(path):
        by_locus = _read_database_fasta(path, locus_ids)
    else:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"reference_id", "locus_id", "sequence"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(
                    "Sequence database TSV requires reference_id, locus_id, and sequence columns"
                )
            for row in reader:
                locus_id = row["locus_id"]
                if locus_id in locus_ids:
                    by_locus.setdefault(locus_id, []).append(
                        (row["reference_id"], row["sequence"].upper())
                    )
    for locus_id, records in by_locus.items():
        names = [name for name, _sequence in records]
        if not records:
            raise ValueError(f"No reference sequences found for locus {locus_id!r}")
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate reference ids found for locus {locus_id!r}")
        if any(not name or not sequence for name, sequence in records):
            raise ValueError(f"Blank reference id or sequence found for locus {locus_id!r}")
    if not by_locus:
        raise ValueError("Sequence database contains no loci matching the supplied panel")
    return by_locus


def _resolve_reference_database_layout(
    database_path: str | Path,
) -> tuple[Path, Path | None]:
    """Resolve sequence data and optional reusable trees from a reference build."""
    path = Path(database_path)
    bundled_database = path / "database"
    bundled_phylogeny = path / "phylogeny"
    if path.is_dir() and bundled_database.is_dir():
        return bundled_database, bundled_phylogeny if bundled_phylogeny.is_dir() else None
    sibling_phylogeny = path.parent / "phylogeny"
    if path.is_dir() and path.name == "database" and sibling_phylogeny.is_dir():
        return path, sibling_phylogeny
    return path, None


def _copy_reference_artifact(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != destination.resolve():
        opener = gzip.open if source.suffix.lower() == ".gz" else open
        with opener(source, "rb") as input_handle, destination.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle)


def read_reference_metadata(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    metadata_path = Path(path)
    if not metadata_path.exists():
        raise ValueError(f"Reference metadata path does not exist: {metadata_path}")
    with metadata_path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        delimiter = "\t" if "\t" in first_line else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or "reference_id" not in reader.fieldnames:
            raise ValueError("Reference metadata requires a reference_id column")
        rows = list(reader)
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        reference_id = str(row.get("reference_id", "")).strip()
        if not reference_id:
            continue
        if reference_id in result:
            raise ValueError(f"Duplicate reference metadata for {reference_id!r}")
        result[reference_id] = {str(key): str(value or "") for key, value in row.items()}
    return result


def _parse_newick(text: str) -> _Node:
    """Parse the branch lengths and leaf names needed from a RAxML Newick tree."""
    text = text.strip()
    position = 0

    def skip_space() -> None:
        nonlocal position
        while position < len(text) and text[position].isspace():
            position += 1

    def label() -> str:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] in {"'", '"'}:
            quote = text[position]
            position += 1
            parts = []
            while position < len(text):
                if text[position] == quote:
                    if position + 1 < len(text) and text[position + 1] == quote:
                        parts.append(quote)
                        position += 2
                        continue
                    position += 1
                    return "".join(parts)
                parts.append(text[position])
                position += 1
            raise ValueError("Unterminated quoted name in RAxML-NG Newick tree")
        start = position
        while position < len(text) and text[position] not in ",():;":
            position += 1
        return text[start:position].strip()

    def branch_length() -> tuple[float, int | None]:
        nonlocal position
        skip_space()
        if position >= len(text) or text[position] != ":":
            return 0.0, None
        position += 1
        start = position
        while position < len(text) and text[position] not in ",);{":
            position += 1
        value = text[start:position].strip()
        try:
            length = float(value)
        except ValueError as exc:
            raise ValueError(f"Invalid branch length {value!r} in RAxML-NG tree") from exc
        edge_num = None
        if position < len(text) and text[position] == "{":
            position += 1
            edge_start = position
            while position < len(text) and text[position] != "}":
                position += 1
            if position >= len(text):
                raise ValueError("Unterminated edge number in jplace tree")
            try:
                edge_num = int(text[edge_start:position])
            except ValueError as exc:
                raise ValueError("Invalid edge number in jplace tree") from exc
            position += 1
        return length, edge_num

    def subtree() -> tuple[_Node, float]:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] == "(":
            position += 1
            children: list[tuple[_Node, float]] = []
            while True:
                children.append(subtree())
                skip_space()
                if position < len(text) and text[position] == ",":
                    position += 1
                    continue
                if position >= len(text) or text[position] != ")":
                    raise ValueError("Malformed RAxML-NG Newick tree")
                position += 1
                break
            # RAxML may put a support value or internal-node label here. It is
            # irrelevant for patristic distance but must still be consumed.
            internal_name = label()
            node = _Node(internal_name or None, children)
            length, node.edge_num = branch_length()
            return node, length
        name = label()
        if not name:
            raise ValueError("Blank leaf name in RAxML-NG Newick tree")
        node = _Node(name, [])
        length, node.edge_num = branch_length()
        return node, length

    root, _root_length = subtree()
    skip_space()
    if position < len(text) and text[position] == ";":
        position += 1
    skip_space()
    if position != len(text):
        raise ValueError("Unexpected trailing content in RAxML-NG Newick tree")
    return root


def _placement_patristic_distances(
    root: _Node,
    placement_edge: int,
    distal_length: float,
    pendant_length: float,
) -> dict[str, float]:
    adjacency: dict[int, list[tuple[_Node, float]]] = {}
    target: tuple[_Node, _Node, float] | None = None

    def connect(node: _Node) -> None:
        nonlocal target
        adjacency.setdefault(id(node), [])
        for child, length in node.children:
            adjacency[id(node)].append((child, length))
            adjacency.setdefault(id(child), []).append((node, length))
            if child.edge_num == placement_edge:
                target = (node, child, length)
            connect(child)

    connect(root)
    if target is None:
        raise ValueError(f"Placement edge {placement_edge} is absent from the jplace tree")
    parent, child, edge_length = target
    distal = min(max(0.0, distal_length), edge_length)
    distances: dict[str, float] = {}

    def walk(node: _Node, blocked: _Node, total: float) -> None:
        if not node.children and node.name is not None:
            distances[node.name] = total
        for neighbor, length in adjacency[id(node)]:
            if neighbor is not blocked:
                walk(neighbor, node, total + length)

    walk(child, parent, pendant_length + distal)
    walk(parent, child, pendant_length + edge_length - distal)
    return distances


def _reference_tree_scale(root: _Node) -> float:
    """Median non-zero patristic distance among reference tips."""
    _names, matrix = _tip_patristic_distance_matrix(root)
    positive = matrix[np.triu_indices(len(matrix), k=1)]
    positive = positive[positive > 0]
    return float(np.median(positive)) if positive.size else 1.0


def _tip_patristic_distance_matrix(root: _Node) -> tuple[list[str], np.ndarray]:
    """Build a dense tip-distance matrix using NumPy's compiled matrix kernels."""
    edge_lengths: list[float] = []
    names: list[str] = []
    tip_paths: list[list[int]] = []

    def collect(node: _Node, path: list[int]) -> None:
        if not node.children and node.name is not None:
            names.append(str(node.name))
            tip_paths.append(path)
        for child, length in node.children:
            edge_index = len(edge_lengths)
            edge_lengths.append(float(length))
            collect(child, [*path, edge_index])

    collect(root, [])
    if len(names) != len(set(names)):
        raise ValueError("Reference tree contains duplicate tip names")
    incidence = np.zeros((len(names), len(edge_lengths)), dtype=np.float64)
    for tip_index, path in enumerate(tip_paths):
        incidence[tip_index, path] = 1.0
    weighted_paths = incidence * np.asarray(edge_lengths, dtype=np.float64)
    root_distances = weighted_paths.sum(axis=1)
    shared_distances = weighted_paths @ incidence.T
    matrix = (
        root_distances[:, None]
        + root_distances[None, :]
        - 2.0 * shared_distances
    )
    np.maximum(matrix, 0.0, out=matrix)
    return names, matrix


def _tip_patristic_distances(root: _Node) -> dict[tuple[str, str], float]:
    """Return the public pair-keyed view of the compiled dense calculation."""
    names, matrix = _tip_patristic_distance_matrix(root)
    return {
        tuple(sorted((left, right))): float(matrix[left_index, right_index])
        for left_index, left in enumerate(names)
        for right_index, right in enumerate(names[left_index + 1 :], left_index + 1)
    }


def _newick_label(label: str) -> str:
    # Newick escapes a quote inside a quoted label by doubling it. Preserve
    # sample IDs exactly instead of silently rewriting meaningful characters.
    return "'" + str(label).replace("'", "''") + "'"


def neighbor_joining_tree(
    labels: list[str], distances: dict[tuple[str, str], float]
) -> str:
    """Build a deterministic Newick neighbor-joining tree from a distance matrix."""
    matrix = np.zeros((len(labels), len(labels)), dtype=np.float64)

    def key(left: str, right: str) -> tuple[str, str]:
        return tuple(sorted((left, right)))

    for left_index, left in enumerate(labels):
        for right_index, right in enumerate(
            labels[left_index + 1 :], start=left_index + 1
        ):
            pair = key(left, right)
            if pair not in distances:
                raise ValueError(f"Missing distance between {left!r} and {right!r}")
            value = max(0.0, float(distances[pair]))
            matrix[left_index, right_index] = value
            matrix[right_index, left_index] = value
    return _neighbor_joining_tree_from_matrix(labels, matrix)


def neighbor_joining_tree_from_matrix(
    labels: list[str], distances: np.ndarray
) -> str:
    """Build deterministic Newick directly from a finite dense distance matrix."""
    return _neighbor_joining_tree_from_matrix(labels, distances)


def _neighbor_joining_tree_from_matrix(
    labels: list[str], distances: np.ndarray
) -> str:
    """Run neighbor joining on a dense matrix through compiled NumPy kernels."""
    if not labels:
        return ";\n"
    if len(labels) == 1:
        return f"({_newick_label(labels[0])}:0.00000000);\n"
    active = list(labels)
    clusters = {label: _newick_label(label) for label in labels}
    matrix = np.asarray(distances, dtype=np.float64).copy()
    if matrix.shape != (len(labels), len(labels)):
        raise ValueError("Neighbor-joining distance matrix shape does not match labels")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Neighbor-joining distances must be finite")
    np.maximum(matrix, 0.0, out=matrix)
    np.fill_diagonal(matrix, 0.0)

    node_number = 0
    while len(active) > 2:
        size = len(active)
        row_sums = matrix.sum(axis=1)
        q_matrix = (size - 2) * matrix - row_sums[:, None] - row_sums[None, :]
        q_matrix[np.tril_indices(size)] = np.inf
        left_index, right_index = np.unravel_index(
            int(np.argmin(q_matrix)), q_matrix.shape
        )
        left = active[left_index]
        right = active[right_index]
        pair_distance = float(matrix[left_index, right_index])
        left_length = 0.5 * pair_distance + (
            float(row_sums[left_index]) - float(row_sums[right_index])
        ) / (2 * (size - 2))
        right_length = pair_distance - left_length
        left_length = max(0.0, left_length)
        right_length = max(0.0, right_length)
        node_number += 1
        joined = f"__NJ_{node_number}"
        clusters[joined] = (
            f"({clusters[left]}:{left_length:.8f},"
            f"{clusters[right]}:{right_length:.8f})"
        )
        retained = [
            index for index in range(size) if index not in {left_index, right_index}
        ]
        joined_distances = np.maximum(
            0.0,
            (
                matrix[left_index, retained]
                + matrix[right_index, retained]
                - pair_distance
            )
            / 2,
        )
        next_matrix = np.zeros((size - 1, size - 1), dtype=float)
        next_matrix[:-1, :-1] = matrix[np.ix_(retained, retained)]
        next_matrix[-1, :-1] = joined_distances
        next_matrix[:-1, -1] = joined_distances
        matrix = next_matrix
        active = [active[index] for index in retained]
        active.append(joined)

    left, right = active
    final_distance = float(matrix[0, 1]) / 2
    return (
        f"({clusters[left]}:{max(0.0, final_distance):.8f},"
        f"{clusters[right]}:{max(0.0, final_distance):.8f});\n"
    )


def read_epa_ng_placement(jplace_path: str | Path) -> tuple[dict[str, float], dict[str, float | int]]:
    best_distances, _expected_distances, metadata = read_epa_ng_placement_statistics(
        jplace_path
    )
    return best_distances, metadata


def read_epa_ng_placement_statistics(
    jplace_path: str | Path,
) -> tuple[dict[str, float], dict[str, float], dict[str, float | int]]:
    """Return best-placement and LWR-weighted expected tip distances."""
    with Path(jplace_path).open() as handle:
        document = json.load(handle)
    fields = list(document.get("fields", []))
    placements = document.get("placements", [])
    if not placements:
        raise ValueError(f"EPA-ng produced no placements in {jplace_path}")
    required = {"edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"}
    if not required.issubset(fields):
        raise ValueError(f"EPA-ng jplace is missing required fields: {sorted(required - set(fields))}")
    indexes = {field: fields.index(field) for field in required}
    candidates = placements[0].get("p", [])
    if not candidates:
        raise ValueError(f"EPA-ng produced an empty placement record in {jplace_path}")
    raw_weights = [max(0.0, float(row[indexes["like_weight_ratio"]])) for row in candidates]
    weight_total = sum(raw_weights)
    if weight_total <= 0:
        weights = [1.0 / len(candidates)] * len(candidates)
    else:
        weights = [weight / weight_total for weight in raw_weights]
    best_index = max(range(len(candidates)), key=lambda index: weights[index])
    best = candidates[best_index]
    edge_num = int(best[indexes["edge_num"]])
    distal_length = float(best[indexes["distal_length"]])
    pendant_length = float(best[indexes["pendant_length"]])
    root = _parse_newick(str(document["tree"]))
    best_distances = _placement_patristic_distances(
        root, edge_num, distal_length, pendant_length
    )
    expected_distances = {reference_id: 0.0 for reference_id in best_distances}
    for row, weight in zip(candidates, weights):
        candidate_distances = _placement_patristic_distances(
            root,
            int(row[indexes["edge_num"]]),
            float(row[indexes["distal_length"]]),
            float(row[indexes["pendant_length"]]),
        )
        if set(candidate_distances) != set(expected_distances):
            raise ValueError("EPA-ng placement candidates refer to inconsistent trees")
        for reference_id, distance in candidate_distances.items():
            expected_distances[reference_id] += weight * distance
    entropy = -sum(weight * math.log(weight) for weight in weights if weight > 0)
    metadata: dict[str, float | int] = {
        "placement_edge": edge_num,
        "like_weight_ratio": float(best[indexes["like_weight_ratio"]]),
        "pendant_length": pendant_length,
        "distal_length": distal_length,
        "placement_count": len(candidates),
        "placement_entropy": entropy,
        "reference_tree_scale": _reference_tree_scale(root),
    }
    return best_distances, expected_distances, metadata


def _write_exact_match_outputs(
    outdir: str | Path,
    sample_id: str,
    match_type: str,
    reference_ids: list[str],
    locus_ids: list[str],
    locus_by_id: dict[str, Locus],
    query_components: dict[str, MarkerComponents],
    index_rows: list[dict[str, str]],
    reference_metadata: dict[str, dict[str, str]],
    snp_weight: float,
    repeat_weight: float,
    progress: ProgressReporter | None,
    whole_genome_metrics: dict[str, dict[str, float | int | bool]] | None = None,
    dnadiff_result_path: Path | None = None,
    dnadiff_available: bool = False,
    dnadiff_applicable: bool = False,
) -> dict[str, Path]:
    """Write normal placement outputs for a deterministic indexed identity hit."""
    whole_genome_metrics = whole_genome_metrics or {}
    if len(reference_ids) > 1 and whole_genome_metrics:
        def tie_key(reference_id: str) -> tuple[bool, float, float, float, float]:
            values = whole_genome_metrics.get(reference_id)
            if values is None:
                return (True, float("inf"), float("inf"), float("inf"), float("inf"))
            ref_af = float(values["align_fraction_ref"])
            query_af = float(values["align_fraction_query"])
            return (
                not bool(values["exact_genome_match"]),
                float(values["snps"]),
                float(values["indel_bases"]),
                -min(ref_af, query_af),
                -max(ref_af, query_af),
            )

        reference_ids = sorted(reference_ids, key=lambda item: (tie_key(item), item))
        rank_by_reference: dict[str, int] = {}
        previous_key: tuple[bool, float, float, float, float] | None = None
        previous_rank = 0
        for position, reference_id in enumerate(reference_ids, start=1):
            current_key = tie_key(reference_id)
            if current_key != previous_key:
                previous_rank = position
                previous_key = current_key
            rank_by_reference[reference_id] = previous_rank
    else:
        rank_by_reference = {reference_id: 1 for reference_id in reference_ids}

    if len(reference_ids) <= 1:
        tie_break_status = "NOT_NEEDED"
        ranking_warning = ""
    elif not dnadiff_applicable:
        tie_break_status = "NOT_APPLICABLE"
        ranking_warning = ""
    elif len(whole_genome_metrics) == len(reference_ids):
        tie_break_status = "APPLIED"
        ranking_warning = ""
    elif dnadiff_available:
        tie_break_status = "INCOMPLETE"
        ranking_warning = "DNADIFF_TIE_BREAK_INCOMPLETE"
    else:
        tie_break_status = "UNAVAILABLE"
        ranking_warning = "DNADIFF_TIE_BREAK_UNAVAILABLE"
    output = Path(outdir) / "phylogeny"
    output.mkdir(parents=True, exist_ok=True)
    row_lookup = {
        (str(row["locus_id"]), str(row["reference_id"])): row
        for row in index_rows
    }
    query_name = f"QUERY__{_SAFE_FILE.sub('_', sample_id).strip('_') or 'sample'}"
    detail_rows: list[dict] = []
    marker_rows: list[dict] = []
    locus_marker_rows: list[dict] = []
    status_rows: list[dict] = []
    for locus_id in locus_ids:
        components = query_components[locus_id]
        marker_rows.append(
            _marker_component_row(
                sample_id,
                locus_by_id[locus_id],
                "query",
                query_name,
                components,
            )
        )
        status_rows.append(
            {
                "locus_id": locus_id,
                "reference_sequences": len(reference_ids),
                "query_sequence": "yes",
                "status": match_type,
            }
        )
        for reference_id in reference_ids:
            index_row = row_lookup[(locus_id, reference_id)]
            marker_rows.append(
                {
                    "sample_id": sample_id,
                    "locus_id": locus_id,
                    "record_type": "reference",
                    "reference_id": reference_id,
                    "repeat_count_raw": index_row.get("repeat_count_raw", ""),
                    "repeat_count": index_row.get("repeat_count", ""),
                    "repeat_sequence": "",
                    "repeat_haplotype": "",
                    "repeat_region_start": "",
                    "repeat_region_end": "",
                    "snp_sequence_length": index_row.get("snp_sequence_length", ""),
                    "masking_method": "indexed_exact_match",
                }
            )
            detail_rows.append(
                {
                    "sample_id": sample_id,
                    "locus_id": locus_id,
                    "reference_id": reference_id,
                    "phylogenetic_distance": "0.00000000",
                    "likelihood_weighted_phylogenetic_distance": "0.00000000",
                    "placement_edge": "",
                    "like_weight_ratio": "1.00000000",
                    "pendant_length": "0.00000000",
                    "distal_length": "0.00000000",
                    "placement_count": 0,
                    "placement_entropy": "0.00000000",
                    "reference_tree_scale": "1.00000000",
                }
            )
            reference_count = index_row.get("repeat_count_raw", "")
            locus_marker_rows.append(
                {
                    "sample_id": sample_id,
                    "locus_id": locus_id,
                    "reference_id": reference_id,
                    "query_repeat_count": ""
                    if components.repeat_count_raw is None
                    else f"{components.repeat_count_raw:.6f}",
                    "reference_repeat_count": reference_count,
                    "repeat_count_delta": "0.000000",
                    "normalized_repeat_distance": "0.00000000",
                    "likelihood_weighted_snp_distance": "0.00000000",
                    "reference_tree_scale": "1.00000000",
                    "placement_normalized_snp_distance": "0.00000000",
                    "direct_snp_distance": "0.00000000",
                    "direct_snp_scale": "1.00000000",
                    "normalized_direct_snp_distance": "0.00000000",
                    "exact_snp_match": "yes",
                    "normalized_snp_distance": "0.00000000",
                    "placement_entropy": "0.00000000",
                }
            )

    detail_path = output / "locus_phylogenetic_distances.tsv"
    marker_components_path = output / "marker_components.tsv"
    locus_marker_path = output / "locus_marker_distances.tsv"
    status_path = output / "locus_status.tsv"
    _write_tsv(detail_rows, detail_path, PLACEMENT_FIELDS)
    _write_tsv(marker_rows, marker_components_path, MARKER_COMPONENT_FIELDS)
    _write_tsv(locus_marker_rows, locus_marker_path, LOCUS_MARKER_DISTANCE_FIELDS)
    _write_tsv(status_rows, status_path, STATUS_FIELDS)

    summary_rows = [
        {
            "sample_id": sample_id,
            "reference_id": reference_id,
            "total_phylogenetic_distance": "0.00000000",
            "total_likelihood_weighted_distance": "0.00000000",
            "compared_loci": len(locus_ids),
            "mean_phylogenetic_distance": "0.00000000",
            "mean_likelihood_weighted_distance": "0.00000000",
            "distance_gap_to_next": "0.00000000"
            if index + 1 < len(reference_ids)
            else "",
            "relative_distance_gap_to_next": "0.00000000"
            if index + 1 < len(reference_ids)
            else "",
            "rank": rank_by_reference[reference_id],
        }
        for index, reference_id in enumerate(reference_ids)
    ]
    summary_path = output / "phylogenetic_matches.tsv"
    _write_tsv(summary_rows, summary_path, SUMMARY_FIELDS)

    combined_rows = []
    for index, reference_id in enumerate(reference_ids):
        metadata = reference_metadata.get(reference_id, {})
        combined_rows.append(
            {
                "sample_id": sample_id,
                "reference_id": reference_id,
                "total_likelihood_weighted_snp_distance": "0.00000000",
                "total_placement_normalized_snp_distance": "0.00000000",
                "total_normalized_direct_snp_distance": "0.00000000",
                "total_normalized_snp_distance": "0.00000000",
                "total_repeat_count_distance": "0.00000000",
                "total_normalized_repeat_distance": "0.00000000",
                "snp_weight": f"{snp_weight:.6f}",
                "repeat_weight": f"{repeat_weight:.6f}",
                "combined_marker_distance": "0.00000000",
                "compared_loci": len(locus_ids),
                "repeat_compared_loci": sum(
                    query_components[locus_id].repeat_count_raw is not None
                    for locus_id in locus_ids
                ),
                "exact_snp_loci": len(locus_ids),
                "exact_marker_loci": len(locus_ids),
                "match_status": match_type,
                "ranking_warning": ranking_warning,
                "whole_genome_exact_match": ""
                if reference_id not in whole_genome_metrics
                else (
                    "yes"
                    if whole_genome_metrics[reference_id]["exact_genome_match"]
                    else "no"
                ),
                "whole_genome_snps": ""
                if reference_id not in whole_genome_metrics
                else whole_genome_metrics[reference_id]["snps"],
                "whole_genome_indel_bases": ""
                if reference_id not in whole_genome_metrics
                else whole_genome_metrics[reference_id]["indel_bases"],
                "whole_genome_align_fraction_ref": ""
                if reference_id not in whole_genome_metrics
                else f"{float(whole_genome_metrics[reference_id]['align_fraction_ref']):.8f}",
                "whole_genome_align_fraction_query": ""
                if reference_id not in whole_genome_metrics
                else f"{float(whole_genome_metrics[reference_id]['align_fraction_query']):.8f}",
                "tie_break_method": "canonical_genome_hash_then_dnadiff_snps_indels_alignment"
                if len(reference_ids) > 1 and dnadiff_applicable
                else "",
                "tie_break_status": tie_break_status,
                "distance_gap_to_next": "0.00000000"
                if index + 1 < len(reference_ids)
                else "",
                "relative_distance_gap_to_next": "0.00000000"
                if index + 1 < len(reference_ids)
                else "",
                "collection_date": metadata.get("collection_date", ""),
                "latitude": metadata.get("latitude", ""),
                "longitude": metadata.get("longitude", ""),
                "location": metadata.get("location", ""),
                "source": metadata.get("source", ""),
                "rank": rank_by_reference[reference_id],
            }
        )
    combined_path = output / "combined_marker_matches.tsv"
    _write_tsv(combined_rows, combined_path, COMBINED_MARKER_FIELDS)

    closest_reference_bands_path = output / "closest_reference_bands.tsv"
    closest_reference_id = reference_ids[0]
    closest_band_rows = [
        {
            "reference_id": closest_reference_id,
            "locus_id": locus_id,
            "product_size_bp": row_lookup[(locus_id, closest_reference_id)][
                "product_size_bp"
            ],
            "repeat_count": row_lookup[(locus_id, closest_reference_id)][
                "repeat_count"
            ],
        }
        for locus_id in locus_ids
    ]
    _write_tsv(
        closest_band_rows,
        closest_reference_bands_path,
        CLOSEST_REFERENCE_BAND_FIELDS,
    )

    query_tree_label = sample_id if sample_id not in reference_ids else query_name
    tree_labels = [*reference_ids, query_tree_label]
    combined_tree_path = output / "combined_markers.tree"
    combined_tree_path.write_text(
        _neighbor_joining_tree_from_matrix(
            tree_labels, np.zeros((len(tree_labels), len(tree_labels)), dtype=float)
        )
    )
    if progress is not None:
        tie_detail = (
            "; resolved by canonical genome identity and whole-genome dnadiff"
            if tie_break_status == "APPLIED"
            else ""
        )
        progress.step(
            f"Resolved {len(reference_ids):,} tied {match_type.lower().replace('_', ' ')} "
            f"reference(s) from the sequence index{tie_detail}; skipped MAFFT and EPA-ng"
        )
    paths = {
        "phylogeny": output,
        "phylogenetic_distances": detail_path,
        "phylogenetic_matches": summary_path,
        "phylogenetic_status": status_path,
        "marker_components": marker_components_path,
        "locus_marker_distances": locus_marker_path,
        "combined_marker_matches": combined_path,
        "combined_marker_tree": combined_tree_path,
        "closest_reference_bands": closest_reference_bands_path,
    }
    if dnadiff_result_path is not None and dnadiff_result_path.exists():
        paths["whole_genome_dnadiff"] = dnadiff_result_path
    return paths


def _add_taxon_assignment_outputs(
    paths: dict[str, Path],
    *,
    sample_id: str,
    target_taxon_id: str | None,
    calibration_path: str | Path | None,
    reference_metadata: dict[str, dict[str, str]],
    panel_sha256: str,
    database_signature: str,
    expected_loci: int,
    alpha: float | None,
    min_loci: int | None,
    min_locus_fraction: float,
    bootstrap_replicates: int,
    min_bootstrap_support: float,
    max_mean_placement_entropy: float | None,
    min_median_placement_lwr: float | None,
    taxon_identification: bool | None = None,
    identification_k: int = 3,
    identification_minimum_margin: float = 0.1,
) -> dict[str, Path]:
    reference_taxa, _taxon_names = _reference_taxa(reference_metadata)
    automatic_enabled = taxon_identification is not False and bool(reference_taxa)
    if automatic_enabled:
        assignment = assign_best_taxon(
            sample_id=sample_id,
            locus_marker_rows=_read_tsv_dicts(paths["locus_marker_distances"]),
            reference_metadata=reference_metadata,
            expected_loci=expected_loci,
            snp_weight=_read_marker_weight(paths.get("combined_marker_matches"), "snp_weight"),
            repeat_weight=_read_marker_weight(paths.get("combined_marker_matches"), "repeat_weight"),
            k=identification_k,
            minimum_loci=min_loci or 3,
            minimum_locus_fraction=min_locus_fraction,
            minimum_relative_margin=identification_minimum_margin,
        )
        output = paths["phylogeny"]
        summary_path = output / "taxonomic_identification.tsv"
        evidence_path = output / "taxonomic_identification_evidence.tsv"
        _write_tsv([assignment.summary], summary_path, TAXONOMIC_SUMMARY_FIELDS)
        _write_tsv(list(assignment.evidence), evidence_path, TAXONOMIC_EVIDENCE_FIELDS)
        paths = {
            **paths,
            "taxonomic_identification": summary_path,
            "taxonomic_identification_evidence": evidence_path,
        }
    if target_taxon_id is None and calibration_path is None:
        return paths
    if not target_taxon_id or calibration_path is None:
        raise ValueError(
            "Taxon assignment requires both target_taxon_id and taxon_calibration_path"
        )
    calibration = TaxonCalibration.read(calibration_path)
    if calibration.panel_sha256 and calibration.panel_sha256 != panel_sha256:
        raise ValueError(
            "Taxon calibration panel signature does not match the active MLVA panel"
        )
    if (
        calibration.database_signature
        and calibration.database_signature != database_signature
    ):
        raise ValueError(
            "Taxon calibration database signature does not match the active "
            "reference sequence database"
        )
    marker_rows = _read_tsv_dicts(paths["locus_marker_distances"])
    placement_rows = _read_tsv_dicts(paths["phylogenetic_distances"])
    assignment = assign_target_taxon(
        sample_id=sample_id,
        target_taxon_id=str(target_taxon_id),
        locus_marker_rows=marker_rows,
        placement_rows=placement_rows,
        reference_metadata=reference_metadata,
        calibration=calibration,
        alpha=alpha,
        min_loci=min_loci,
        min_locus_fraction=min_locus_fraction,
        bootstrap_replicates=bootstrap_replicates,
        min_bootstrap_support=min_bootstrap_support,
        max_mean_placement_entropy=max_mean_placement_entropy,
        min_median_placement_lwr=min_median_placement_lwr,
        expected_loci=expected_loci,
    )
    output = paths["phylogeny"]
    summary_path = output / "taxon_assignment.tsv"
    candidates_path = output / "taxon_assignment_candidates.tsv"
    loci_path = output / "taxon_assignment_loci.tsv"
    _write_tsv([assignment.summary], summary_path, TAXON_ASSIGNMENT_FIELDS)
    _write_tsv(list(assignment.candidates), candidates_path, TAXON_CANDIDATE_FIELDS)
    _write_tsv(list(assignment.loci), loci_path, TAXON_LOCUS_FIELDS)
    return {
        **paths,
        "taxon_assignment": summary_path,
        "taxon_assignment_candidates": candidates_path,
        "taxon_assignment_loci": loci_path,
    }


def _read_marker_weight(path: str | Path | None, field: str) -> float:
    if path is None:
        return 1.0
    rows = _read_tsv_dicts(path)
    value = _finite(rows[0].get(field)) if rows else None
    return 1.0 if value is None else value


def _read_tsv_dicts(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_phylogenetic_placement(
    query_sequences: dict[str, str],
    database_path: str | Path,
    outdir: str | Path,
    sample_id: str,
    locus_ids: set[str] | list[Locus],
    threads: int,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    epa_ng_bin: str = "epa-ng",
    raxml_model: str = "DNA",
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    reference_metadata_path: str | Path | None = None,
    progress: ProgressReporter | None = None,
    exact_match_fast_path: bool = True,
    query_assembly_path: str | Path | None = None,
    dnadiff_bin: str = "dnadiff",
    target_taxon_id: str | None = None,
    taxon_calibration_path: str | Path | None = None,
    taxon_alpha: float | None = None,
    taxon_min_loci: int | None = None,
    taxon_min_locus_fraction: float = 0.8,
    taxon_bootstrap_replicates: int = 2000,
    taxon_min_bootstrap_support: float = 0.95,
    taxon_max_mean_placement_entropy: float | None = None,
    taxon_min_median_placement_lwr: float | None = None,
    taxon_identification: bool | None = None,
    taxon_k: int = 3,
    taxon_minimum_margin: float = 0.1,
) -> dict[str, Path]:
    if snp_weight < 0 or repeat_weight < 0 or snp_weight + repeat_weight <= 0:
        raise ValueError("SNP and repeat weights must be non-negative with a positive total")
    locus_by_id = (
        {locus.locus_id: locus for locus in locus_ids}
        if not isinstance(locus_ids, set)
        else {}
    )
    requested_locus_ids = set(locus_by_id) if locus_by_id else set(locus_ids)
    sequence_database_path, reusable_phylogeny = _resolve_reference_database_layout(
        database_path
    )
    references: dict[str, list[tuple[str, str]]] | None = None
    if reference_metadata_path is None and sequence_database_path.is_dir():
        automatic_metadata = sequence_database_path / "reference_metadata.tsv"
        if automatic_metadata.exists():
            reference_metadata_path = automatic_metadata
    reference_metadata = read_reference_metadata(reference_metadata_path)
    if taxon_identification is True and not _reference_taxa(reference_metadata)[0]:
        raise ValueError(
            "--taxon-identification requires reference metadata with taxon_id "
            "(taxid is normalized during reference construction)"
        )
    automatic_taxon_identification = (
        taxon_identification is not False and bool(_reference_taxa(reference_metadata)[0])
    )
    if (
        locus_by_id
        and exact_match_fast_path
        and target_taxon_id is None
        and not automatic_taxon_identification
    ):
        sequence_index_path: Path | None = (
            sequence_database_path / "reference_sequence_index.tsv"
            if sequence_database_path.is_dir()
            else None
        )
        index_rows = (
            _read_reference_sequence_index(sequence_index_path)
            if sequence_index_path is not None
            else []
        )
        if not index_rows:
            references = read_sequence_database(
                sequence_database_path, requested_locus_ids
            )
            index_rows = reference_sequence_index_rows(
                references,
                list(locus_by_id.values()),
                sequence_database_path,
            )
        elif (
            index_rows[0].get("index_version") != "1"
            or index_rows[0].get("panel_sha256")
            != _panel_digest(list(locus_by_id.values()))
            or index_rows[0].get("database_signature")
            != _database_stat_digest(sequence_database_path)
        ):
            references = read_sequence_database(
                sequence_database_path, requested_locus_ids
            )
            index_rows = reference_sequence_index_rows(
                references,
                list(locus_by_id.values()),
                sequence_database_path,
            )
        match_type, exact_references, exact_loci, exact_query_components = (
            _exact_reference_group(
                query_sequences,
                locus_by_id,
                index_rows,
            )
        )
        if match_type:
            whole_genome_metrics: dict[
                str, dict[str, float | int | bool]
            ] = {}
            dnadiff_result_path: Path | None = None
            dnadiff_available = False
            if query_assembly_path is not None and len(exact_references) > 1:
                mapping_path = sequence_database_path / "reference_assemblies.tsv"
                dnadiff_available = mapping_path.exists()
                if dnadiff_available:
                    whole_genome_metrics, dnadiff_result_path = _dnadiff_tie_break(
                        query_assembly_path,
                        exact_references,
                        sequence_database_path,
                        Path(outdir) / "phylogeny",
                        threads,
                        dnadiff_bin,
                    )
            exact_paths = _write_exact_match_outputs(
                outdir,
                sample_id,
                match_type,
                exact_references,
                exact_loci,
                locus_by_id,
                exact_query_components,
                index_rows,
                reference_metadata,
                snp_weight,
                repeat_weight,
                progress,
                whole_genome_metrics,
                dnadiff_result_path,
                dnadiff_available,
                query_assembly_path is not None,
            )
            return _add_taxon_assignment_outputs(
                exact_paths,
                sample_id=sample_id,
                target_taxon_id=target_taxon_id,
                calibration_path=taxon_calibration_path,
                reference_metadata=reference_metadata,
                panel_sha256=_panel_digest(list(locus_by_id.values())),
                database_signature=_database_stat_digest(sequence_database_path),
                expected_loci=len(requested_locus_ids),
                alpha=taxon_alpha,
                min_loci=taxon_min_loci,
                min_locus_fraction=taxon_min_locus_fraction,
                bootstrap_replicates=taxon_bootstrap_replicates,
                min_bootstrap_support=taxon_min_bootstrap_support,
                max_mean_placement_entropy=taxon_max_mean_placement_entropy,
                min_median_placement_lwr=taxon_min_median_placement_lwr,
                taxon_identification=taxon_identification,
                identification_k=taxon_k,
                identification_minimum_margin=taxon_minimum_margin,
            )
    if references is None:
        references = read_sequence_database(sequence_database_path, requested_locus_ids)
    mafft = check_mafft(mafft_bin)
    epa_ng = check_epa_ng(epa_ng_bin)
    raxml_ng: str | None = None
    output = Path(outdir) / "phylogeny"
    output.mkdir(parents=True, exist_ok=True)
    detail_rows: list[dict] = []
    status_rows: list[dict] = []
    safe_sample = _SAFE_FILE.sub("_", sample_id).strip("_") or "sample"
    query_name = f"QUERY__{safe_sample}"
    placed_loci: list[str] = []
    marker_rows: list[dict] = []
    query_components_by_locus: dict[str, MarkerComponents] = {}
    reference_components_by_locus: dict[str, dict[str, MarkerComponents]] = {}
    reference_distance_matrices: dict[str, tuple[list[str], np.ndarray]] = {}
    direct_snp_distances_by_locus: dict[str, dict[str, float]] = {}
    placement_members_by_locus: dict[str, dict[str, list[str]]] = {}
    reference_tree_jobs: dict[str, _ReferenceTreeJob] = {}
    placement_jobs: dict[str, _PlacementJob] = {}
    sequence_artifacts: list[Path] = []

    for locus_id in sorted(references):
        query_sequence = query_sequences.get(locus_id, "")
        locus = locus_by_id.get(locus_id)
        locus_references = references[locus_id]
        if locus is not None:
            reference_components = {
                reference_id: decompose_marker_sequence(locus, sequence)
                for reference_id, sequence in locus_references
            }
            reference_components_by_locus[locus_id] = reference_components
            for reference_id, components in reference_components.items():
                marker_rows.append(
                    _marker_component_row(
                        sample_id, locus, "reference", reference_id, components
                    )
                )
            if query_sequence:
                query_components = decompose_marker_sequence(locus, query_sequence)
                query_components_by_locus[locus_id] = query_components
                marker_rows.append(
                    _marker_component_row(
                        sample_id, locus, "query", query_name, query_components
                    )
                )
                query_sequence = query_components.snp_sequence
            haplotype_members: dict[str, list[str]] = {}
            for reference_id, _sequence in locus_references:
                haplotype_members.setdefault(
                    reference_components[reference_id].snp_sequence, []
                ).append(reference_id)
            if len(haplotype_members) >= 2:
                placement_members = {
                    sorted(members)[0]: sorted(members)
                    for members in haplotype_members.values()
                }
                locus_references = [
                    (
                        representative,
                        reference_components[representative].snp_sequence,
                    )
                    for representative in sorted(placement_members)
                ]
            else:
                placement_members = {
                    reference_id: [reference_id]
                    for reference_id, _sequence in locus_references
                }
                locus_references = [
                    (reference_id, reference_components[reference_id].snp_sequence)
                    for reference_id, _sequence in locus_references
                ]
            placement_members_by_locus[locus_id] = placement_members
        else:
            placement_members_by_locus[locus_id] = {
                reference_id: [reference_id]
                for reference_id, _sequence in locus_references
            }
        if query_sequence and any(reference_id == query_name for reference_id, _sequence in references[locus_id]):
            raise ValueError(
                f"Reference id {query_name!r} is reserved for the query sequence"
            )
        safe_locus = _SAFE_FILE.sub("_", locus_id).strip("_") or "locus"
        locus_dir = output / safe_locus
        locus_dir.mkdir(parents=True, exist_ok=True)
        reference_fasta = locus_dir / "references.fasta"
        reference_alignment = locus_dir / "references.aligned.fasta"
        reference_tree = locus_dir / "reference_tree.nwk"
        reference_prefix = locus_dir / "reference"
        reference_model = Path(f"{reference_prefix}.raxml.bestModel")
        query_fasta = locus_dir / "query.fasta"
        placed_alignment = locus_dir / "query_placed.aligned.fasta"
        query_alignment = locus_dir / "query.aligned.fasta"
        sequence_artifacts.extend(
            [
                reference_fasta,
                reference_alignment,
                query_fasta,
                placed_alignment,
                query_alignment,
            ]
        )
        epa_outdir = locus_dir / "epa-ng"
        _write_fasta(locus_references, reference_fasta)
        if reusable_phylogeny is not None:
            reusable_locus_dir = reusable_phylogeny / safe_locus
            reusable_alignment = reusable_locus_dir / "references.aligned.fasta.gz"
            if not reusable_alignment.is_file():
                reusable_alignment = reusable_locus_dir / "references.aligned.fasta"
            reusable_tree = reusable_locus_dir / "reference_tree.nwk"
            reusable_model = reusable_locus_dir / "reference.raxml.bestModel"
            missing_artifacts = [
                path
                for path in (reusable_alignment, reusable_tree, reusable_model)
                if not path.is_file()
            ]
            if missing_artifacts:
                status_rows.append(
                    {
                        "locus_id": locus_id,
                        "reference_sequences": len(references[locus_id]),
                        "query_sequence": "yes" if query_sequence else "no",
                        "status": "REFERENCE_TREE_UNAVAILABLE",
                    }
                )
                continue
            _copy_reference_artifact(reusable_alignment, reference_alignment)
            _copy_reference_artifact(reusable_tree, reference_tree)
            _copy_reference_artifact(reusable_model, reference_model)
            reusable_names = {name for name, _sequence in _read_fasta(reference_alignment)}
            placement_members = placement_members_by_locus[locus_id]
            if reusable_names == {
                reference_id for reference_id, _sequence in references[locus_id]
            }:
                placement_members_by_locus[locus_id] = {
                    reference_id: [reference_id] for reference_id in reusable_names
                }
            elif reusable_names != set(placement_members):
                raise RuntimeError(
                    f"Reusable reference alignment for {locus_id} has tip identifiers "
                    "that do not match either the raw references or collapsed SNP haplotypes"
                )
        else:
            if raxml_ng is None:
                raxml_ng = check_raxml_ng(raxml_ng_bin)
            reference_tree_jobs[locus_id] = _ReferenceTreeJob(
                locus_id=locus_id,
                reference_fasta=reference_fasta,
                reference_alignment=reference_alignment,
                reference_prefix=reference_prefix,
                reference_tree=reference_tree,
                reference_model=reference_model,
            )
        if not query_sequence:
            status_rows.append(
                {
                    "locus_id": locus_id,
                    "reference_sequences": len(references[locus_id]),
                    "query_sequence": "no",
                    "status": "REFERENCE_TREE_ONLY",
                }
            )
            continue
        _write_fasta([(query_name, query_sequence)], query_fasta)
        placement_jobs[locus_id] = _PlacementJob(
            locus_id=locus_id,
            query_name=query_name,
            query_fasta=query_fasta,
            reference_alignment=reference_alignment,
            reference_tree=reference_tree,
            reference_model=reference_model,
            placed_alignment=placed_alignment,
            query_alignment=query_alignment,
            epa_outdir=epa_outdir,
        )

    if reference_tree_jobs:
        cpu_budget = resolve_threads(threads)
        worker_count = min(cpu_budget, len(reference_tree_jobs))
        native_threads = max(1, cpu_budget // worker_count)
        if progress is not None:
            progress.step(
                f"Building {len(reference_tree_jobs):,} independent reference "
                f"locus trees with {worker_count} worker(s) and "
                f"{native_threads} native thread(s) per worker"
            )
        if worker_count == 1:
            for locus_id, job in reference_tree_jobs.items():
                _run_reference_tree_job(
                    job, mafft, str(raxml_ng), raxml_model, native_threads
                )
        else:
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(
                        _run_reference_tree_job,
                        job,
                        mafft,
                        str(raxml_ng),
                        raxml_model,
                        native_threads,
                    ): locus_id
                    for locus_id, job in reference_tree_jobs.items()
                }
                for future in as_completed(futures):
                    locus_id = futures[future]
                    try:
                        future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            f"Reference tree construction failed for locus "
                            f"{locus_id!r}"
                        ) from exc

    placement_results: dict[str, Path] = {}
    if placement_jobs:
        cpu_budget = resolve_threads(threads)
        worker_count = min(cpu_budget, len(placement_jobs))
        native_threads = max(1, cpu_budget // worker_count)
        if progress is not None:
            progress.step(
                f"Running {len(placement_jobs):,} independent EPA-ng locus "
                f"placements with {worker_count} worker(s) and "
                f"{native_threads} native thread(s) per worker"
            )
        if worker_count == 1:
            for completed, (locus_id, job) in enumerate(
                placement_jobs.items(), start=1
            ):
                placement_results[locus_id] = _run_placement_job(
                    job, mafft, epa_ng, native_threads
                )
                if progress is not None:
                    progress.count(
                        "Completed EPA-ng loci",
                        completed,
                        len(placement_jobs),
                        force=True,
                    )
        else:
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(
                        _run_placement_job, job, mafft, epa_ng, native_threads
                    ): locus_id
                    for locus_id, job in placement_jobs.items()
                }
                for completed, future in enumerate(as_completed(futures), start=1):
                    locus_id = futures[future]
                    try:
                        placement_results[locus_id] = future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            f"Phylogenetic placement failed for locus {locus_id!r}"
                        ) from exc
                    if progress is not None:
                        progress.count(
                            "Completed EPA-ng loci",
                            completed,
                            len(placement_jobs),
                            force=True,
                        )

    # Parse completed placements in locus order for deterministic tables.
    if placement_results and progress is not None:
        progress.step("Computing reference-tip distances and placement summaries")
    for locus_id in sorted(placement_results):
        jplace_path = placement_results[locus_id]
        (
            placement_distances,
            expected_placement_distances,
            placement_metadata,
        ) = read_epa_ng_placement_statistics(jplace_path)
        with Path(jplace_path).open() as jplace_handle:
            jplace_tree = _parse_newick(str(json.load(jplace_handle)["tree"]))
        tree_names, tree_distance_matrix = _tip_patristic_distance_matrix(jplace_tree)
        placement_members = placement_members_by_locus[locus_id]
        expected_tree_references = set(placement_members)
        if set(placement_distances) != expected_tree_references:
            missing = sorted(expected_tree_references - set(placement_distances))
            unexpected = sorted(set(placement_distances) - expected_tree_references)
            raise RuntimeError(
                f"EPA-ng tree/reference mismatch for {locus_id}: missing={missing}, unexpected={unexpected}"
            )
        representative_by_reference = {
            reference_id: representative
            for representative, members in placement_members.items()
            for reference_id in members
        }
        raw_reference_names = sorted(representative_by_reference)
        tree_indexes = {name: index for index, name in enumerate(tree_names)}
        representative_indexes = np.asarray(
            [tree_indexes[representative_by_reference[name]] for name in raw_reference_names],
            dtype=np.intp,
        )
        reference_distance_matrices[locus_id] = (
            raw_reference_names,
            tree_distance_matrix[
                np.ix_(representative_indexes, representative_indexes)
            ],
        )
        aligned_records = dict(_read_fasta(placement_jobs[locus_id].placed_alignment))
        aligned_query = aligned_records.get(query_name)
        if aligned_query is None:
            raise RuntimeError(
                f"MAFFT placement alignment for {locus_id} lacks query {query_name!r}"
            )
        missing_aligned_references = expected_tree_references - set(aligned_records)
        if missing_aligned_references:
            raise RuntimeError(
                f"MAFFT placement alignment for {locus_id} lacks references: "
                + ", ".join(sorted(missing_aligned_references))
            )
        direct_snp_distances_by_locus[locus_id] = {
            reference_id: _aligned_snp_distance(
                aligned_query,
                aligned_records[representative_by_reference[reference_id]],
            )
            for reference_id in raw_reference_names
        }
        for representative, distance in placement_distances.items():
            for reference_id in placement_members[representative]:
                detail_rows.append(
                    {
                        "sample_id": sample_id,
                        "locus_id": locus_id,
                        "reference_id": reference_id,
                        "phylogenetic_distance": f"{distance:.8f}",
                        "likelihood_weighted_phylogenetic_distance": (
                            f"{expected_placement_distances[representative]:.8f}"
                        ),
                        "placement_edge": placement_metadata["placement_edge"],
                        "like_weight_ratio": f'{placement_metadata["like_weight_ratio"]:.8f}',
                        "pendant_length": f'{placement_metadata["pendant_length"]:.8f}',
                        "distal_length": f'{placement_metadata["distal_length"]:.8f}',
                        "placement_count": placement_metadata["placement_count"],
                        "placement_entropy": f'{placement_metadata["placement_entropy"]:.8f}',
                        "reference_tree_scale": f'{placement_metadata["reference_tree_scale"]:.8f}',
                    }
                )
        placed_loci.append(locus_id)
        status_rows.append(
            {
                "locus_id": locus_id,
                "reference_sequences": len(references[locus_id]),
                "query_sequence": "yes",
                "status": "PLACED",
            }
        )
        if progress is not None:
            progress.count(
                "Processed placement summaries",
                len(placed_loci),
                len(placement_results),
                force=True,
            )

    detail_path = output / "locus_phylogenetic_distances.tsv"
    _write_tsv(detail_rows, detail_path, PLACEMENT_FIELDS)
    marker_components_path = output / "marker_components.tsv"
    _write_tsv(marker_rows, marker_components_path, MARKER_COMPONENT_FIELDS)

    repeat_scale_by_locus: dict[str, float] = {}
    for locus_id, components_by_reference in reference_components_by_locus.items():
        counts = [
            float(components.repeat_count_raw)
            for components in components_by_reference.values()
            if components.repeat_count_raw is not None
        ]
        repeat_scale_by_locus[locus_id] = max(
            0.5, statistics.pstdev(counts) if len(counts) > 1 else 0.0
        )

    direct_snp_scale_by_locus: dict[str, float] = {}
    for locus_id, distances_by_reference in direct_snp_distances_by_locus.items():
        positive = [distance for distance in distances_by_reference.values() if distance > 0]
        direct_snp_scale_by_locus[locus_id] = (
            statistics.median(positive) if positive else 1.0
        )

    locus_marker_rows: list[dict] = []
    for row in detail_rows:
        locus_id = str(row["locus_id"])
        reference_id = str(row["reference_id"])
        query_components = query_components_by_locus.get(locus_id)
        reference_components = reference_components_by_locus.get(locus_id, {}).get(
            reference_id
        )
        query_count = (
            query_components.repeat_count_raw if query_components is not None else None
        )
        reference_count = (
            reference_components.repeat_count_raw
            if reference_components is not None
            else None
        )
        if query_count is not None and reference_count is not None:
            repeat_delta: str | float = abs(query_count - reference_count)
            normalized_repeat: str | float = (
                repeat_delta / repeat_scale_by_locus[locus_id]
            )
        else:
            repeat_delta = ""
            normalized_repeat = ""
        weighted_snp = float(row["likelihood_weighted_phylogenetic_distance"])
        tree_scale = max(float(row["reference_tree_scale"]), 1e-12)
        placement_normalized_snp = weighted_snp / tree_scale
        direct_snp = direct_snp_distances_by_locus[locus_id][reference_id]
        direct_snp_scale = max(direct_snp_scale_by_locus[locus_id], 1e-12)
        normalized_direct_snp = direct_snp / direct_snp_scale
        exact_snp_match = (
            query_components is not None
            and reference_components is not None
            and query_components.snp_sequence == reference_components.snp_sequence
        )
        if query_components is None or reference_components is None:
            exact_snp_match = direct_snp == 0
        normalized_snp = (
            0.0
            if exact_snp_match
            else (placement_normalized_snp + normalized_direct_snp) / 2.0
        )
        locus_marker_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "reference_id": reference_id,
                "query_repeat_count": "" if query_count is None else f"{query_count:.6f}",
                "reference_repeat_count": ""
                if reference_count is None
                else f"{reference_count:.6f}",
                "repeat_count_delta": ""
                if repeat_delta == ""
                else f"{float(repeat_delta):.6f}",
                "normalized_repeat_distance": ""
                if normalized_repeat == ""
                else f"{float(normalized_repeat):.8f}",
                "likelihood_weighted_snp_distance": f"{weighted_snp:.8f}",
                "reference_tree_scale": f"{tree_scale:.8f}",
                "placement_normalized_snp_distance": f"{placement_normalized_snp:.8f}",
                "direct_snp_distance": f"{direct_snp:.8f}",
                "direct_snp_scale": f"{direct_snp_scale:.8f}",
                "normalized_direct_snp_distance": f"{normalized_direct_snp:.8f}",
                "exact_snp_match": "yes" if exact_snp_match else "no",
                "normalized_snp_distance": f"{normalized_snp:.8f}",
                "placement_entropy": row["placement_entropy"],
            }
        )
    locus_marker_path = output / "locus_marker_distances.tsv"
    _write_tsv(
        locus_marker_rows, locus_marker_path, LOCUS_MARKER_DISTANCE_FIELDS
    )
    totals: dict[str, list[float]] = {}
    weighted_totals: dict[str, list[float]] = {}
    for row in detail_rows:
        reference_id = str(row["reference_id"])
        totals.setdefault(reference_id, []).append(float(row["phylogenetic_distance"]))
        weighted_totals.setdefault(reference_id, []).append(
            float(row["likelihood_weighted_phylogenetic_distance"])
        )
    # A partial reference must not win merely because absent loci contribute
    # nothing. Rank only references represented in every successfully placed
    # locus, so the total is comparable across candidates.
    complete_totals = {
        reference_id: values
        for reference_id, values in totals.items()
        if len(values) == len(placed_loci)
    }
    summary_rows = []
    ordered = sorted(
        complete_totals.items(),
        key=lambda item: (
            sum(weighted_totals[item[0]]),
            sum(item[1]),
            item[0],
        ),
    )
    for rank, (reference_id, values) in enumerate(ordered, start=1):
        weighted_values = weighted_totals[reference_id]
        weighted_total = sum(weighted_values)
        if rank < len(ordered):
            next_reference = ordered[rank][0]
            next_total = sum(weighted_totals[next_reference])
            distance_gap: str | float = f"{next_total - weighted_total:.8f}"
            relative_gap: str | float = f"{(next_total - weighted_total) / max(next_total, 1e-12):.8f}"
        else:
            distance_gap = ""
            relative_gap = ""
        summary_rows.append(
            {
                "sample_id": sample_id,
                "reference_id": reference_id,
                "total_phylogenetic_distance": f"{sum(values):.8f}",
                "total_likelihood_weighted_distance": f"{weighted_total:.8f}",
                "compared_loci": len(values),
                "mean_phylogenetic_distance": f"{sum(values) / len(values):.8f}",
                "mean_likelihood_weighted_distance": (
                    f"{weighted_total / len(weighted_values):.8f}"
                ),
                "distance_gap_to_next": distance_gap,
                "relative_distance_gap_to_next": relative_gap,
                "rank": rank,
            }
        )
    summary_path = output / "phylogenetic_matches.tsv"
    _write_tsv(summary_rows, summary_path, SUMMARY_FIELDS)

    complete_references = set(complete_totals)
    combined_totals: dict[str, dict[str, float | int]] = {}
    for row in locus_marker_rows:
        reference_id = str(row["reference_id"])
        if reference_id not in complete_references:
            continue
        values = combined_totals.setdefault(
            reference_id,
            {
                "snp": 0.0,
                "placement_normalized_snp": 0.0,
                "normalized_direct_snp": 0.0,
                "normalized_snp": 0.0,
                "repeat": 0.0,
                "normalized_repeat": 0.0,
                "loci": 0,
                "repeat_loci": 0,
                "exact_snp_loci": 0,
                "exact_marker_loci": 0,
            },
        )
        values["snp"] += float(row["likelihood_weighted_snp_distance"])
        values["placement_normalized_snp"] += float(
            row["placement_normalized_snp_distance"]
        )
        values["normalized_direct_snp"] += float(
            row["normalized_direct_snp_distance"]
        )
        values["normalized_snp"] += float(row["normalized_snp_distance"])
        values["loci"] += 1
        if row["exact_snp_match"] == "yes":
            values["exact_snp_loci"] += 1
        if row["repeat_count_delta"] != "":
            values["repeat"] += float(row["repeat_count_delta"])
            values["normalized_repeat"] += float(row["normalized_repeat_distance"])
            values["repeat_loci"] += 1
        repeat_matches = row["repeat_count_delta"] == "" or math.isclose(
            float(row["repeat_count_delta"]), 0.0, abs_tol=1e-12
        )
        if row["exact_snp_match"] == "yes" and repeat_matches:
            values["exact_marker_loci"] += 1

    required_repeat_loci = {
        locus_id
        for locus_id in placed_loci
        if query_components_by_locus.get(locus_id) is not None
        and query_components_by_locus[locus_id].repeat_count_raw is not None
    }
    eligible_combined_references = [
        reference_id
        for reference_id, values in combined_totals.items()
        if int(values["repeat_loci"]) == len(required_repeat_loci)
    ]
    combined_order = sorted(
        eligible_combined_references,
        key=lambda reference_id: (
            snp_weight * float(combined_totals[reference_id]["normalized_snp"])
            + repeat_weight
            * float(combined_totals[reference_id]["normalized_repeat"]),
            reference_id,
        ),
    )
    combined_rows: list[dict] = []
    combined_scores = {
        reference_id: (
            snp_weight * float(combined_totals[reference_id]["normalized_snp"])
            + repeat_weight
            * float(combined_totals[reference_id]["normalized_repeat"])
        )
        for reference_id in combined_order
    }
    placement_only_scores = {
        reference_id: (
            snp_weight
            * float(combined_totals[reference_id]["placement_normalized_snp"])
            + repeat_weight
            * float(combined_totals[reference_id]["normalized_repeat"])
        )
        for reference_id in combined_order
    }
    minimum_placement_only_score = min(placement_only_scores.values(), default=0.0)
    previous_score: float | None = None
    previous_rank = 0
    for index, reference_id in enumerate(combined_order):
        values = combined_totals[reference_id]
        score = combined_scores[reference_id]
        rank = (
            previous_rank
            if previous_score is not None
            and math.isclose(score, previous_score, rel_tol=1e-12, abs_tol=1e-12)
            else index + 1
        )
        previous_score = score
        previous_rank = rank
        if index + 1 < len(combined_order):
            next_score = combined_scores[combined_order[index + 1]]
            gap = next_score - score
            relative_gap: str | float = gap / max(next_score, 1e-12)
        else:
            gap = ""
            relative_gap = ""
        exact_marker_match = int(values["exact_marker_loci"]) == int(values["loci"])
        ranking_warning = (
            "EXACT_MATCH_OVERRIDES_PLACEMENT"
            if exact_marker_match
            and placement_only_scores[reference_id]
            > minimum_placement_only_score + 1e-12
            else ""
        )
        combined_rows.append(
            {
                "sample_id": sample_id,
                "reference_id": reference_id,
                "total_likelihood_weighted_snp_distance": f'{float(values["snp"]):.8f}',
                "total_placement_normalized_snp_distance": f'{float(values["placement_normalized_snp"]):.8f}',
                "total_normalized_direct_snp_distance": f'{float(values["normalized_direct_snp"]):.8f}',
                "total_normalized_snp_distance": f'{float(values["normalized_snp"]):.8f}',
                "total_repeat_count_distance": f'{float(values["repeat"]):.8f}',
                "total_normalized_repeat_distance": f'{float(values["normalized_repeat"]):.8f}',
                "snp_weight": f"{snp_weight:.6f}",
                "repeat_weight": f"{repeat_weight:.6f}",
                "combined_marker_distance": f"{score:.8f}",
                "compared_loci": values["loci"],
                "repeat_compared_loci": values["repeat_loci"],
                "exact_snp_loci": values["exact_snp_loci"],
                "exact_marker_loci": values["exact_marker_loci"],
                "match_status": "EXACT_MARKER_MATCH"
                if exact_marker_match
                else "NON_EXACT",
                "ranking_warning": ranking_warning,
                "distance_gap_to_next": "" if gap == "" else f"{float(gap):.8f}",
                "relative_distance_gap_to_next": ""
                if relative_gap == ""
                else f"{float(relative_gap):.8f}",
                "collection_date": reference_metadata.get(reference_id, {}).get(
                    "collection_date", ""
                ),
                "latitude": reference_metadata.get(reference_id, {}).get("latitude", ""),
                "longitude": reference_metadata.get(reference_id, {}).get("longitude", ""),
                "location": reference_metadata.get(reference_id, {}).get("location", ""),
                "source": reference_metadata.get(reference_id, {}).get("source", ""),
                "rank": rank,
            }
        )
    combined_path = output / "combined_marker_matches.tsv"
    _write_tsv(combined_rows, combined_path, COMBINED_MARKER_FIELDS)
    closest_reference_bands_path = output / "closest_reference_bands.tsv"
    closest_reference_id = ""
    if combined_rows:
        closest_reference_id = str(combined_rows[0]["reference_id"])
    elif summary_rows:
        closest_reference_id = str(summary_rows[0]["reference_id"])
    closest_reference_band_rows = []
    if closest_reference_id:
        for locus_id in sorted(references):
            reference_sequence = next(
                (
                    sequence
                    for reference_id, sequence in references[locus_id]
                    if reference_id == closest_reference_id
                ),
                "",
            )
            if not reference_sequence:
                continue
            components = reference_components_by_locus.get(locus_id, {}).get(
                closest_reference_id
            )
            closest_reference_band_rows.append(
                {
                    "reference_id": closest_reference_id,
                    "locus_id": locus_id,
                    "product_size_bp": len(reference_sequence),
                    "repeat_count": ""
                    if components is None or components.repeat_count is None
                    else components.repeat_count,
                }
            )
    _write_tsv(
        closest_reference_band_rows,
        closest_reference_bands_path,
        CLOSEST_REFERENCE_BAND_FIELDS,
    )
    query_tree_label = sample_id if sample_id not in eligible_combined_references else query_name
    tree_labels = [*sorted(eligible_combined_references), query_tree_label]
    query_locus_rows = {
        (str(row["locus_id"]), str(row["reference_id"])): row
        for row in locus_marker_rows
    }
    tree_scale_by_locus = {
        str(row["locus_id"]): max(float(row["reference_tree_scale"]), 1e-12)
        for row in detail_rows
    }

    combined_matrix = np.zeros(
        (len(tree_labels), len(tree_labels)), dtype=np.float64
    )
    reference_labels = tree_labels[:-1]
    if progress is not None:
        progress.step(
            f"Combining marker distances for {len(tree_labels):,} tree labels"
        )
    for locus_number, locus_id in enumerate(placed_loci, start=1):
        locus_names, locus_distances = reference_distance_matrices[locus_id]
        locus_indexes = {name: index for index, name in enumerate(locus_names)}
        try:
            reference_indexes = np.asarray(
                [locus_indexes[label] for label in reference_labels], dtype=np.intp
            )
        except KeyError as exc:
            raise RuntimeError(
                f"Combined-tree reference {exc.args[0]!r} is absent from locus "
                f"{locus_id!r}"
            ) from exc
        if reference_labels:
            reference_snp = locus_distances[
                np.ix_(reference_indexes, reference_indexes)
            ] / tree_scale_by_locus[locus_id]
            combined_matrix[:-1, :-1] += snp_weight * reference_snp
            query_snp = np.asarray(
                [
                    float(
                        query_locus_rows[(locus_id, reference_id)][
                            "normalized_snp_distance"
                        ]
                    )
                    for reference_id in reference_labels
                ],
                dtype=np.float64,
            )
            combined_matrix[-1, :-1] += snp_weight * query_snp
            combined_matrix[:-1, -1] += snp_weight * query_snp

        if locus_id in repeat_scale_by_locus:
            reference_components = reference_components_by_locus.get(locus_id, {})
            query_components = query_components_by_locus.get(locus_id)
            repeat_counts = np.asarray(
                [
                    (
                        reference_components[reference_id].repeat_count_raw
                        if reference_id in reference_components
                        else None
                    )
                    for reference_id in reference_labels
                ]
                + [
                    query_components.repeat_count_raw
                    if query_components is not None
                    else None
                ],
                dtype=np.float64,
            )
            valid_repeats = np.isfinite(repeat_counts)
            repeat_distances = np.abs(
                repeat_counts[:, None] - repeat_counts[None, :]
            )
            repeat_distances[
                ~(valid_repeats[:, None] & valid_repeats[None, :])
            ] = 0.0
            combined_matrix += (
                repeat_weight
                * repeat_distances
                / repeat_scale_by_locus[locus_id]
            )
        if progress is not None:
            progress.count(
                "Combined marker loci",
                locus_number,
                len(placed_loci),
                force=locus_number == len(placed_loci),
            )
    np.fill_diagonal(combined_matrix, 0.0)
    combined_tree_path = output / "combined_markers.tree"
    if progress is not None:
        progress.step(
            f"Building combined neighbor-joining tree for {len(tree_labels):,} labels"
        )
    combined_tree_path.write_text(
        _neighbor_joining_tree_from_matrix(tree_labels, combined_matrix)
    )
    if progress is not None:
        progress.step("Finished phylogenetic placement summaries")
    status_path = output / "locus_status.tsv"
    _write_tsv(status_rows, status_path, STATUS_FIELDS)
    for sequence_path in sequence_artifacts:
        if sequence_path.is_file():
            gzip_output_file(sequence_path)
    paths = {
        "phylogeny": output,
        "phylogenetic_distances": detail_path,
        "phylogenetic_matches": summary_path,
        "phylogenetic_status": status_path,
        "marker_components": marker_components_path,
        "locus_marker_distances": locus_marker_path,
        "combined_marker_matches": combined_path,
        "combined_marker_tree": combined_tree_path,
        "closest_reference_bands": closest_reference_bands_path,
    }
    return _add_taxon_assignment_outputs(
        paths,
        sample_id=sample_id,
        target_taxon_id=target_taxon_id,
        calibration_path=taxon_calibration_path,
        reference_metadata=reference_metadata,
        panel_sha256=(
            _panel_digest(list(locus_by_id.values())) if locus_by_id else ""
        ),
        database_signature=_database_stat_digest(sequence_database_path),
        expected_loci=len(requested_locus_ids),
        alpha=taxon_alpha,
        min_loci=taxon_min_loci,
        min_locus_fraction=taxon_min_locus_fraction,
        bootstrap_replicates=taxon_bootstrap_replicates,
        min_bootstrap_support=taxon_min_bootstrap_support,
        max_mean_placement_entropy=taxon_max_mean_placement_entropy,
        min_median_placement_lwr=taxon_min_median_placement_lwr,
        taxon_identification=taxon_identification,
        identification_k=taxon_k,
        identification_minimum_margin=taxon_minimum_margin,
    )


def build_reference_phylogenies(
    database_path: str | Path,
    outdir: str | Path,
    loci: list[Locus],
    threads: int,
    *,
    min_references: int = 3,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    raxml_model: str = "DNA",
    progress: ProgressReporter | None = None,
) -> dict[str, Path]:
    """Build reusable, repeat-masked reference alignments and locus trees.

    The raw amplicons remain in the database.  Trees use the SNP component
    when a repeat can be bounded from panel flanks or motif; this prevents
    variable repeat length from being interpreted as many independent SNPs.
    """
    if min_references < 2:
        raise ValueError("min_references must be at least 2")
    locus_by_id = {locus.locus_id: locus for locus in loci}
    references = read_sequence_database(database_path, set(locus_by_id))
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    database = Path(database_path)
    sequence_index_path = (
        database / "reference_sequence_index.tsv"
        if database.is_dir()
        else output / "reference_sequence_index.tsv"
    )
    _write_tsv(
        reference_sequence_index_rows(references, loci, database),
        sequence_index_path,
        REFERENCE_SEQUENCE_INDEX_FIELDS,
    )
    mafft = check_mafft(mafft_bin)
    raxml_ng = check_raxml_ng(raxml_ng_bin)
    status_rows: list[dict] = []
    component_rows: list[dict] = []
    haplotype_rows: list[dict] = []
    sequence_artifacts: list[Path] = []

    sorted_locus_ids = sorted(locus_by_id)
    for locus_number, locus_id in enumerate(sorted_locus_ids, start=1):
        records = references.get(locus_id, [])
        safe_locus = _SAFE_FILE.sub("_", locus_id).strip("_") or "locus"
        locus_dir = output / safe_locus
        portable_tree = output / f"{safe_locus}.tree"
        portable_tree.unlink(missing_ok=True)
        components_by_reference = {
            reference_id: decompose_marker_sequence(locus_by_id[locus_id], sequence)
            for reference_id, sequence in records
        }
        haplotype_members: dict[str, list[str]] = {}
        for reference_id, components in components_by_reference.items():
            haplotype_members.setdefault(components.snp_sequence, []).append(reference_id)
            component_rows.append(
                _marker_component_row(
                    "REFERENCE_DATABASE",
                    locus_by_id[locus_id],
                    "reference",
                    reference_id,
                    components,
                )
            )
        use_collapsed_haplotypes = len(haplotype_members) >= 2
        masked_records: list[tuple[str, str]] = []
        for snp_sequence, members in sorted(
            haplotype_members.items(), key=lambda item: sorted(item[1])[0]
        ):
            sorted_members = sorted(members)
            representative = sorted_members[0]
            for reference_id in sorted_members:
                haplotype_rows.append(
                    {
                        "locus_id": locus_id,
                        "haplotype_id": representative,
                        "reference_id": reference_id,
                        "snp_sha256": _sequence_digest(snp_sequence),
                    }
                )
            if use_collapsed_haplotypes:
                masked_records.append((representative, snp_sequence))
            else:
                masked_records.extend(
                    (reference_id, snp_sequence) for reference_id in sorted_members
                )
        if len(records) < min_references:
            status_rows.append(
                {
                    "locus_id": locus_id,
                    "reference_sequences": len(records),
                    "tree_haplotypes": len(masked_records),
                    "status": "INSUFFICIENT_REFERENCES",
                    "tree": "",
                }
            )
            if progress is not None:
                progress.count(
                    "Processed tree loci", locus_number, len(sorted_locus_ids), force=True
                )
            continue
        if progress is not None:
            progress.step(
                f"Building tree for {locus_id} ({locus_number}/{len(sorted_locus_ids)}; "
                f"{len(records):,} references, {len(masked_records):,} SNP haplotypes)"
            )
        locus_dir.mkdir(parents=True, exist_ok=True)
        reference_fasta = locus_dir / "references.fasta"
        alignment = locus_dir / "references.aligned.fasta"
        sequence_artifacts.extend([reference_fasta, alignment])
        tree = locus_dir / "reference_tree.nwk"
        prefix = locus_dir / "reference"
        _write_fasta(masked_records, reference_fasta)
        _run_mafft(
            build_mafft_reference_command(reference_fasta, threads, mafft),
            alignment,
            f"reference alignment for {locus_id}",
        )
        _run_raxml_ng(
            build_raxml_ng_command(alignment, prefix, threads, raxml_ng, raxml_model),
            prefix,
            tree,
            f"reference tree search for {locus_id}",
            progress,
        )
        shutil.copyfile(tree, portable_tree)
        status_rows.append(
            {
                "locus_id": locus_id,
                "reference_sequences": len(records),
                "tree_haplotypes": len(masked_records),
                "status": "BUILT",
                "tree": str(portable_tree),
            }
        )
        if progress is not None:
            progress.count(
                "Processed tree loci", locus_number, len(sorted_locus_ids), force=True
            )

    status_path = output / "reference_tree_status.tsv"
    components_path = output / "reference_marker_components.tsv"
    haplotypes_path = output / "reference_haplotype_groups.tsv"
    _write_tsv(
        status_rows,
        status_path,
        ["locus_id", "reference_sequences", "tree_haplotypes", "status", "tree"],
    )
    _write_tsv(component_rows, components_path, MARKER_COMPONENT_FIELDS)
    _write_tsv(haplotype_rows, haplotypes_path, REFERENCE_HAPLOTYPE_FIELDS)
    for sequence_path in sequence_artifacts:
        if sequence_path.is_file():
            gzip_output_file(sequence_path)
    return {
        "phylogeny": output,
        "tree_status": status_path,
        "marker_components": components_path,
        "sequence_index": sequence_index_path,
        "haplotype_groups": haplotypes_path,
    }


def dominant_read_query_sequences(
    features: list[RepeatFeature], asv_rows: list[dict]
) -> dict[str, str]:
    """Return the dominant full amplicon observed at each FASTQ locus."""
    feature_by_key = {(feature.locus_id, feature.read_id): feature for feature in features}
    dominant: dict[str, dict] = {}
    for row in asv_rows:
        locus_id = str(row.get("locus_id", ""))
        current = dominant.get(locus_id)
        if current is None or int(row.get("support_reads") or 0) > int(
            current.get("support_reads") or 0
        ):
            dominant[locus_id] = row
    sequences = {}
    for locus_id, row in dominant.items():
        feature = feature_by_key.get((locus_id, str(row.get("representative_read_id") or "")))
        if feature is not None:
            sequence = (feature.amplicon_sequence or feature.repeat_sequence).upper()
            if sequence:
                sequences[locus_id] = sequence
    return sequences
