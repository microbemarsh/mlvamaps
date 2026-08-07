from __future__ import annotations

import csv
import math
import re
import shutil
import statistics
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .concurrency import resolve_threads
from .io import write_tsv
from .models import Locus
from .phylogeny import (
    _aligned_snp_distance,
    _iupac_find,
    _parse_newick,
    _read_fasta,
    _reverse_complement,
    _run_mafft,
    _run_raxml_ng,
    _tip_patristic_distance_matrix,
    _write_fasta,
    build_mafft_reference_command,
    build_raxml_ng_command,
    check_mafft,
    check_raxml_ng,
    decompose_marker_sequence,
    neighbor_joining_tree_from_matrix,
)


SEQUENCE_STATUS_FIELDS = [
    "sample_id",
    "locus_id",
    "status",
    "source",
    "masking_method",
    "sequence_length",
    "details",
]

LOCUS_STATUS_FIELDS = [
    "locus_id",
    "samples",
    "snp_haplotypes",
    "tree_scale",
    "status",
    "alignment",
    "tree",
]

LOCUS_SNP_FIELDS = [
    "locus_id",
    "sample_1",
    "sample_2",
    "snp_patristic_distance",
    "locus_tree_scale",
    "normalized_snp_distance",
    "repeat_count_distance",
    "locus_repeat_scale",
    "normalized_repeat_distance",
    "combined_locus_distance",
]

COMBINED_PAIRWISE_FIELDS = [
    "sample_1",
    "sample_2",
    "loci_compared",
    "fraction_loci_compared",
    "mean_normalized_snp_distance",
    "mean_normalized_repeat_distance",
    "snp_weight",
    "repeat_weight",
    "combined_marker_distance",
    "comparison_status",
]

COMBINED_OUTPUT_NAMES = (
    "combined_marker_sequence_status.tsv",
    "locus_tree_status.tsv",
    "locus_snp_distances.tsv",
    "combined_marker_pairwise_distances.tsv",
    "combined_marker_distance_matrix.tsv",
    "combined_marker_nj.tree",
    "combined_marker_metadata.tsv",
)

_SAFE_FILE = re.compile(r"[^A-Za-z0-9_.-]+")


@dataclass(frozen=True)
class RecoveredSequence:
    sequence: str
    source: str
    masking_method: str


def _safe_name(value: str) -> str:
    return _SAFE_FILE.sub("_", value).strip("_") or "locus"


def _finite_repeat(value: object) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _read_delimited(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def _precomputed_masked_sequence(sample_dir: Path, locus_id: str) -> str:
    locus_dir = sample_dir / "phylogeny" / _safe_name(locus_id)
    for name in ("query.fasta.gz", "query.fasta"):
        path = locus_dir / name
        if not path.is_file():
            continue
        records = _read_fasta(path)
        if len(records) == 1:
            return records[0][1].replace("-", "").upper()
    return ""


def _retained_sample_sequences(
    sample_dir: Path,
) -> tuple[dict[str, str], dict[str, list[tuple[int, str, str]]]]:
    assembly_sequences: dict[str, str] = {}
    for name in ("assembly_amplicons.fasta.gz", "assembly_amplicons.fasta"):
        path = sample_dir / name
        if path.is_file():
            assembly_sequences = dict(_read_fasta(path))
            break
    candidates: dict[str, list[tuple[int, str, str]]] = {}
    matches = sample_dir / "local_assembly_pcr" / "matches.tsv"
    for row in _read_delimited(matches):
        locus_id = str(row.get("primer_name", ""))
        if locus_id and row.get("product_seq"):
            candidates.setdefault(locus_id, []).append(
                (0, str(matches), str(row["product_seq"]))
            )
    for name in ("local_locus_products.fasta.gz", "local_locus_products.fasta"):
        path = sample_dir / name
        if not path.is_file():
            continue
        for record_id, sequence in _read_fasta(path):
            marker = "_local_primary"
            contig_marker = "_contig_"
            if record_id.endswith(marker):
                locus_id = record_id[: -len(marker)]
            elif contig_marker in record_id:
                locus_id = record_id.split(contig_marker, 1)[0]
            else:
                continue
            candidates.setdefault(locus_id, []).append((1, str(path), sequence))
        break
    return assembly_sequences, candidates


def _masked_candidate(locus: Locus, sequence: str):
    """Crop retained contigs to their primer-bounded product before masking."""
    components = decompose_marker_sequence(locus, sequence)
    oriented = components.oriented_sequence
    reverse_site = _reverse_complement(locus.reverse_primer)
    forward = _iupac_find(locus.forward_primer, oriented)
    reverse = _iupac_find(
        reverse_site,
        oriented,
        forward + len(locus.forward_primer) if forward >= 0 else 0,
    )
    if forward >= 0 and reverse > forward:
        components = decompose_marker_sequence(
            locus, oriented[forward : reverse + len(reverse_site)]
        )
    return components


def recover_masked_sequences(
    samples: list[object],
    locus_order: list[str],
    loci_by_id: dict[str, Locus],
) -> tuple[dict[tuple[str, str], RecoveredSequence], list[dict[str, object]]]:
    """Recover only selected complete products and return repeat-masked SNP sequences."""
    recovered: dict[tuple[str, str], RecoveredSequence] = {}
    status_rows: list[dict[str, object]] = []
    for sample in samples:
        sample_id = str(getattr(sample, "sample_id"))
        sample_dir = Path(getattr(sample, "path")).parent
        assembly_sequences, retained_candidates = _retained_sample_sequences(sample_dir)
        rows_by_locus = {
            str(row.get("locus_id", "")): row for row in getattr(sample, "rows")
        }
        for locus_id in locus_order:
            call = rows_by_locus.get(locus_id, {})
            if _finite_repeat(call.get("repeat_count")) is None:
                status_rows.append(
                    _sequence_status(sample_id, locus_id, "NO_EXACT_REPEAT_CALL")
                )
                continue
            sequence = _precomputed_masked_sequence(sample_dir, locus_id)
            if sequence:
                item = RecoveredSequence(
                    sequence, "precomputed_phylogeny_query", "precomputed_repeat_mask"
                )
                recovered[(sample_id, locus_id)] = item
                status_rows.append(_sequence_status(sample_id, locus_id, "RECOVERED", item))
                continue
            locus = loci_by_id.get(locus_id)
            if locus is None:
                status_rows.append(
                    _sequence_status(
                        sample_id,
                        locus_id,
                        "PANEL_REQUIRED",
                        details=(
                            "no precomputed repeat-masked query was found; provide a rich "
                            "--loci panel to recover and mask retained amplicons"
                        ),
                    )
                )
                continue

            evidence = str(call.get("evidence", ""))
            evidence_sequence = assembly_sequences.get(evidence, "")
            candidates: list[tuple[int, str, str]] = []
            if evidence_sequence:
                candidates.append((0, "assembly_amplicons:evidence", evidence_sequence))
            candidates.extend(retained_candidates.get(locus_id, []))
            accepted: dict[str, tuple[str, str]] = {}
            unmasked = 0
            expected_size = _finite_repeat(call.get("product_size_bp"))
            accepted_priority: int | None = None
            for priority, source, candidate in candidates:
                if accepted_priority is not None and priority > accepted_priority:
                    break
                components = _masked_candidate(locus, candidate)
                if components.masking_method.startswith("unmasked"):
                    unmasked += 1
                    continue
                if (
                    source != "assembly_amplicons:evidence"
                    and expected_size is not None
                    and not math.isclose(
                        expected_size, len(components.oriented_sequence), abs_tol=0.5
                    )
                ):
                    continue
                accepted_priority = priority
                accepted.setdefault(
                    components.snp_sequence,
                    (source, components.masking_method),
                )
                if source == "assembly_amplicons:evidence":
                    accepted = {
                        components.snp_sequence: (source, components.masking_method)
                    }
                    break
            if len(accepted) == 1:
                masked, (source, method) = next(iter(accepted.items()))
                item = RecoveredSequence(masked, source, method)
                recovered[(sample_id, locus_id)] = item
                status_rows.append(_sequence_status(sample_id, locus_id, "RECOVERED", item))
            elif len(accepted) > 1:
                status_rows.append(
                    _sequence_status(
                        sample_id,
                        locus_id,
                        "AMBIGUOUS_RETAINED_PRODUCTS",
                        details=f"{len(accepted)} distinct accepted masked sequences",
                    )
                )
            else:
                details = "no complete selected amplicon sequence was recoverable"
                if unmasked:
                    details += f"; {unmasked} candidate(s) had no safely bounded VNTR tract"
                status_rows.append(
                    _sequence_status(
                        sample_id, locus_id, "SEQUENCE_UNAVAILABLE", details=details
                    )
                )
    return recovered, status_rows


def _sequence_status(
    sample_id: str,
    locus_id: str,
    status: str,
    item: RecoveredSequence | None = None,
    *,
    details: str = "",
) -> dict[str, object]:
    return {
        "sample_id": sample_id,
        "locus_id": locus_id,
        "status": status,
        "source": "" if item is None else item.source,
        "masking_method": "" if item is None else item.masking_method,
        "sequence_length": "" if item is None else len(item.sequence),
        "details": details,
    }


def _write_matrix(path: Path, sample_ids: list[str], matrix: np.ndarray) -> None:
    write_tsv(
        [
            {
                "sample_id": sample_id,
                **{
                    other: f"{float(matrix[row, column]):.8f}"
                    for column, other in enumerate(sample_ids)
                },
            }
            for row, sample_id in enumerate(sample_ids)
        ],
        path,
        ["sample_id", *sample_ids],
    )


def _complete_subset(
    distances: np.ndarray, callable_counts: np.ndarray, sample_ids: list[str]
) -> list[int]:
    active = list(range(len(sample_ids)))
    while len(active) > 1:
        unsupported = ~np.isfinite(distances[np.ix_(active, active)])
        np.fill_diagonal(unsupported, False)
        counts = unsupported.sum(axis=1)
        if int(counts.max(initial=0)) == 0:
            break
        worst = int(counts.max())
        choices = [position for position, count in enumerate(counts) if int(count) == worst]
        remove = min(
            choices,
            key=lambda position: (
                int(callable_counts[active[position]]),
                sample_ids[active[position]].casefold(),
            ),
        )
        active.pop(remove)
    return active


def _copy_metadata_subset(source: Path, destination: Path, retained: set[str]) -> None:
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [row for row in reader if str(row.get("sample_id", "")) in retained]
    write_tsv(rows, destination, fields)


def export_combined_markers(
    samples: list[object],
    locus_order: list[str],
    repeat_calls: np.ndarray,
    metadata_path: Path,
    outdir: str | Path,
    *,
    loci: list[Locus] | None = None,
    min_pairwise_loci: int = 1,
    min_pairwise_fraction: float = 0.0,
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    threads: int = 1,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    raxml_model: str = "DNA",
    force: bool = False,
) -> dict[str, Path | int | str]:
    """Build retrospective per-locus SNP trees and a repeat-aware combined tree."""
    if snp_weight < 0 or repeat_weight < 0 or snp_weight + repeat_weight <= 0:
        raise ValueError("SNP and repeat weights must be non-negative with a positive total")
    output = Path(outdir)
    threads = resolve_threads(threads)
    tree_dir = output / "locus_trees"
    if tree_dir.exists() and force:
        shutil.rmtree(tree_dir)
    tree_dir.mkdir(parents=True, exist_ok=True)
    sample_ids = [str(getattr(sample, "sample_id")) for sample in samples]
    sample_count = len(sample_ids)
    applicable = np.zeros((sample_count, len(locus_order)), dtype=bool)
    locus_indexes = {locus_id: index for index, locus_id in enumerate(locus_order)}
    for sample_index, sample in enumerate(samples):
        for locus_id in getattr(sample, "locus_order", locus_order):
            if locus_id in locus_indexes:
                applicable[sample_index, locus_indexes[locus_id]] = True
    loci_by_id = {locus.locus_id: locus for locus in (loci or [])}
    recovered, sequence_status = recover_masked_sequences(
        samples, locus_order, loci_by_id
    )
    write_tsv(
        sequence_status,
        output / "combined_marker_sequence_status.tsv",
        SEQUENCE_STATUS_FIELDS,
    )

    snp_sum = np.zeros((sample_count, sample_count), dtype=np.float64)
    repeat_sum = np.zeros((sample_count, sample_count), dtype=np.float64)
    loci_compared = np.zeros((sample_count, sample_count), dtype=np.int32)
    locus_status: list[dict[str, object]] = []
    mafft: str | None = None
    raxml_ng: str | None = None
    built_loci = 0
    snp_path = output / "locus_snp_distances.tsv"
    snp_temporary = snp_path.with_name(snp_path.name + ".tmp")
    with snp_temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOCUS_SNP_FIELDS, delimiter="\t")
        writer.writeheader()
        for locus_index, locus_id in enumerate(locus_order):
            members = [
                index
                for index, sample_id in enumerate(sample_ids)
                if (sample_id, locus_id) in recovered
                and math.isfinite(float(repeat_calls[index, locus_index]))
            ]
            safe_locus = _safe_name(locus_id)
            locus_dir = tree_dir / safe_locus
            locus_dir.mkdir(parents=True, exist_ok=True)
            alignment_path = locus_dir / "samples.repeat_masked.aligned.fasta"
            tree_path = locus_dir / "samples.tree"
            if not members:
                locus_status.append(
                    _locus_status(locus_id, 0, 0, "NO_RECOVERED_SEQUENCES")
                )
                continue
            haplotype_members: dict[str, list[int]] = {}
            for index in members:
                sequence = recovered[(sample_ids[index], locus_id)].sequence
                haplotype_members.setdefault(sequence, []).append(index)
            haplotypes = sorted(
                haplotype_members.items(), key=lambda item: sample_ids[min(item[1])]
            )
            haplotype_records = [
                (f"H{number:06d}", sequence)
                for number, (sequence, _members) in enumerate(haplotypes, start=1)
            ]
            haplotype_matrix = np.zeros((len(haplotypes), len(haplotypes)))
            inference_status = "NO_SNP_VARIATION"
            if len(haplotypes) >= 2:
                if mafft is None:
                    mafft = check_mafft(mafft_bin)
                input_path = locus_dir / "samples.repeat_masked.fasta"
                _write_fasta(haplotype_records, input_path)
                _run_mafft(
                    build_mafft_reference_command(input_path, threads, mafft),
                    alignment_path,
                    f"retrospective sample alignment for {locus_id}",
                )
                aligned = dict(_read_fasta(alignment_path))
                if len(haplotypes) < 4:
                    for left in range(len(haplotypes)):
                        for right in range(left + 1, len(haplotypes)):
                            distance = _aligned_snp_distance(
                                aligned[haplotype_records[left][0]],
                                aligned[haplotype_records[right][0]],
                            )
                            haplotype_matrix[left, right] = distance
                            haplotype_matrix[right, left] = distance
                    inference_status = (
                        "TWO_HAPLOTYPE_DISTANCE"
                        if len(haplotypes) == 2
                        else "THREE_HAPLOTYPE_DISTANCE"
                    )
                else:
                    if raxml_ng is None:
                        raxml_ng = check_raxml_ng(raxml_ng_bin)
                    prefix = locus_dir / "samples"
                    raw_tree = locus_dir / "haplotypes.raxml.tree"
                    _run_raxml_ng(
                        build_raxml_ng_command(
                            alignment_path, prefix, threads, raxml_ng, raxml_model
                        ),
                        prefix,
                        raw_tree,
                        f"retrospective sample tree for {locus_id}",
                    )
                    names, raw_matrix = _tip_patristic_distance_matrix(
                        _parse_newick(raw_tree.read_text(encoding="utf-8"))
                    )
                    indexes = {name: index for index, name in enumerate(names)}
                    order = np.asarray(
                        [indexes[name] for name, _sequence in haplotype_records],
                        dtype=np.intp,
                    )
                    haplotype_matrix = raw_matrix[np.ix_(order, order)]
                    inference_status = "RAXML_NG"
            else:
                _write_fasta(haplotype_records, alignment_path)

            positive = haplotype_matrix[np.triu_indices(len(haplotypes), k=1)]
            positive = positive[positive > 0]
            tree_scale = float(np.median(positive)) if positive.size else 1.0
            haplotype_by_member = {
                member: haplotype_index
                for haplotype_index, (_sequence, indexes) in enumerate(haplotypes)
                for member in indexes
            }
            expanded = np.zeros((len(members), len(members)), dtype=np.float64)
            repeat_values = repeat_calls[members, locus_index]
            repeat_scale = max(
                0.5,
                statistics.pstdev(float(value) for value in repeat_values)
                if len(repeat_values) > 1
                else 0.0,
            )
            for left_position, left in enumerate(members):
                for right_position, right in enumerate(
                    members[left_position + 1 :], start=left_position + 1
                ):
                    distance = float(
                        haplotype_matrix[
                            haplotype_by_member[left], haplotype_by_member[right]
                        ]
                    )
                    expanded[left_position, right_position] = distance
                    expanded[right_position, left_position] = distance
                    normalized = distance / max(tree_scale, 1e-12)
                    repeat_distance = abs(
                        float(repeat_calls[left, locus_index])
                        - float(repeat_calls[right, locus_index])
                    )
                    normalized_repeat = repeat_distance / repeat_scale
                    snp_sum[left, right] += normalized
                    snp_sum[right, left] += normalized
                    writer.writerow(
                        {
                            "locus_id": locus_id,
                            "sample_1": sample_ids[left],
                            "sample_2": sample_ids[right],
                            "snp_patristic_distance": f"{distance:.8f}",
                            "locus_tree_scale": f"{tree_scale:.8f}",
                            "normalized_snp_distance": f"{normalized:.8f}",
                            "repeat_count_distance": f"{repeat_distance:.8f}",
                            "locus_repeat_scale": f"{repeat_scale:.8f}",
                            "normalized_repeat_distance": f"{normalized_repeat:.8f}",
                            "combined_locus_distance": (
                                f"{snp_weight * normalized + repeat_weight * normalized_repeat:.8f}"
                            ),
                        }
                    )

            repeat_distances = np.abs(repeat_values[:, None] - repeat_values[None, :])
            normalized_repeats = repeat_distances / repeat_scale
            for left_position, left in enumerate(members):
                for right_position, right in enumerate(
                    members[left_position + 1 :], start=left_position + 1
                ):
                    repeat_sum[left, right] += normalized_repeats[left_position, right_position]
                    repeat_sum[right, left] += normalized_repeats[left_position, right_position]
                    loci_compared[left, right] += 1
                    loci_compared[right, left] += 1
            tree_path.write_text(
                neighbor_joining_tree_from_matrix(
                    [sample_ids[index] for index in members], expanded
                ),
                encoding="utf-8",
            )
            built_loci += 1
            locus_status.append(
                _locus_status(
                    locus_id,
                    len(members),
                    len(haplotypes),
                    inference_status,
                    tree_scale,
                    alignment_path,
                    tree_path,
                )
            )
    snp_temporary.replace(snp_path)
    write_tsv(locus_status, output / "locus_tree_status.tsv", LOCUS_STATUS_FIELDS)

    combined = np.full((sample_count, sample_count), np.nan, dtype=np.float64)
    np.fill_diagonal(combined, 0.0)
    np.fill_diagonal(loci_compared, len(locus_order))
    pairwise_rows: list[dict[str, object]] = []
    for left in range(sample_count):
        for right in range(left + 1, sample_count):
            compared = int(loci_compared[left, right])
            shared_applicable = int(
                np.count_nonzero(applicable[left] & applicable[right])
            )
            fraction = compared / shared_applicable if shared_applicable else 0.0
            supported = (
                compared >= min_pairwise_loci and fraction >= min_pairwise_fraction
            )
            mean_snp = snp_sum[left, right] / compared if compared else math.nan
            mean_repeat = repeat_sum[left, right] / compared if compared else math.nan
            distance = (
                snp_weight * mean_snp + repeat_weight * mean_repeat
                if supported
                else math.nan
            )
            if supported:
                combined[left, right] = combined[right, left] = distance
            pairwise_rows.append(
                {
                    "sample_1": sample_ids[left],
                    "sample_2": sample_ids[right],
                    "loci_compared": compared,
                    "fraction_loci_compared": f"{fraction:.8f}",
                    "mean_normalized_snp_distance": ""
                    if not compared
                    else f"{mean_snp:.8f}",
                    "mean_normalized_repeat_distance": ""
                    if not compared
                    else f"{mean_repeat:.8f}",
                    "snp_weight": f"{snp_weight:.6f}",
                    "repeat_weight": f"{repeat_weight:.6f}",
                    "combined_marker_distance": ""
                    if not supported
                    else f"{distance:.8f}",
                    "comparison_status": "sufficient"
                    if supported
                    else "insufficient_overlap",
                }
            )
    write_tsv(
        pairwise_rows,
        output / "combined_marker_pairwise_distances.tsv",
        COMBINED_PAIRWISE_FIELDS,
    )
    callable_counts = np.asarray(
        [
            sum((sample_id, locus_id) in recovered for locus_id in locus_order)
            for sample_id in sample_ids
        ],
        dtype=np.int32,
    )
    eligible = np.flatnonzero(callable_counts > 0)
    if eligible.size:
        retained_local = _complete_subset(
            combined[np.ix_(eligible, eligible)],
            callable_counts[eligible],
            [sample_ids[index] for index in eligible],
        )
        retained = [int(eligible[index]) for index in retained_local]
    else:
        retained = []
    retained_ids = [sample_ids[index] for index in retained]
    retained_matrix = combined[np.ix_(retained, retained)]
    matrix_path = output / "combined_marker_distance_matrix.tsv"
    tree_path = output / "combined_marker_nj.tree"
    metadata_output = output / "combined_marker_metadata.tsv"
    _write_matrix(matrix_path, retained_ids, retained_matrix)
    if retained_ids:
        tree_path.write_text(
            neighbor_joining_tree_from_matrix(retained_ids, retained_matrix),
            encoding="utf-8",
        )
    else:
        tree_path.unlink(missing_ok=True)
    _copy_metadata_subset(metadata_path, metadata_output, set(retained_ids))
    return {
        "locus_trees": tree_dir,
        "sequence_status": output / "combined_marker_sequence_status.tsv",
        "locus_status": output / "locus_tree_status.tsv",
        "locus_snp_distances": snp_path,
        "pairwise": output / "combined_marker_pairwise_distances.tsv",
        "distance_matrix": matrix_path,
        "tree": tree_path if retained_ids else "",
        "metadata": metadata_output,
        "loci_built": built_loci,
        "tree_samples": len(retained_ids),
    }


def _locus_status(
    locus_id: str,
    samples: int,
    haplotypes: int,
    status: str,
    scale: float | str = "",
    alignment: Path | str = "",
    tree: Path | str = "",
) -> dict[str, object]:
    return {
        "locus_id": locus_id,
        "samples": samples,
        "snp_haplotypes": haplotypes,
        "tree_scale": "" if scale == "" else f"{float(scale):.8f}",
        "status": status,
        "alignment": str(alignment),
        "tree": str(tree),
    }
