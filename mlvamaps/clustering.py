from __future__ import annotations

import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import parasail

from .io import gzip_output_file
from .models import Locus, RepeatFeature


_DNA_MATRIX = parasail.matrix_create("ACGTN", 0, -1)
_SIZE_SUFFIX = re.compile(r";size=\d+;?$")


def _check_vsearch(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"VSEARCH executable {executable!r} was not found. Install vsearch "
            "from Bioconda or pass --vsearch-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"Could not run VSEARCH at {path}: {detail}")
    return path


def _reset_work_dir(work_dir: Path) -> None:
    if work_dir.exists():
        for path in work_dir.iterdir():
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
    work_dir.mkdir(parents=True, exist_ok=True)


def build_vsearch_derep_command(
    input_path: Path,
    unique_path: Path,
    uc_path: Path,
    executable: str = "vsearch",
) -> list[str]:
    """Build the native exact-dereplication command for one locus."""
    return [
        executable,
        "--fastx_uniques",
        str(input_path),
        "--fastaout",
        str(unique_path),
        "--uc",
        str(uc_path),
        "--sizeout",
        "--strand",
        "plus",
        "--minseqlength",
        "1",
        "--quiet",
    ]


def build_vsearch_cluster_command(
    unique_path: Path,
    centroid_path: Path,
    uc_path: Path,
    threads: int,
    min_identity: float,
    executable: str = "vsearch",
) -> list[str]:
    """Build abundance-sorted, gap-aware global clustering for one locus."""
    return [
        executable,
        "--cluster_size",
        str(unique_path),
        "--id",
        str(min_identity),
        "--iddef",
        "1",
        "--sizein",
        "--sizeout",
        "--strand",
        "plus",
        "--qmask",
        "none",
        # VNTRs are low-complexity and indels disrupt long exact words. Use a
        # sensitive seed, then let VSEARCH's SIMD global alignment decide.
        "--wordlength",
        "3",
        "--minwordmatches",
        "1",
        # Apply the same affine penalties at sequence ends and internally so
        # an indel is not replaced by a cheap terminal shift.
        "--gapopen",
        "4",
        "--gapext",
        "2",
        "--centroids",
        str(centroid_path),
        "--uc",
        str(uc_path),
        "--threads",
        str(threads),
        "--minseqlength",
        "1",
        "--quiet",
    ]


def _run_vsearch(command: list[str], stage: str) -> None:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"VSEARCH {stage} failed (exit {result.returncode}): {detail}"
        )


def _write_locus_fastqs(
    features: list[RepeatFeature], work_dir: Path
) -> tuple[list[tuple[str, Path]], dict[str, RepeatFeature]]:
    by_locus: dict[str, list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        by_locus[feature.locus_id].append(feature)

    inputs: list[tuple[str, Path]] = []
    feature_by_safe_id: dict[str, RepeatFeature] = {}
    input_dir = work_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    read_idx = 0
    for locus_idx, (locus_id, locus_features) in enumerate(sorted(by_locus.items())):
        path = input_dir / f"locus_{locus_idx:04d}.fastq"
        with path.open("w") as handle:
            for feature in locus_features:
                safe_id = f"mlva_read_{read_idx}"
                read_idx += 1
                sequence = feature.amplicon_sequence
                if not sequence:
                    raise RuntimeError(
                        f"Read {feature.read_id!r} has no primer-trimmed amplicon for VSEARCH"
                    )
                quality = feature.amplicon_quality or ("I" * len(sequence))
                if len(quality) != len(sequence):
                    raise RuntimeError(
                        f"Sequence/quality lengths differ for read {feature.read_id!r}"
                    )
                handle.write(f"@{safe_id}\n{sequence}\n+\n{quality}\n")
                feature_by_safe_id[safe_id] = feature
        inputs.append((locus_id, path))
    return inputs, feature_by_safe_id


def _uc_label(value: str) -> str:
    return _SIZE_SUFFIX.sub("", value.split()[0])


def _read_derep_uc(path: Path) -> dict[str, list[str]]:
    """Map every dereplicated sequence label to its original read labels."""
    members: dict[str, list[str]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or fields[0] not in {"S", "H"} or len(fields) < 10:
                continue
            query = _uc_label(fields[8])
            centroid = query if fields[0] == "S" else _uc_label(fields[9])
            members[centroid].append(query)
    return dict(members)


def _read_cluster_uc(path: Path) -> dict[str, list[str]]:
    """Map each VSEARCH centroid to its dereplicated sequence labels."""
    clusters: dict[str, list[str]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or fields[0] not in {"S", "H"} or len(fields) < 10:
                continue
            query = _uc_label(fields[8])
            centroid = query if fields[0] == "S" else _uc_label(fields[9])
            clusters[centroid].append(query)
    return dict(clusters)


def _cluster_one_locus(
    locus_id: str,
    input_path: Path,
    work_dir: Path,
    threads: int,
    min_identity: float,
    executable: str,
) -> list[tuple[str, list[str]]]:
    locus_dir = work_dir / input_path.stem
    locus_dir.mkdir(parents=True, exist_ok=True)
    uniques = locus_dir / "uniques.fasta"
    derep_uc = locus_dir / "derep.uc"
    centroids = locus_dir / "centroids.fasta"
    clusters_uc = locus_dir / "clusters.uc"

    _run_vsearch(
        build_vsearch_derep_command(input_path, uniques, derep_uc, executable),
        f"exact dereplication for locus {locus_id}",
    )
    _run_vsearch(
        build_vsearch_cluster_command(
            uniques, centroids, clusters_uc, threads, min_identity, executable
        ),
        f"global clustering for locus {locus_id}",
    )

    reads_by_unique = _read_derep_uc(derep_uc)
    unique_by_centroid = _read_cluster_uc(clusters_uc)
    expanded: list[tuple[str, list[str]]] = []
    for centroid, unique_ids in unique_by_centroid.items():
        read_ids = [
            read_id
            for unique_id in unique_ids
            for read_id in reads_by_unique.get(unique_id, [])
        ]
        if not read_ids:
            raise RuntimeError(
                f"VSEARCH cluster at locus {locus_id} could not be expanded to reads"
            )
        expanded.append((centroid, read_ids))
    gzip_output_file(uniques)
    gzip_output_file(centroids)
    return expanded


def _alignment_metrics(query: str, target: str) -> dict[str, int | str]:
    """Return exact edit metrics from a SIMD Parasail global traceback."""
    if not query or not target:
        aligned_query = query or ("-" * len(target))
        aligned_target = target or ("-" * len(query))
        matches = 0
    else:
        result = parasail.nw_trace_scan_32(query, target, 1, 1, _DNA_MATRIX)
        traceback = result.traceback
        aligned_query = traceback.query
        aligned_target = traceback.ref
        matches = traceback.comp.count("|")

    insertions = aligned_target.count("-")
    deletions = aligned_query.count("-")
    substitutions = len(aligned_query) - insertions - deletions - matches
    return {
        "aligned_repeat_sequence": aligned_query,
        "aligned_representative_sequence": aligned_target,
        "insertions_vs_representative": insertions,
        "deletions_vs_representative": deletions,
        "substitutions_vs_representative": substitutions,
        "edit_distance_to_representative": insertions + deletions + substitutions,
    }


def _alignment_metrics_pair(pair: tuple[str, str]) -> dict[str, int | str]:
    return _alignment_metrics(*pair)


def cluster_vntr_asvs(
    features: list[RepeatFeature],
    loci: list[Locus],
    work_dir: str | Path,
    threads: int,
    min_cluster_size: int = 2,
    min_identity: float = 0.97,
    executable: str = "vsearch",
) -> tuple[list[dict], list[tuple[str, str]], list[dict]]:
    """Cluster each locus with VSEARCH and analyze observed centroid sequences."""
    if not features:
        return [], [], []
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be at least 1")
    if not 0.0 <= min_identity <= 1.0:
        raise ValueError("min_identity must be between 0 and 1")
    executable = _check_vsearch(executable)
    work_dir = Path(work_dir)
    _reset_work_dir(work_dir)
    inputs, feature_by_id = _write_locus_fastqs(features, work_dir)

    known_loci = {locus.locus_id for locus in loci}
    locus_totals = Counter(feature.locus_id for feature in features)
    table: list[dict] = []
    fasta: list[tuple[str, str]] = []
    memberships: list[dict] = []
    alignment_executor: ThreadPoolExecutor | None = None

    for locus_id, input_path in inputs:
        if locus_id not in known_loci:
            raise RuntimeError(f"Features contain unknown MLVA locus {locus_id!r}")
        raw_clusters = _cluster_one_locus(
            locus_id, input_path, work_dir, threads, min_identity, executable
        )
        retained = [item for item in raw_clusters if len(item[1]) >= min_cluster_size]
        retained.sort(key=lambda item: (-len(item[1]), item[0]))

        for cluster_number, (centroid_id, member_ids) in enumerate(retained, start=1):
            member_features = [feature_by_id[read_id] for read_id in member_ids]
            representative = feature_by_id.get(centroid_id)
            if representative is None or representative not in member_features:
                raise RuntimeError(
                    f"VSEARCH centroid {centroid_id!r} at locus {locus_id} is not an observed read"
                )
            variant_id = f"{locus_id}_ASV{cluster_number}"
            representative_sequence = representative.repeat_sequence
            unique_repeat_sequences = sorted(
                {feature.repeat_sequence for feature in member_features}
            )
            alignment_pairs = [
                (sequence, representative_sequence) for sequence in unique_repeat_sequences
            ]
            if threads > 1 and len(alignment_pairs) > 1:
                if alignment_executor is None:
                    alignment_executor = ThreadPoolExecutor(max_workers=threads)
                aligned_metrics = alignment_executor.map(
                    _alignment_metrics_pair, alignment_pairs
                )
                metrics_by_sequence = dict(zip(unique_repeat_sequences, aligned_metrics))
            else:
                metrics_by_sequence = {
                    sequence: _alignment_metrics(sequence, representative_sequence)
                    for sequence in unique_repeat_sequences
                }

            metrics_by_read: dict[str, dict[str, int | str]] = {}
            for feature in member_features:
                metrics = metrics_by_sequence[feature.repeat_sequence]
                metrics_by_read[feature.read_id] = metrics
                memberships.append(
                    {
                        "sample_id": "",
                        "read_id": feature.read_id,
                        "locus_id": locus_id,
                        "variant_id": variant_id,
                        "repeat_count": feature.nearest_integer_repeat_count,
                        "repeat_sequence": feature.repeat_sequence,
                        **metrics,
                    }
                )

            cluster_depth = len(member_features)
            edit_distances = [
                int(value["edit_distance_to_representative"])
                for value in metrics_by_read.values()
            ]
            table.append(
                {
                    "sample_id": "",
                    "locus_id": locus_id,
                    "variant_id": variant_id,
                    "repeat_count": representative.nearest_integer_repeat_count,
                    "support_reads": cluster_depth,
                    "unique_sequences": len(unique_repeat_sequences),
                    "frequency": round(cluster_depth / locus_totals[locus_id], 6),
                    "representative_read_id": representative.read_id,
                    "representative_pattern": representative.repeat_pattern,
                    "representative_sequence": representative_sequence,
                    "representative_length_bp": len(representative_sequence),
                    "reads_with_indels": sum(
                        1
                        for value in metrics_by_read.values()
                        if value["insertions_vs_representative"]
                        or value["deletions_vs_representative"]
                    ),
                    "total_insertions": sum(
                        int(value["insertions_vs_representative"])
                        for value in metrics_by_read.values()
                    ),
                    "total_deletions": sum(
                        int(value["deletions_vs_representative"])
                        for value in metrics_by_read.values()
                    ),
                    "total_substitutions": sum(
                        int(value["substitutions_vs_representative"])
                        for value in metrics_by_read.values()
                    ),
                    "mean_edit_distance_to_representative": round(
                        sum(edit_distances) / len(edit_distances), 4
                    )
                    if edit_distances
                    else 0,
                    "max_edit_distance_to_representative": max(
                        edit_distances, default=0
                    ),
                }
            )
            fasta.append((variant_id, representative_sequence))

    if alignment_executor is not None:
        alignment_executor.shutdown()
    for _locus_id, input_path in inputs:
        if input_path.is_file():
            gzip_output_file(input_path)
    return table, fasta, memberships
