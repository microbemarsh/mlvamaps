from __future__ import annotations

import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import parasail

from .models import Locus, RepeatFeature


MIN_SAVONT_VERSION = (0, 6, 1)
_DNA_MATRIX = parasail.matrix_create("ACGTN", 0, -1)


def _version_tuple(value: str) -> tuple[int, ...]:
    match = re.search(r"(\d+)\.(\d+)\.(\d+)", value)
    return tuple(map(int, match.groups())) if match else ()


def _check_savont(executable: str) -> None:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"Savont executable {executable!r} was not found. Install savont>=0.6.1 "
            "from Bioconda or pass --savont-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    version = _version_tuple(result.stdout + result.stderr)
    if not version or version < MIN_SAVONT_VERSION:
        found = ".".join(map(str, version)) if version else "unknown"
        raise RuntimeError(f"mlva-seer requires savont>=0.6.1; found {found} at {path}")


def _is_pooled_quantification_index_panic(detail: str) -> bool:
    return (
        "panicked at" in detail
        and "index out of bounds" in detail
        and "Per-sample quantification" in detail
    )


def _clear_partial_savont_outputs(work_dir: Path) -> None:
    """Remove failed Savont outputs while retaining MLVA Seer's input FASTQs."""
    for path in work_dir.iterdir():
        if path.name == "inputs":
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def build_savont_command(
    input_paths: list[Path],
    output_dir: Path,
    threads: int,
    min_read_length: int,
    max_read_length: int,
    min_cluster_size: int,
    executable: str = "savont",
) -> list[str]:
    """Build the single pooled Savont ASV command used by the read pipeline."""
    return [
        executable,
        "asv",
        *map(str, input_paths),
        "--output-dir",
        str(output_dir),
        "--threads",
        str(threads),
        "--pooled-samples",
        "--single-strand",
        "--min-read-length",
        str(min_read_length),
        "--max-read-length",
        str(max_read_length),
        "--min-cluster-size",
        str(min_cluster_size),
    ]


def _write_locus_fastqs(
    features: list[RepeatFeature], work_dir: Path
) -> tuple[list[Path], dict[str, RepeatFeature]]:
    by_locus: dict[str, list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        by_locus[feature.locus_id].append(feature)

    inputs: list[Path] = []
    feature_by_safe_id: dict[str, RepeatFeature] = {}
    input_dir = work_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    read_idx = 0
    for locus_idx, (locus_id, locus_features) in enumerate(sorted(by_locus.items())):
        stem = f"locus_{locus_idx:04d}"
        path = input_dir / f"{stem}.fastq"
        with path.open("w") as handle:
            for feature in locus_features:
                safe_id = f"mlva_read_{read_idx}"
                read_idx += 1
                sequence = feature.amplicon_sequence
                if not sequence:
                    raise RuntimeError(
                        f"Read {feature.read_id!r} has no primer-trimmed amplicon for Savont"
                    )
                quality = feature.amplicon_quality or ("I" * len(sequence))
                if len(quality) != len(sequence):
                    raise RuntimeError(f"Sequence/quality lengths differ for read {feature.read_id!r}")
                handle.write(f"@{safe_id}\n{sequence}\n+\n{quality}\n")
                feature_by_safe_id[safe_id] = feature
        inputs.append(path)
    return inputs, feature_by_safe_id


def _read_savont_clusters(path: Path) -> dict[int, list[str]]:
    clusters: dict[int, list[str]] = {}
    current: int | None = None
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            match = re.match(r"final_cluster_(\d+)$", fields[0])
            if match:
                current = int(match.group(1))
                clusters[current] = []
            elif current is not None and fields[0]:
                clusters[current].append(fields[0].split()[0])
    return clusters


def _select_representative(cluster: list[RepeatFeature]) -> RepeatFeature:
    """Select a real, well-supported read sequence to represent a cluster."""
    sequence_counts = Counter(feature.repeat_sequence for feature in cluster)
    sequence_quality: dict[str, float] = defaultdict(float)
    for feature in cluster:
        sequence_quality[feature.repeat_sequence] += feature.mean_qscore
    representative_sequence = min(
        sequence_counts,
        key=lambda sequence: (
            -sequence_counts[sequence],
            -(sequence_quality[sequence] / sequence_counts[sequence]),
            sequence,
        ),
    )
    candidates = [
        feature for feature in cluster if feature.repeat_sequence == representative_sequence
    ]
    return min(
        candidates,
        key=lambda feature: (
            -feature.mean_qscore,
            -feature.flank_quality_score,
            feature.read_id,
        ),
    )


def _alignment_metrics(query: str, target: str) -> dict[str, int | str]:
    """Return exact edit metrics from a SIMD Parasail global traceback."""
    if not query or not target:
        aligned_query = query or ("-" * len(target))
        aligned_target = target or ("-" * len(query))
        matches = 0
    else:
        # With match=0, mismatch=-1, and open/extend=1, the optimal score is
        # the negative Levenshtein distance. The scan_32 implementation uses
        # SIMD lanes while retaining an exact end-to-end traceback.
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
    executable: str = "savont",
) -> tuple[list[dict], list[tuple[str, str]], list[dict]]:
    """Use Savont cluster membership, then analyze clusters using observed representatives."""
    if not features:
        return [], [], []
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be at least 1")
    _check_savont(executable)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    _clear_partial_savont_outputs(work_dir)
    shutil.rmtree(work_dir / "inputs", ignore_errors=True)
    inputs, feature_by_id = _write_locus_fastqs(features, work_dir)
    lengths = [len(feature.amplicon_sequence) for feature in features]
    command = build_savont_command(
        inputs, work_dir, threads, min(lengths), max(lengths), min_cluster_size, executable
    )
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    pooled_failure_detail = ""
    if result.returncode:
        pooled_failure_detail = (result.stderr or result.stdout).strip()
        if _is_pooled_quantification_index_panic(pooled_failure_detail):
            _clear_partial_savont_outputs(work_dir)
            fallback_command = [argument for argument in command if argument != "--pooled-samples"]
            result = subprocess.run(
                fallback_command, capture_output=True, text=True, check=False
            )
            if not result.returncode:
                (work_dir / "pooled_quantification_fallback.log").write_text(
                    pooled_failure_detail + "\n"
                )
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        if pooled_failure_detail:
            detail = (
                "pooled run failed with a Savont per-sample index panic; "
                f"non-pooled retry also failed:\n{detail}"
            )
        raise RuntimeError(f"Savont ASV calling failed (exit {result.returncode}): {detail}")

    cluster_path = work_dir / "final_clusters.tsv"
    if not cluster_path.exists():
        raise RuntimeError(f"Savont completed without required output: {cluster_path}")

    clusters = _read_savont_clusters(work_dir / "final_clusters.tsv")
    known_loci = {locus.locus_id for locus in loci}
    locus_totals = Counter(feature.locus_id for feature in features)
    per_locus_number: Counter[str] = Counter()
    table: list[dict] = []
    fasta: list[tuple[str, str]] = []
    memberships: list[dict] = []

    alignment_executor: ThreadPoolExecutor | None = None
    for cluster_idx in sorted(clusters):
        member_features = [feature_by_id[read_id] for read_id in clusters.get(cluster_idx, []) if read_id in feature_by_id]
        member_loci = {feature.locus_id for feature in member_features}
        if len(member_loci) > 1:
            raise RuntimeError(
                f"Savont ASV {cluster_idx} mixed MLVA loci: {', '.join(sorted(member_loci))}"
            )
        if not member_loci:
            continue
        locus_id = next(iter(member_loci))
        if locus_id not in known_loci:
            raise RuntimeError(f"Savont ASV {cluster_idx} has unknown MLVA locus {locus_id!r}")
        per_locus_number[locus_id] += 1
        variant_id = f"{locus_id}_ASV{per_locus_number[locus_id]}"
        representative = _select_representative(member_features)
        representative_sequence = representative.repeat_sequence
        unique_repeat_sequences = sorted({feature.repeat_sequence for feature in member_features})
        alignment_pairs = [(sequence, representative_sequence) for sequence in unique_repeat_sequences]
        if threads > 1 and len(alignment_pairs) > 1:
            if alignment_executor is None:
                alignment_executor = ThreadPoolExecutor(max_workers=threads)
            aligned_metrics = alignment_executor.map(_alignment_metrics_pair, alignment_pairs)
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
            memberships.append({
                "sample_id": "",
                "read_id": feature.read_id,
                "locus_id": locus_id,
                "variant_id": variant_id,
                "repeat_count": feature.nearest_integer_repeat_count,
                "repeat_sequence": feature.repeat_sequence,
                **metrics,
            })
        locus_depth = len(member_features)
        edit_distances = [int(value["edit_distance_to_representative"]) for value in metrics_by_read.values()]
        table.append({
            "sample_id": "",
            "locus_id": locus_id,
            "variant_id": variant_id,
            "repeat_count": representative.nearest_integer_repeat_count,
            "support_reads": locus_depth,
            "unique_sequences": len({feature.repeat_sequence for feature in member_features}),
            "frequency": round(locus_depth / locus_totals[locus_id], 6) if locus_totals[locus_id] else 0,
            "representative_read_id": representative.read_id,
            "representative_pattern": representative.repeat_pattern,
            "representative_sequence": representative_sequence,
            "representative_length_bp": len(representative_sequence),
            "reads_with_indels": sum(1 for value in metrics_by_read.values() if value["insertions_vs_representative"] or value["deletions_vs_representative"]),
            "total_insertions": sum(int(value["insertions_vs_representative"]) for value in metrics_by_read.values()),
            "total_deletions": sum(int(value["deletions_vs_representative"]) for value in metrics_by_read.values()),
            "total_substitutions": sum(int(value["substitutions_vs_representative"]) for value in metrics_by_read.values()),
            "mean_edit_distance_to_representative": round(sum(edit_distances) / len(edit_distances), 4) if edit_distances else 0,
            "max_edit_distance_to_representative": max(edit_distances, default=0),
        })
        fasta.append((variant_id, representative_sequence))
    if alignment_executor is not None:
        alignment_executor.shutdown()
    return table, fasta, memberships
