from __future__ import annotations

import csv
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pywfa import WavefrontAligner

from .models import Locus, RepeatFeature
from .sequence import find_best, revcomp


MIN_SAVONT_VERSION = (0, 6, 1)


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
) -> tuple[list[Path], dict[str, RepeatFeature], dict[str, str]]:
    by_locus: dict[str, list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        by_locus[feature.locus_id].append(feature)

    inputs: list[Path] = []
    feature_by_safe_id: dict[str, RepeatFeature] = {}
    sample_to_locus: dict[str, str] = {}
    input_dir = work_dir / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    read_idx = 0
    for locus_idx, (locus_id, locus_features) in enumerate(sorted(by_locus.items())):
        stem = f"locus_{locus_idx:04d}"
        path = input_dir / f"{stem}.fastq"
        sample_to_locus[stem] = locus_id
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
    return inputs, feature_by_safe_id, sample_to_locus


def _read_feature_table(path: Path) -> tuple[list[str], dict[int, list[float]]]:
    depths: dict[int, list[float]] = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "#OTU ID":
            raise RuntimeError(f"Unexpected Savont feature table format in {path}")
        samples = header[1:]
        for row in reader:
            match = re.match(r"final_consensus_(\d+)_depth_", row[0])
            if match:
                depths[int(match.group(1))] = [float(value) for value in row[1:]]
    return samples, depths


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


def _read_fasta(path: Path) -> dict[int, str]:
    records: dict[int, str] = {}
    current: int | None = None
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(parts)
                match = re.match(r">final_consensus_(\d+)_depth_", line)
                current = int(match.group(1)) if match else None
                parts = []
            elif current is not None:
                parts.append(line.upper())
    if current is not None:
        records[current] = "".join(parts)
    return records


def _repeat_from_amplicon(sequence: str, locus: Locus) -> str:
    start = 0
    end = len(sequence)
    if locus.forward_primer:
        pos, _ = find_best(locus.forward_primer, sequence, 5)
        start = (pos + len(locus.forward_primer)) if pos is not None else len(locus.forward_primer)
    reverse_binding = revcomp(locus.reverse_primer) if locus.reverse_primer else ""
    if reverse_binding:
        pos, _ = find_best(reverse_binding, sequence[start:], 5)
        end = start + pos if pos is not None else max(start, len(sequence) - len(reverse_binding))
    if locus.left_flank_sequence:
        pos, _ = find_best(locus.left_flank_sequence, sequence[start:end], 3)
        if pos is not None:
            start += pos + len(locus.left_flank_sequence)
    if locus.right_flank_sequence:
        pos, _ = find_best(locus.right_flank_sequence, sequence[start:end], 3)
        if pos is not None:
            end = start + pos
    return sequence[start:end]


def _alignment_metrics(query: str, target: str) -> dict[str, int | str]:
    """Return exact edit metrics from an end-to-end WFA2 traceback."""
    if not query or not target:
        aligned_query = list(query or ("-" * len(target)))
        aligned_target = list(target or ("-" * len(query)))
    else:
        aligner = WavefrontAligner(
            pattern=query,
            distance="affine",
            match=0,
            mismatch=1,
            gap_opening=0,
            gap_extension=1,
            scope="full",
            span="end-to-end",
            heuristic=None,
        )
        result = aligner(target, clip_cigar=False)
        if result.status != 0:
            raise RuntimeError(f"WFA2 global alignment failed with status {result.status}")
        query_index = 0
        target_index = 0
        aligned_query = []
        aligned_target = []
        for operation, length in result.cigartuples:
            if operation in (0, 7, 8):  # M, =, X
                aligned_query.extend(query[query_index : query_index + length])
                aligned_target.extend(target[target_index : target_index + length])
                query_index += length
                target_index += length
            elif operation == 2:  # D: pattern/query base opposite a text gap
                aligned_query.extend(query[query_index : query_index + length])
                aligned_target.extend("-" * length)
                query_index += length
            elif operation == 1:  # I: text/target base opposite a pattern gap
                aligned_query.extend("-" * length)
                aligned_target.extend(target[target_index : target_index + length])
                target_index += length
            else:
                raise RuntimeError(f"Unexpected WFA2 CIGAR operation {operation}")
        if query_index != len(query) or target_index != len(target):
            raise RuntimeError("WFA2 returned an incomplete end-to-end traceback")

    pairs = list(zip(aligned_query, aligned_target))
    insertions = sum(
        1 for query_base, target_base in pairs if query_base != "-" and target_base == "-"
    )
    deletions = sum(
        1 for query_base, target_base in pairs if query_base == "-" and target_base != "-"
    )
    substitutions = sum(
        1
        for query_base, target_base in pairs
        if query_base != "-" and target_base != "-" and query_base != target_base
    )
    return {
        "aligned_repeat_sequence": "".join(aligned_query),
        "aligned_consensus_sequence": "".join(aligned_target),
        "insertions_vs_consensus": insertions,
        "deletions_vs_consensus": deletions,
        "substitutions_vs_consensus": substitutions,
        "edit_distance_to_consensus": insertions + deletions + substitutions,
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
    """Call MLVA ASVs with Savont and annotate per-read indels against its consensuses."""
    if not features:
        return [], [], []
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be at least 1")
    _check_savont(executable)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    inputs, feature_by_id, sample_to_locus = _write_locus_fastqs(features, work_dir)
    lengths = [len(feature.amplicon_sequence) for feature in features]
    command = build_savont_command(
        inputs, work_dir, threads, min(lengths), max(lengths), min_cluster_size, executable
    )
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"Savont ASV calling failed (exit {result.returncode}): {detail}")

    required = [work_dir / name for name in ("final_asvs.fasta", "feature-table.tsv", "final_clusters.tsv")]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise RuntimeError("Savont completed without required output(s): " + ", ".join(missing))

    samples, depths = _read_feature_table(work_dir / "feature-table.tsv")
    clusters = _read_savont_clusters(work_dir / "final_clusters.tsv")
    consensuses = _read_fasta(work_dir / "final_asvs.fasta")
    sample_indexes = {sample_to_locus[name]: idx for idx, name in enumerate(samples) if name in sample_to_locus}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    locus_totals = Counter(feature.locus_id for feature in features)
    per_locus_number: Counter[str] = Counter()
    table: list[dict] = []
    fasta: list[tuple[str, str]] = []
    memberships: list[dict] = []

    alignment_executor: ProcessPoolExecutor | None = None
    for cluster_idx in sorted(consensuses):
        member_features = [feature_by_id[read_id] for read_id in clusters.get(cluster_idx, []) if read_id in feature_by_id]
        member_loci = {feature.locus_id for feature in member_features}
        nonzero = {
            locus_id for locus_id, sample_idx in sample_indexes.items()
            if cluster_idx in depths
            and sample_idx < len(depths[cluster_idx])
            and depths[cluster_idx][sample_idx] > 0
        }
        observed_loci = member_loci | nonzero
        if len(observed_loci) > 1:
            raise RuntimeError(
                f"Savont ASV {cluster_idx} mixed MLVA loci: {', '.join(sorted(observed_loci))}"
            )
        if not member_loci:
            if len(nonzero) != 1:
                continue
            locus_id = next(iter(nonzero))
        else:
            locus_id = next(iter(member_loci))
        locus = locus_by_id[locus_id]
        per_locus_number[locus_id] += 1
        variant_id = f"{locus_id}_ASV{per_locus_number[locus_id]}"
        repeat_consensus = _repeat_from_amplicon(consensuses[cluster_idx], locus)
        count_counts = Counter(feature.nearest_integer_repeat_count for feature in member_features)
        repeat_count = min(count_counts, key=lambda count: (-count_counts[count], count)) if count_counts else round(len(repeat_consensus) / max(len(locus.repeat_motif), 1))
        pattern_counts = Counter(feature.repeat_pattern for feature in member_features)
        consensus_pattern = min(pattern_counts, key=lambda pattern: (-pattern_counts[pattern], pattern)) if pattern_counts else ""
        unique_repeat_sequences = sorted({feature.repeat_sequence for feature in member_features})
        alignment_pairs = [(sequence, repeat_consensus) for sequence in unique_repeat_sequences]
        if threads > 1 and len(alignment_pairs) > 1:
            if alignment_executor is None:
                alignment_executor = ProcessPoolExecutor(max_workers=threads)
            aligned_metrics = alignment_executor.map(_alignment_metrics_pair, alignment_pairs)
            metrics_by_sequence = dict(zip(unique_repeat_sequences, aligned_metrics))
        else:
            metrics_by_sequence = {
                sequence: _alignment_metrics(sequence, repeat_consensus)
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
        locus_depth = float(depths.get(cluster_idx, [0.0] * len(samples))[sample_indexes[locus_id]]) if locus_id in sample_indexes else float(len(member_features))
        edit_distances = [int(value["edit_distance_to_consensus"]) for value in metrics_by_read.values()]
        table.append({
            "sample_id": "",
            "locus_id": locus_id,
            "variant_id": variant_id,
            "repeat_count": repeat_count,
            "support_reads": int(round(locus_depth)),
            "unique_sequences": len({feature.repeat_sequence for feature in member_features}),
            "frequency": round(locus_depth / locus_totals[locus_id], 6) if locus_totals[locus_id] else 0,
            "consensus_pattern": consensus_pattern,
            "consensus_sequence": repeat_consensus,
            "consensus_length_bp": len(repeat_consensus),
            "reads_with_indels": sum(1 for value in metrics_by_read.values() if value["insertions_vs_consensus"] or value["deletions_vs_consensus"]),
            "total_insertions": sum(int(value["insertions_vs_consensus"]) for value in metrics_by_read.values()),
            "total_deletions": sum(int(value["deletions_vs_consensus"]) for value in metrics_by_read.values()),
            "total_substitutions": sum(int(value["substitutions_vs_consensus"]) for value in metrics_by_read.values()),
            "mean_edit_distance_to_consensus": round(sum(edit_distances) / len(edit_distances), 4) if edit_distances else 0,
            "max_edit_distance_to_consensus": max(edit_distances, default=0),
        })
        fasta.append((variant_id, repeat_consensus))
    if alignment_executor is not None:
        alignment_executor.shutdown()
    return table, fasta, memberships
