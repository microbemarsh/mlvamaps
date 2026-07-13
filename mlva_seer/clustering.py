from __future__ import annotations

from collections import Counter, defaultdict

from .models import RepeatFeature

try:
    import edlib
except ImportError:  # pragma: no cover - optional accelerated backend
    edlib = None

try:
    import spoars
except ImportError:  # pragma: no cover - dependency error is reported when clustering is used
    spoars = None


def _edit_distance(a: str, b: str) -> int:
    if edlib is not None:
        return int(edlib.align(a, b, mode="NW", task="distance")["editDistance"])
    if len(a) > len(b):
        a, b = b, a
    previous = list(range(len(a) + 1))
    for b_idx, b_base in enumerate(b, start=1):
        current = [b_idx]
        for a_idx, a_base in enumerate(a, start=1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[a_idx] + 1,
                    previous[a_idx - 1] + (a_base != b_base),
                )
            )
        previous = current
    return previous[-1]


def _aligned_identity(a: str, b: str) -> float:
    length = max(len(a), len(b))
    if length == 0:
        return 1.0
    return 1 - (_edit_distance(a, b) / length)


def _poa_consensus(
    sequence_counts: Counter[str],
) -> tuple[str, dict[str, dict[str, int | str]]]:
    if spoars is None:
        raise RuntimeError(
            "VNTR consensus generation requires spoars>=0.1.3; install mlva-seer with its current dependencies"
        )
    sequences = sorted(sequence_counts, key=lambda sequence: (-sequence_counts[sequence], sequence))
    weights = [sequence_counts[sequence] for sequence in sequences]
    graph = spoars.poa(sequences, alignment_type="global", weights=weights)
    consensus = graph.consensus()
    alignment = graph.msa(include_consensus=True)
    aligned_consensus = alignment[-1]
    metrics: dict[str, dict[str, int | str]] = {}
    for sequence, aligned_sequence in zip(sequences, alignment[:-1]):
        insertions = sum(
            sequence_base != "-" and consensus_base == "-"
            for sequence_base, consensus_base in zip(aligned_sequence, aligned_consensus)
        )
        deletions = sum(
            sequence_base == "-" and consensus_base != "-"
            for sequence_base, consensus_base in zip(aligned_sequence, aligned_consensus)
        )
        substitutions = sum(
            sequence_base != "-" and consensus_base != "-" and sequence_base != consensus_base
            for sequence_base, consensus_base in zip(aligned_sequence, aligned_consensus)
        )
        metrics[sequence] = {
            "aligned_repeat_sequence": aligned_sequence,
            "aligned_consensus_sequence": aligned_consensus,
            "insertions_vs_consensus": insertions,
            "deletions_vs_consensus": deletions,
            "substitutions_vs_consensus": substitutions,
            "edit_distance_to_consensus": insertions + deletions + substitutions,
        }
    return consensus, metrics


def _cluster_partition(
    features: list[RepeatFeature], min_identity: float
) -> list[list[RepeatFeature]]:
    by_sequence: dict[str, list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        by_sequence[feature.repeat_sequence].append(feature)
    sequence_bins = sorted(by_sequence.values(), key=lambda items: (-len(items), items[0].repeat_sequence))

    clusters: list[list[RepeatFeature]] = []
    representatives: list[str] = []
    for sequence_bin in sequence_bins:
        sequence = sequence_bin[0].repeat_sequence
        ranked_matches = [
            (_aligned_identity(sequence, representative), cluster_idx)
            for cluster_idx, representative in enumerate(representatives)
        ]
        eligible = [match for match in ranked_matches if match[0] >= min_identity]
        if eligible:
            _identity, cluster_idx = max(eligible, key=lambda item: (item[0], -item[1]))
            clusters[cluster_idx].extend(sequence_bin)
        else:
            representatives.append(sequence)
            clusters.append(list(sequence_bin))
    return clusters


def cluster_vntr_asvs(
    features: list[RepeatFeature], min_identity: float = 0.85
) -> tuple[list[dict], list[tuple[str, str]], list[dict]]:
    """Cluster VNTR reads by locus/count and generate indel-aware POA consensuses."""
    if not 0 <= min_identity <= 1:
        raise ValueError("min_identity must be between 0 and 1")
    grouped: dict[tuple[str, int], list[RepeatFeature]] = defaultdict(list)
    locus_depths: Counter[str] = Counter()
    for feature in features:
        grouped[(feature.locus_id, feature.nearest_integer_repeat_count)].append(feature)
        locus_depths[feature.locus_id] += 1

    clusters_by_locus: dict[str, list[list[RepeatFeature]]] = defaultdict(list)
    for (locus_id, _repeat_count), partition in sorted(grouped.items()):
        clusters_by_locus[locus_id].extend(_cluster_partition(partition, min_identity))

    table: list[dict] = []
    fasta: list[tuple[str, str]] = []
    memberships: list[dict] = []
    for locus_id, locus_clusters in sorted(clusters_by_locus.items()):
        locus_clusters.sort(
            key=lambda cluster: (
                -len(cluster),
                cluster[0].nearest_integer_repeat_count,
                min(feature.repeat_sequence for feature in cluster),
            )
        )
        depth = locus_depths[locus_id]
        for idx, cluster in enumerate(locus_clusters, start=1):
            variant_id = f"{locus_id}_ASV{idx}"
            sequence_counts = Counter(feature.repeat_sequence for feature in cluster)
            consensus, sequence_metrics = _poa_consensus(sequence_counts)
            representative = max(
                cluster,
                key=lambda feature: (
                    sequence_counts[feature.repeat_sequence],
                    feature.repeat_sequence,
                ),
            )
            pattern_counts = Counter(feature.repeat_pattern for feature in cluster)
            consensus_pattern = min(
                pattern_counts,
                key=lambda pattern: (-pattern_counts[pattern], pattern),
            )
            total_insertions = sum(
                int(sequence_metrics[sequence]["insertions_vs_consensus"]) * count
                for sequence, count in sequence_counts.items()
            )
            total_deletions = sum(
                int(sequence_metrics[sequence]["deletions_vs_consensus"]) * count
                for sequence, count in sequence_counts.items()
            )
            total_substitutions = sum(
                int(sequence_metrics[sequence]["substitutions_vs_consensus"]) * count
                for sequence, count in sequence_counts.items()
            )
            edit_distances = [
                int(sequence_metrics[feature.repeat_sequence]["edit_distance_to_consensus"])
                for feature in cluster
            ]
            reads_with_indels = sum(
                1
                for feature in cluster
                if int(sequence_metrics[feature.repeat_sequence]["insertions_vs_consensus"])
                or int(sequence_metrics[feature.repeat_sequence]["deletions_vs_consensus"])
            )
            row = {
                "sample_id": "",
                "locus_id": locus_id,
                "variant_id": variant_id,
                "repeat_count": representative.nearest_integer_repeat_count,
                "support_reads": len(cluster),
                "unique_sequences": len(sequence_counts),
                "frequency": round(len(cluster) / depth, 6) if depth else 0,
                "consensus_pattern": consensus_pattern,
                "consensus_sequence": consensus,
                "consensus_length_bp": len(consensus),
                "reads_with_indels": reads_with_indels,
                "total_insertions": total_insertions,
                "total_deletions": total_deletions,
                "total_substitutions": total_substitutions,
                "mean_edit_distance_to_consensus": round(sum(edit_distances) / len(edit_distances), 4),
                "max_edit_distance_to_consensus": max(edit_distances, default=0),
            }
            table.append(row)
            fasta.append((variant_id, consensus))
            for feature in cluster:
                metrics = sequence_metrics[feature.repeat_sequence]
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
    return table, fasta, memberships
