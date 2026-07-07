from __future__ import annotations

from collections import defaultdict

from .models import RepeatFeature
from .sequence import identity, majority_consensus


def cluster_vntr_asvs(features: list[RepeatFeature], min_identity: float = 0.85) -> tuple[list[dict], list[tuple[str, str]]]:
    grouped: dict[str, list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        grouped[feature.locus_id].append(feature)

    table = []
    fasta = []
    for locus_id, locus_features in grouped.items():
        clusters: list[list[RepeatFeature]] = []
        for feature in locus_features:
            placed = False
            for cluster in clusters:
                representative = cluster[0]
                same_count = feature.nearest_integer_repeat_count == representative.nearest_integer_repeat_count
                similar = identity(feature.repeat_sequence, representative.repeat_sequence) >= min_identity
                if same_count and similar:
                    cluster.append(feature)
                    placed = True
                    break
            if not placed:
                clusters.append([feature])
        clusters.sort(key=len, reverse=True)
        depth = len(locus_features)
        for idx, cluster in enumerate(clusters, start=1):
            variant_id = f"{locus_id}_ASV{idx}"
            consensus = majority_consensus([feature.repeat_sequence for feature in cluster])
            representative = cluster[0]
            row = {
                "sample_id": "",
                "locus_id": locus_id,
                "variant_id": variant_id,
                "repeat_count": representative.nearest_integer_repeat_count,
                "support_reads": len(cluster),
                "frequency": round(len(cluster) / depth, 6) if depth else 0,
                "consensus_pattern": representative.repeat_pattern,
                "consensus_sequence": consensus,
            }
            table.append(row)
            fasta.append((variant_id, consensus))
    return table, fasta
