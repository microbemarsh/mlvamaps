from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

from .calling import legacy_round_repeat_count, normalize_allele
from .clustering import _alignment_metrics, _alignment_metrics_pair
from .concurrency import DEFAULT_THREADS, resolve_threads
from .models import RepeatFeature


def _mapped_group_label(feature: RepeatFeature, recruitment: dict | None) -> str:
    """Use the competitive product mapping as the repeat-variant partition."""
    if recruitment is not None:
        candidate = recruitment.get("candidate_allele")
        if candidate not in ("", None):
            try:
                return str(normalize_allele(float(candidate)))
            except (TypeError, ValueError):
                pass
    # Primer-only mode has no competitive mapping label. Keep its reads
    # together so POA, rather than raw read measurements, defines the allele.
    return "PRIMARY"


def mapped_read_variant_groups(
    features: list[RepeatFeature],
    recruitment_rows: list[dict],
    sample_id: str,
    threads: int = DEFAULT_THREADS,
) -> tuple[list[dict], list[tuple[str, str]], list[dict]]:
    """Build locus/allele evidence groups from competitive read mappings.

    The returned schemas intentionally match the historical ASV tables so the
    confidence and mixture callers can consume mapping evidence without a
    sequence-clustering step.
    """
    recruitment_by_read = {
        str(row["read_id"]): row
        for row in recruitment_rows
        if row.get("genotype_informative") == "yes"
    }
    grouped: dict[tuple[str, str], list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        grouped[
            (
                feature.locus_id,
                _mapped_group_label(
                    feature, recruitment_by_read.get(feature.read_id)
                ),
            )
        ].append(feature)

    locus_totals = Counter(feature.locus_id for feature in features)
    table: list[dict] = []
    fasta: list[tuple[str, str]] = []
    memberships: list[dict] = []
    groups_by_locus: dict[str, list[tuple[str, list[RepeatFeature]]]] = (
        defaultdict(list)
    )
    for (locus_id, label), group_features in grouped.items():
        groups_by_locus[locus_id].append((label, group_features))

    thread_count = resolve_threads(threads)
    alignment_executor = (
        ThreadPoolExecutor(max_workers=thread_count)
        if thread_count > 1
        else None
    )
    metrics_cache: dict[tuple[str, str], dict[str, int | str]] = {}
    try:
        for locus_id, locus_groups in sorted(groups_by_locus.items()):
            locus_groups.sort(
                key=lambda item: (
                    -len(item[1]),
                    item[0],
                )
            )
            for group_number, (label, group_features) in enumerate(
                locus_groups, start=1
            ):
                representative = max(
                    group_features,
                    key=lambda feature: (
                        feature.flank_quality_score,
                        feature.mean_qscore,
                        -feature.mismatch_count_in_repeat_region,
                        feature.read_id,
                    ),
                )
                variant_id = f"{locus_id}_MAP{group_number}"
                representative_sequence = representative.repeat_sequence
                unique_sequences = sorted(
                    {feature.repeat_sequence for feature in group_features}
                )
                missing_pairs = [
                    (sequence, representative_sequence)
                    for sequence in unique_sequences
                    if (sequence, representative_sequence) not in metrics_cache
                ]
                if alignment_executor is not None and len(missing_pairs) > 1:
                    missing_metrics = alignment_executor.map(
                        _alignment_metrics_pair, missing_pairs
                    )
                    metrics_cache.update(zip(missing_pairs, missing_metrics))
                else:
                    for pair in missing_pairs:
                        metrics_cache[pair] = _alignment_metrics(*pair)
                metrics_by_read = {
                    feature.read_id: metrics_cache[
                        (feature.repeat_sequence, representative_sequence)
                    ]
                    for feature in group_features
                }
                for feature in group_features:
                    metrics = metrics_by_read[feature.read_id]
                    memberships.append(
                        {
                            "sample_id": sample_id,
                            "read_id": feature.read_id,
                            "locus_id": locus_id,
                            "variant_id": variant_id,
                            "repeat_count": (
                                label
                                if label != "PRIMARY"
                                else legacy_round_repeat_count(
                                    feature.raw_repeat_count_estimate
                                )
                            ),
                            "repeat_sequence": feature.repeat_sequence,
                            **metrics,
                        }
                    )
                edit_distances = [
                    int(metrics["edit_distance_to_representative"])
                    for metrics in metrics_by_read.values()
                ]
                repeat_count = (
                    normalize_allele(float(label))
                    if label != "PRIMARY"
                    else legacy_round_repeat_count(
                        representative.raw_repeat_count_estimate
                    )
                )
                row = {
                    "sample_id": sample_id,
                    "locus_id": locus_id,
                    "variant_id": variant_id,
                    "repeat_count": repeat_count,
                    "support_reads": len(group_features),
                    "unique_sequences": len(unique_sequences),
                    "frequency": round(
                        len(group_features) / max(locus_totals[locus_id], 1), 6
                    ),
                    "representative_read_id": representative.read_id,
                    "representative_pattern": representative.repeat_pattern,
                    "representative_sequence": representative_sequence,
                    "representative_length_bp": len(representative_sequence),
                    "reads_with_indels": sum(
                        bool(
                            metrics["insertions_vs_representative"]
                            or metrics["deletions_vs_representative"]
                        )
                        for metrics in metrics_by_read.values()
                    ),
                    "total_insertions": sum(
                        int(metrics["insertions_vs_representative"])
                        for metrics in metrics_by_read.values()
                    ),
                    "total_deletions": sum(
                        int(metrics["deletions_vs_representative"])
                        for metrics in metrics_by_read.values()
                    ),
                    "total_substitutions": sum(
                        int(metrics["substitutions_vs_representative"])
                        for metrics in metrics_by_read.values()
                    ),
                    "mean_edit_distance_to_representative": round(
                        sum(edit_distances) / max(len(edit_distances), 1), 4
                    ),
                    "max_edit_distance_to_representative": max(
                        edit_distances, default=0
                    ),
                }
                table.append(row)
                fasta.append((variant_id, representative_sequence))
    finally:
        if alignment_executor is not None:
            alignment_executor.shutdown()
    return table, fasta, memberships
