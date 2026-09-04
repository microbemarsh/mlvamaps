from __future__ import annotations

import csv
import json
import logging
import math
import random
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


CALIBRATION_SCHEMA_VERSION = "1.0"
CHANNELS = ("repeat", "snp", "joint")
LOGGER = logging.getLogger(__name__)

TAXON_LOCUS_DISCRIMINATION_FIELDS = [
    "locus_id",
    "taxonomic_weight",
    "normalized_information_gain",
    "reference_coverage",
    "references_compared",
    "taxa_compared",
    "marker_signatures",
]

TAXON_ASSIGNMENT_FIELDS = [
    "sample_id",
    "target_taxon_id",
    "target_taxon_name",
    "decision",
    "decision_reason",
    "prediction_set",
    "target_joint_p_value",
    "target_snp_p_value",
    "target_repeat_p_value",
    "target_bootstrap_support",
    "bootstrap_tie_fraction",
    "target_distance",
    "best_alternative_taxon_id",
    "best_alternative_taxon_name",
    "best_alternative_joint_p_value",
    "best_alternative_distance",
    "distance_margin",
    "bootstrap_margin_lower",
    "bootstrap_margin_upper",
    "callable_loci",
    "placed_loci",
    "target_favoring_loci",
    "alternative_favoring_loci",
    "unresolved_loci",
    "mean_placement_entropy",
    "median_best_lwr",
    "qc_status",
    "qc_flags",
    "calibration_version",
]

TAXON_CANDIDATE_FIELDS = [
    "sample_id",
    "taxon_id",
    "taxon_name",
    "reference_count",
    "repeat_distance",
    "snp_distance",
    "joint_distance",
    "repeat_nonconformity",
    "snp_nonconformity",
    "joint_nonconformity",
    "repeat_p_value",
    "snp_p_value",
    "joint_p_value",
    "repeat_accepted",
    "snp_accepted",
    "joint_accepted",
    "nearest_references",
    "rank",
]

CALIBRATION_SCORE_FIELDS = [
    "reference_id",
    "taxon_id",
    "channel",
    "within_distance",
    "outside_distance",
    "nonconformity",
]

TAXON_LOCUS_FIELDS = [
    "sample_id",
    "locus_id",
    "target_taxon_id",
    "alternative_taxon_id",
    "target_repeat_distance",
    "alternative_repeat_distance",
    "target_snp_distance",
    "alternative_snp_distance",
    "target_joint_distance",
    "alternative_joint_distance",
    "joint_distance_margin",
    "placement_entropy",
    "best_placement_lwr",
    "interpretation",
]

TAXONOMIC_SUMMARY_FIELDS = [
    "sample_id",
    "best_taxon",
    "best_species",
    "taxon_score",
    "second_best_taxon",
    "second_best_score",
    "score_margin",
    "informative_loci",
    "expected_loci",
    "locus_recovery_fraction",
    "taxonomic_status",
    "status_reason",
    "assignment",
    "assignment_taxon_id",
    "assignment_rank",
    "assignment_status",
    "confidence",
    "closest_taxon",
    "closest_taxon_id",
    "closest_distance",
    "runner_up_taxon",
    "runner_up_taxon_id",
    "runner_up_distance",
    "distance_margin",
    "relative_margin",
    "loci_recovered",
    "informative_loci_recovered",
    "discriminative_loci_recovered",
    "loci_supporting_assignment",
    "conflicting_loci",
    "bootstrap_support",
    "bootstrap_replicates",
    "input_mode",
]

TAXONOMIC_EVIDENCE_FIELDS = [
    "sample_id",
    "taxon_id",
    "species",
    "rank",
    "score",
    "distance",
    "references_compared",
    "informative_loci",
    "locus_recovery_fraction",
    "is_best_taxon",
    "compatibility",
    "bootstrap_win_fraction",
]

TAXONOMIC_LOCUS_EVIDENCE_FIELDS = [
    "sample_id",
    "locus_id",
    "recovered",
    "taxonomic_weight",
    "discriminative",
    "favored_taxon_id",
    "favored_taxon",
    "closest_taxon_distance",
    "runner_up_distance",
    "distance_margin",
    "supports_assignment",
    "conflicts_with_assignment",
    "depth",
    "consensus_strength",
    "quality_status",
]


@dataclass(frozen=True)
class TaxonCalibration:
    """Persisted label-conditional nonconformity distributions.

    A conformal p-value is a compatibility measure under the reference cohort,
    not a posterior probability that a taxon is present.
    """

    panel_sha256: str
    database_signature: str
    taxon_counts: dict[str, int]
    scores: dict[str, dict[str, tuple[float, ...]]]
    k: int = 3
    alpha: float = 0.05
    snp_weight: float = 1.0
    repeat_weight: float = 1.0
    minimum_loci: int = 3
    schema_version: str = CALIBRATION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != CALIBRATION_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported taxon calibration schema {self.schema_version!r}"
            )
        if self.k < 1:
            raise ValueError("Taxon calibration k must be at least 1")
        if not 0 < self.alpha < 1:
            raise ValueError("Taxon calibration alpha must be between 0 and 1")
        if self.minimum_loci < 1:
            raise ValueError("Taxon calibration minimum_loci must be at least 1")
        if self.snp_weight < 0 or self.repeat_weight < 0:
            raise ValueError("Taxon calibration marker weights cannot be negative")
        if self.snp_weight + self.repeat_weight <= 0:
            raise ValueError("At least one taxon calibration marker weight must be positive")
        for channel in CHANNELS:
            if channel not in self.scores:
                raise ValueError(f"Taxon calibration lacks {channel!r} scores")
            for taxon_id, values in self.scores[channel].items():
                if taxon_id not in self.taxon_counts:
                    raise ValueError(
                        f"Calibration scores refer to unknown taxon {taxon_id!r}"
                    )
                if not values or any(not math.isfinite(value) or value < 0 for value in values):
                    raise ValueError(
                        f"Calibration scores for {taxon_id!r}/{channel} must be "
                        "non-empty, finite, and non-negative"
                    )

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "panel_sha256": self.panel_sha256,
            "database_signature": self.database_signature,
            "taxon_counts": dict(sorted(self.taxon_counts.items())),
            "scores": {
                channel: {
                    taxon_id: list(values)
                    for taxon_id, values in sorted(self.scores[channel].items())
                }
                for channel in CHANNELS
            },
            "k": self.k,
            "alpha": self.alpha,
            "snp_weight": self.snp_weight,
            "repeat_weight": self.repeat_weight,
            "minimum_loci": self.minimum_loci,
        }

    @classmethod
    def from_dict(cls, document: Mapping) -> "TaxonCalibration":
        scores_document = document.get("scores", {})
        return cls(
            schema_version=str(document.get("schema_version", "")),
            panel_sha256=str(document.get("panel_sha256", "")),
            database_signature=str(document.get("database_signature", "")),
            taxon_counts={
                str(key): int(value)
                for key, value in dict(document.get("taxon_counts", {})).items()
            },
            scores={
                channel: {
                    str(taxon_id): tuple(float(value) for value in values)
                    for taxon_id, values in dict(
                        scores_document.get(channel, {})
                    ).items()
                }
                for channel in CHANNELS
            },
            k=int(document.get("k", 3)),
            alpha=float(document.get("alpha", 0.05)),
            snp_weight=float(document.get("snp_weight", 1.0)),
            repeat_weight=float(document.get("repeat_weight", 1.0)),
            minimum_loci=int(document.get("minimum_loci", 3)),
        )

    def write(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n")
        return output

    @classmethod
    def read(cls, path: str | Path) -> "TaxonCalibration":
        calibration_path = Path(path)
        try:
            document = json.loads(calibration_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"Could not read taxon calibration {calibration_path}: {exc}"
            ) from exc
        return cls.from_dict(document)


@dataclass(frozen=True)
class TaxonAssignment:
    summary: dict
    candidates: tuple[dict, ...]
    loci: tuple[dict, ...]


@dataclass(frozen=True)
class AutomaticTaxonAssignment:
    summary: dict
    evidence: tuple[dict, ...]
    loci: tuple[dict, ...] = ()


def _finite(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def conformal_p_value(query_score: float, calibration_scores: Sequence[float]) -> float:
    """Return a finite-sample, label-conditional conformal p-value."""
    if not math.isfinite(query_score) or query_score < 0:
        raise ValueError("Query nonconformity must be finite and non-negative")
    values = [float(value) for value in calibration_scores]
    if not values:
        raise ValueError("At least one calibration score is required")
    if any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("Calibration scores must be finite and non-negative")
    return (1 + sum(value >= query_score for value in values)) / (len(values) + 1)


def _mean_k_smallest(values: Iterable[float], k: int) -> float | None:
    ordered = sorted(value for value in values if math.isfinite(value))
    if not ordered:
        return None
    selected = ordered[: min(k, len(ordered))]
    return sum(selected) / len(selected)


def _percentile(values: Sequence[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def _reference_taxa(
    metadata: Mapping[str, Mapping[str, str]],
) -> tuple[dict[str, str], dict[str, str]]:
    reference_taxon: dict[str, str] = {}
    taxon_names: dict[str, str] = {}
    for reference_id, row in metadata.items():
        taxon_id = str(
            row.get("taxon_id") or row.get("taxid") or row.get("ncbi_taxid") or ""
        ).strip()
        if not taxon_id:
            continue
        reference_taxon[str(reference_id)] = taxon_id
        name = str(
            row.get("taxon_name")
            or row.get("species")
            or row.get("organism_name")
            or row.get("scientific_name")
            or ""
        ).strip()
        if name:
            previous = taxon_names.get(taxon_id)
            if previous is not None and previous != name:
                raise ValueError(
                    f"Taxon {taxon_id!r} has conflicting names {previous!r} and {name!r}"
                )
            taxon_names[taxon_id] = name
    return reference_taxon, taxon_names


def _marker_signature(row: Mapping) -> tuple[str, str] | None:
    """Return the stable repeat/SNP signature stored in the reference index."""
    repeat = str(row.get("repeat_count", "")).strip()
    snp = str(row.get("snp_sha256", "")).strip()
    if not repeat and not snp:
        return None
    return repeat, snp


def derive_locus_discrimination_weights(
    reference_index_rows: Iterable[Mapping],
    reference_metadata: Mapping[str, Mapping[str, str]],
) -> list[dict]:
    """Derive deterministic, explainable taxonomic information weights.

    The weight is normalized mutual information between a locus marker signature
    (repeat count plus repeat-masked SNP digest) and taxon label, multiplied by
    reference coverage.  It is zero for a marker shared identically among taxa
    and approaches one when signatures uniquely identify represented taxa.
    """
    reference_taxon, _names = _reference_taxa(reference_metadata)
    labeled_references = set(reference_taxon)
    by_locus: dict[str, list[tuple[str, tuple[str, str]]]] = {}
    for row in reference_index_rows:
        reference_id = str(row.get("reference_id", "")).strip()
        locus_id = str(row.get("locus_id", "")).strip()
        signature = _marker_signature(row)
        if reference_id in reference_taxon and locus_id and signature is not None:
            by_locus.setdefault(locus_id, []).append(
                (reference_taxon[reference_id], signature)
            )
    results = []
    for locus_id in sorted(by_locus):
        observations = by_locus[locus_id]
        total = len(observations)
        taxon_counts: dict[str, int] = {}
        signature_counts: dict[tuple[str, str], int] = {}
        joint_counts: dict[tuple[str, tuple[str, str]], int] = {}
        for taxon_id, signature in observations:
            taxon_counts[taxon_id] = taxon_counts.get(taxon_id, 0) + 1
            signature_counts[signature] = signature_counts.get(signature, 0) + 1
            key = (taxon_id, signature)
            joint_counts[key] = joint_counts.get(key, 0) + 1
        taxon_entropy = -sum(
            (count / total) * math.log(count / total, 2)
            for count in taxon_counts.values()
        )
        mutual_information = sum(
            (count / total)
            * math.log(
                (count * total)
                / (taxon_counts[taxon_id] * signature_counts[signature]),
                2,
            )
            for (taxon_id, signature), count in joint_counts.items()
        )
        information_gain = (
            mutual_information / taxon_entropy if taxon_entropy > 0 else 0.0
        )
        coverage = total / len(labeled_references) if labeled_references else 0.0
        # Singleton signatures cannot demonstrate a reproducible taxon pattern;
        # discount them rather than overfitting isolate-specific SNP hashes.
        replicated_fraction = sum(
            count for count in signature_counts.values() if count >= 2
        ) / total
        results.append(
            {
                "locus_id": locus_id,
                "taxonomic_weight": max(
                    0.0, min(1.0, information_gain * coverage * replicated_fraction)
                ),
                "normalized_information_gain": information_gain,
                "reference_coverage": coverage,
                "references_compared": total,
                "taxa_compared": len(taxon_counts),
                "marker_signatures": len(signature_counts),
            }
        )
    return results


def read_locus_discrimination_weights(
    path: str | Path | None,
) -> dict[str, float]:
    if path is None or not Path(path).is_file():
        return {}
    with Path(path).open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {
            str(row.get("locus_id", "")): float(row.get("taxonomic_weight", 0) or 0)
            for row in rows
            if str(row.get("locus_id", ""))
        }


def _locus_reference_values(
    locus_marker_rows: Iterable[Mapping],
    reference_taxon: Mapping[str, str],
    snp_weight: float,
    repeat_weight: float,
) -> dict[str, dict[str, dict[str, dict[str, float]]]]:
    """Return locus -> taxon -> channel -> reference -> distance."""
    result: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
    for row in locus_marker_rows:
        locus_id = str(row.get("locus_id", ""))
        reference_id = str(row.get("reference_id", ""))
        taxon_id = reference_taxon.get(reference_id)
        if not locus_id or taxon_id is None:
            continue
        repeat = _finite(row.get("normalized_repeat_distance"))
        snp = _finite(row.get("normalized_snp_distance"))
        channels: dict[str, float] = {}
        if repeat is not None:
            channels["repeat"] = repeat
        if snp is not None:
            channels["snp"] = snp
        if snp is not None and repeat is not None:
            channels["joint"] = snp_weight * snp + repeat_weight * repeat
        elif snp is not None and repeat_weight == 0:
            channels["joint"] = snp_weight * snp
        elif repeat is not None and snp_weight == 0:
            channels["joint"] = repeat_weight * repeat
        for channel, value in channels.items():
            result.setdefault(locus_id, {}).setdefault(taxon_id, {}).setdefault(
                channel, {}
            )[reference_id] = value
    return result


def _aggregate_taxon_distances(
    locus_values: Mapping[
        str, Mapping[str, Mapping[str, Mapping[str, float]]]
    ],
    selected_loci: Sequence[str],
    taxa: Sequence[str],
    k: int,
    locus_weights: Mapping[str, float] | None = None,
) -> tuple[dict[str, dict[str, float | None]], dict[str, set[str]]]:
    """Aggregate each reference across loci, then average the k nearest."""
    weights = locus_weights or {locus_id: 1.0 for locus_id in selected_loci}
    distances: dict[str, dict[str, float | None]] = {
        taxon_id: {channel: None for channel in CHANNELS} for taxon_id in taxa
    }
    nearest: dict[str, set[str]] = {taxon_id: set() for taxon_id in taxa}
    for taxon_id in taxa:
        for channel in CHANNELS:
            by_reference: dict[str, list[float]] = {}
            for locus_id in selected_loci:
                reference_values = (
                    locus_values.get(locus_id, {})
                    .get(taxon_id, {})
                    .get(channel, {})
                )
                for reference_id, value in reference_values.items():
                    by_reference.setdefault(reference_id, []).append(value)
            complete = sorted(
                (
                    (
                        sum(
                            value * max(0.0, float(weights.get(locus_id, 1.0)))
                            for locus_id, value in zip(selected_loci, values)
                        )
                        / sum(
                            max(0.0, float(weights.get(locus_id, 1.0)))
                            for locus_id in selected_loci
                        ),
                        reference_id,
                    )
                    for reference_id, values in by_reference.items()
                    if len(values) == len(selected_loci)
                    and sum(
                        max(0.0, float(weights.get(locus_id, 1.0)))
                        for locus_id in selected_loci
                    ) > 0
                ),
                key=lambda item: (item[0], item[1]),
            )
            selected = complete[: min(k, len(complete))]
            if selected:
                distances[taxon_id][channel] = sum(
                    value for value, _reference_id in selected
                ) / len(selected)
                nearest[taxon_id].update(
                    reference_id for _value, reference_id in selected
                )
    return distances, nearest


def _locus_taxon_distances(
    values: Mapping[str, Mapping[str, Mapping[str, Mapping[str, float]]]],
    loci: Sequence[str],
    taxa: Sequence[str],
    k: int,
) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    for locus_id in loci:
        for taxon_id in taxa:
            distance = _mean_k_smallest(
                values.get(locus_id, {}).get(taxon_id, {}).get("joint", {}).values(), k
            )
            if distance is not None:
                result.setdefault(locus_id, {})[taxon_id] = distance
    return result


def assign_best_taxon(
    *,
    sample_id: str,
    locus_marker_rows: Iterable[Mapping],
    reference_metadata: Mapping[str, Mapping[str, str]],
    expected_loci: int,
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    k: int = 3,
    minimum_loci: int = 3,
    minimum_locus_fraction: float = 0.8,
    minimum_relative_margin: float = 0.1,
    locus_weights: Mapping[str, float] | None = None,
    input_mode: str = "assembly",
    bootstrap_replicates: int = 200,
    minimum_bootstrap_support: float = 0.9,
    minimum_discriminative_loci: int = 2,
    maximum_compatible_distance: float = 1.0,
    locus_quality: Mapping[str, Mapping[str, object]] | None = None,
    seed: int = 0,
) -> AutomaticTaxonAssignment:
    """Select the nearest annotated taxon using existing MLVA marker distances.

    For every query locus represented in every candidate taxon, the existing
    normalized marker distance is ``snp_weight * SNP + repeat_weight * repeat``.
    Each reference distance is the mean over those loci, and a taxon's distance
    is the mean of its nearest ``k`` complete references.  The reported score is
    ``locus_recovery_fraction / (1 + taxon_distance)``.  Missing query loci are
    therefore a confidence penalty and never evidence that a taxon is absent.
    A multi-taxon result is SUPPORTED only when locus thresholds pass and
    ``(second_distance - best_distance) / max(second_distance, 1e-12)`` meets
    ``minimum_relative_margin``.
    """
    if expected_loci < 1:
        raise ValueError("expected_loci must be at least 1")
    if k < 1 or minimum_loci < 1:
        raise ValueError("k and minimum_loci must be at least 1")
    if not 0 <= minimum_locus_fraction <= 1:
        raise ValueError("minimum_locus_fraction must be between 0 and 1")
    if not 0 <= minimum_relative_margin <= 1:
        raise ValueError("minimum_relative_margin must be between 0 and 1")
    if snp_weight < 0 or repeat_weight < 0 or snp_weight + repeat_weight <= 0:
        raise ValueError("SNP and repeat weights must be non-negative with a positive total")

    reference_taxon, taxon_names = _reference_taxa(reference_metadata)
    taxa = sorted(set(reference_taxon.values()))
    if not taxa:
        raise ValueError("Reference metadata contains no taxon_id values")
    marker_rows = [dict(row) for row in locus_marker_rows]
    compared_reference_ids = {
        str(row.get("reference_id", "")).strip()
        for row in marker_rows
        if str(row.get("reference_id", "")).strip()
    }
    missing_metadata = sorted(compared_reference_ids - set(reference_taxon))
    if missing_metadata:
        raise ValueError(
            "Taxonomic metadata is missing for compared references: "
            + ", ".join(missing_metadata[:10])
        )
    if input_mode not in {"assembly", "fastq", "illumina"}:
        raise ValueError("input_mode must be assembly, fastq, or illumina")
    values = _locus_reference_values(
        marker_rows, reference_taxon, snp_weight, repeat_weight
    )
    informative = sorted(
        locus_id
        for locus_id, by_taxon in values.items()
        if all(bool(by_taxon.get(taxon_id, {}).get("joint")) for taxon_id in taxa)
    )
    recovery = min(1.0, len(informative) / expected_loci) if expected_loci else 0.0
    supplied_weights = locus_weights is not None
    weights = {
        locus_id: max(0.0, float((locus_weights or {}).get(locus_id, 1.0)))
        for locus_id in informative
    }
    discriminative = [locus_id for locus_id in informative if weights[locus_id] >= 0.05]
    scoring_loci = discriminative if discriminative else ([] if supplied_weights else informative)
    distances, nearest = _aggregate_taxon_distances(
        values, scoring_loci, taxa, k, weights
    )
    ranked = sorted(
        (
            (float(channels["joint"]), taxon_id)
            for taxon_id, channels in distances.items()
            if channels["joint"] is not None
        ),
        key=lambda item: (item[0], item[1]),
    )
    per_locus = _locus_taxon_distances(values, informative, taxa, k)
    bootstrap_wins = {taxon_id: 0 for taxon_id in taxa}
    valid_bootstraps = 0
    if scoring_loci and bootstrap_replicates > 0:
        rng = random.Random(seed)
        for _replicate in range(bootstrap_replicates):
            selected = [rng.choice(scoring_loci) for _ in scoring_loci]
            replicate_distances = _aggregate_bootstrap_joint_distances(
                values, selected, taxa, k, weights
            )
            ordered = sorted(
                (float(distance), taxon_id)
                for taxon_id, distance in replicate_distances.items()
                if distance is not None
            )
            if not ordered:
                continue
            valid_bootstraps += 1
            bootstrap_wins[ordered[0][1]] += 1
    evidence = []
    for rank, (distance, taxon_id) in enumerate(ranked, start=1):
        evidence.append(
            {
                "sample_id": sample_id,
                "taxon_id": taxon_id,
                "species": taxon_names.get(taxon_id, ""),
                "rank": rank,
                "score": recovery / (1.0 + distance),
                "distance": distance,
                "references_compared": len(nearest[taxon_id]),
                "informative_loci": len(informative),
                "locus_recovery_fraction": recovery,
                "is_best_taxon": "yes" if rank == 1 else "no",
                "compatibility": 1.0 / (1.0 + distance),
                "bootstrap_win_fraction": (
                    bootstrap_wins[taxon_id] / valid_bootstraps
                    if valid_bootstraps else ""
                ),
            }
        )

    required_loci = min(expected_loci, minimum_loci)
    fastq_mode = input_mode in {"fastq", "illumina"}
    required_fraction = max(minimum_locus_fraction, 0.9) if fastq_mode else minimum_locus_fraction
    required_margin = minimum_relative_margin * (1.5 if fastq_mode else 1.0)
    required_discriminative = max(minimum_discriminative_loci, 3 if fastq_mode else 2)
    best_distance = ranked[0][0] if ranked else None
    second_distance = ranked[1][0] if len(ranked) > 1 else None
    absolute_margin = (
        second_distance - best_distance
        if best_distance is not None and second_distance is not None else None
    )
    relative_margin = (
        absolute_margin / max(second_distance, 1e-12)
        if absolute_margin is not None and second_distance is not None else None
    )
    best_taxon_id = ranked[0][1] if ranked else ""
    bootstrap_support = (
        bootstrap_wins.get(best_taxon_id, 0) / valid_bootstraps
        if valid_bootstraps else None
    )
    locus_rows = []
    supporting = conflicting = 0
    for locus_id in sorted(set(values) | set((locus_quality or {}))):
        ordered = sorted((distance, taxon_id) for taxon_id, distance in per_locus.get(locus_id, {}).items())
        favored_id = (
            ordered[0][1]
            if ordered and not (
                len(ordered) > 1
                and math.isclose(ordered[0][0], ordered[1][0], rel_tol=1e-12, abs_tol=1e-12)
            )
            else ""
        )
        supports = bool(best_taxon_id and favored_id == best_taxon_id)
        conflicts = bool(best_taxon_id and favored_id and favored_id != best_taxon_id)
        if weights.get(locus_id, 0.0) >= 0.05:
            supporting += int(supports)
            conflicting += int(conflicts)
        quality = dict((locus_quality or {}).get(locus_id, {}))
        locus_rows.append({
            "sample_id": sample_id,
            "locus_id": locus_id,
            "recovered": "yes" if locus_id in informative else "no",
            "taxonomic_weight": weights.get(locus_id, 0.0),
            "discriminative": "yes" if weights.get(locus_id, 0.0) >= 0.05 else "no",
            "favored_taxon_id": favored_id,
            "favored_taxon": taxon_names.get(favored_id, favored_id),
            "closest_taxon_distance": ordered[0][0] if ordered else "",
            "runner_up_distance": ordered[1][0] if len(ordered) > 1 else "",
            "distance_margin": ordered[1][0] - ordered[0][0] if len(ordered) > 1 else "",
            "supports_assignment": "yes" if supports else "no",
            "conflicts_with_assignment": "yes" if conflicts else "no",
            "depth": quality.get("depth", ""),
            "consensus_strength": quality.get("consensus_strength", ""),
            "quality_status": quality.get("status", ""),
        })

    low_quality_discriminative = sum(
        1
        for locus_id in discriminative
        if str((locus_quality or {}).get(locus_id, {}).get("status", "")).upper()
        not in {"", "PASS", "PRESENT", "INCLUDED"}
        or (
            _finite((locus_quality or {}).get(locus_id, {}).get("consensus_strength"))
            is not None
            and float((locus_quality or {})[locus_id]["consensus_strength"]) < 0.75
        )
    )
    if len(informative) < required_loci or recovery < required_fraction:
        status = "INSUFFICIENT_EVIDENCE"
        reason = "INSUFFICIENT_INFORMATIVE_LOCI"
    elif len(discriminative) < required_discriminative:
        status = "AMBIGUOUS"
        reason = "INSUFFICIENT_DISCRIMINATIVE_LOCI"
    elif len(taxa) == 1:
        status = "AMBIGUOUS"
        reason = "NO_ALTERNATIVE_TAXON_IN_REFERENCE_PANEL"
    elif len(ranked) < 2:
        status = "INSUFFICIENT_EVIDENCE"
        reason = "FEWER_THAN_TWO_COMPARABLE_TAXA"
    elif fastq_mode and low_quality_discriminative > max(0, len(discriminative) // 4):
        status = "AMBIGUOUS"
        reason = "FASTQ_LOCUS_QUALITY_INSUFFICIENT"
    else:
        assert best_distance is not None and relative_margin is not None
        if best_distance > maximum_compatible_distance:
            status = "INSUFFICIENT_EVIDENCE"
            reason = "CLOSEST_TAXON_NOT_COMPATIBLE"
        elif relative_margin + 1e-12 < required_margin:
            status = "AMBIGUOUS"
            reason = "TAXON_DISTANCE_MARGIN_BELOW_THRESHOLD"
        elif bootstrap_replicates > 0 and (
            bootstrap_support is None
            or bootstrap_support < minimum_bootstrap_support
        ):
            status = "AMBIGUOUS"
            reason = "BOOTSTRAP_SUPPORT_BELOW_THRESHOLD"
        elif conflicting > supporting:
            status = "AMBIGUOUS"
            reason = "CONFLICTING_LOCUS_EVIDENCE"
        else:
            status = "SUPPORTED"
            reason = "NEAREST_TAXON_SEPARATED"

    best = evidence[0] if evidence else {}
    second = evidence[1] if len(evidence) > 1 else {}
    best_score = _finite(best.get("score"))
    second_score = _finite(second.get("score"))
    if status == "SUPPORTED":
        assignment_name = str(best.get("species") or best.get("taxon_id", ""))
        assignment_id = str(best.get("taxon_id", ""))
        assignment_rank = "species"
        assignment_status = "SPECIES_ASSIGNED"
        confidence = "HIGH" if bootstrap_support is not None and bootstrap_support >= 0.95 and not conflicting else "MODERATE"
    elif best:
        assignment_name = str(best.get("species") or best.get("taxon_id", ""))
        assignment_id = str(best.get("taxon_id", ""))
        assignment_rank = "species"
        assignment_status = "CLOSEST_TAXON_LOW_CONFIDENCE"
        confidence = "LOW"
    else:
        assignment_name = "Unresolved"
        assignment_rank = "unresolved"
        assignment_id = ""
        assignment_status = "UNRESOLVED"
        confidence = "UNRESOLVED"
    summary = {
        "sample_id": sample_id,
        "best_taxon": best.get("taxon_id", ""),
        "best_species": best.get("species", ""),
        "taxon_score": "" if best_score is None else best_score,
        "second_best_taxon": second.get("taxon_id", ""),
        "second_best_score": "" if second_score is None else second_score,
        "score_margin": ""
        if best_score is None or second_score is None
        else best_score - second_score,
        "informative_loci": len(informative),
        "expected_loci": expected_loci,
        "locus_recovery_fraction": recovery,
        "taxonomic_status": status,
        "status_reason": reason,
        "assignment": assignment_name,
        "assignment_taxon_id": assignment_id,
        "assignment_rank": assignment_rank,
        "assignment_status": assignment_status,
        "confidence": confidence,
        "closest_taxon": best.get("species") or best.get("taxon_id", ""),
        "closest_taxon_id": best.get("taxon_id", ""),
        "closest_distance": "" if best_distance is None else best_distance,
        "runner_up_taxon": second.get("species") or second.get("taxon_id", ""),
        "runner_up_taxon_id": second.get("taxon_id", ""),
        "runner_up_distance": "" if second_distance is None else second_distance,
        "distance_margin": "" if absolute_margin is None else absolute_margin,
        "relative_margin": "" if relative_margin is None else relative_margin,
        "loci_recovered": len(informative),
        "informative_loci_recovered": len(informative),
        "discriminative_loci_recovered": len(discriminative),
        "loci_supporting_assignment": supporting,
        "conflicting_loci": conflicting,
        "bootstrap_support": "" if bootstrap_support is None else bootstrap_support,
        "bootstrap_replicates": valid_bootstraps,
        "input_mode": input_mode,
    }
    LOGGER.debug(
        "Taxon candidates=%s winner=%s distance=%s runner_up=%s margin=%s bootstrap=%s support_loci=%s conflicts=%s status=%s backoff=%s",
        [(row["taxon_id"], row["distance"]) for row in evidence], best_taxon_id,
        best_distance, second.get("taxon_id", ""), relative_margin, bootstrap_support,
        supporting, conflicting, status, assignment_rank,
    )
    return AutomaticTaxonAssignment(summary, tuple(evidence), tuple(locus_rows))


def _aggregate_bootstrap_joint_distances(
    locus_values: Mapping[
        str, Mapping[str, Mapping[str, Mapping[str, float]]]
    ],
    selected_loci: Sequence[str],
    taxa: Sequence[str],
    k: int,
    locus_weights: Mapping[str, float] | None = None,
) -> dict[str, float | None]:
    """Aggregate only joint distances for a bootstrap replicate.

    Bootstrap decisions do not use channel-specific distances or nearest-reference
    identities. Keeping this hot path separate avoids rebuilding those unused
    structures for every replicate while retaining reference-across-loci averaging.
    """
    distances: dict[str, float | None] = {}
    locus_count = len(selected_loci)
    weights = locus_weights or {}
    for taxon_id in taxa:
        by_reference: dict[str, list[tuple[str, float]]] = {}
        for locus_id in selected_loci:
            reference_values = (
                locus_values.get(locus_id, {})
                .get(taxon_id, {})
                .get("joint", {})
            )
            for reference_id, value in reference_values.items():
                by_reference.setdefault(reference_id, []).append((locus_id, value))
        complete = sorted(
            sum(value * weights.get(locus_id, 1.0) for locus_id, value in values)
            / sum(weights.get(locus_id, 1.0) for locus_id, _value in values)
            for values in by_reference.values()
            if len(values) == locus_count
            and sum(weights.get(locus_id, 1.0) for locus_id, _value in values) > 0
        )
        selected = complete[: min(k, len(complete))]
        distances[taxon_id] = (
            sum(selected) / len(selected) if selected else None
        )
    return distances


def _nonconformity(
    distances: Mapping[str, Mapping[str, float | None]],
    taxon_id: str,
    channel: str,
    epsilon: float = 1e-12,
) -> float | None:
    within = distances[taxon_id][channel]
    outside_values = [
        values[channel]
        for other_taxon, values in distances.items()
        if other_taxon != taxon_id and values[channel] is not None
    ]
    outside = min(outside_values) if outside_values else None
    if within is None or outside is None:
        return None
    return float(within) / max(float(outside), epsilon)


def build_taxon_calibration(
    *,
    reference_locus_rows: Iterable[Mapping],
    reference_metadata: Mapping[str, Mapping[str, str]],
    panel_sha256: str,
    database_signature: str,
    k: int = 3,
    alpha: float = 0.05,
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    minimum_loci: int = 3,
) -> tuple[TaxonCalibration, list[dict]]:
    """Build label-conditional scores from leave-one-reference-out distances.

    ``reference_locus_rows`` must contain one row per query reference, locus,
    and candidate reference. Required columns are ``query_reference_id``,
    ``reference_id``, ``locus_id``, ``normalized_repeat_distance``, and
    ``normalized_snp_distance``. The query reference is always excluded from
    its own neighborhood.
    """
    if k < 1:
        raise ValueError("Calibration k must be at least 1")
    reference_taxon, _taxon_names = _reference_taxa(reference_metadata)
    taxa = sorted(set(reference_taxon.values()))
    if len(taxa) < 2:
        raise ValueError("Calibration requires at least two labeled taxa")
    rows = [dict(row) for row in reference_locus_rows]
    query_ids = sorted(
        {
            str(row.get("query_reference_id", ""))
            for row in rows
            if str(row.get("query_reference_id", ""))
        }
    )
    missing_metadata = sorted(set(query_ids) - set(reference_taxon))
    if missing_metadata:
        raise ValueError(
            "Calibration queries lack taxon metadata: " + ", ".join(missing_metadata)
        )

    scores: dict[str, dict[str, list[float]]] = {
        channel: {taxon_id: [] for taxon_id in taxa} for channel in CHANNELS
    }
    score_rows: list[dict] = []
    for query_id in query_ids:
        query_taxon = reference_taxon[query_id]
        query_rows = [
            row
            for row in rows
            if str(row.get("query_reference_id", "")) == query_id
            and str(row.get("reference_id", "")) != query_id
        ]
        locus_values = _locus_reference_values(
            query_rows, reference_taxon, snp_weight, repeat_weight
        )
        callable_loci = sorted(
            locus_id
            for locus_id, by_taxon in locus_values.items()
            if all(bool(by_taxon.get(taxon_id, {}).get("joint")) for taxon_id in taxa)
        )
        if len(callable_loci) < minimum_loci:
            continue
        distances, _nearest = _aggregate_taxon_distances(
            locus_values, callable_loci, taxa, k
        )
        for channel in CHANNELS:
            within = distances[query_taxon][channel]
            outside_values = [
                values[channel]
                for taxon_id, values in distances.items()
                if taxon_id != query_taxon and values[channel] is not None
            ]
            outside = min(outside_values) if outside_values else None
            score = _nonconformity(distances, query_taxon, channel)
            if within is None or outside is None or score is None:
                continue
            scores[channel][query_taxon].append(score)
            score_rows.append(
                {
                    "reference_id": query_id,
                    "taxon_id": query_taxon,
                    "channel": channel,
                    "within_distance": within,
                    "outside_distance": outside,
                    "nonconformity": score,
                }
            )

    missing_scores = [
        f"{channel}/{taxon_id}"
        for channel in CHANNELS
        for taxon_id in taxa
        if not scores[channel][taxon_id]
    ]
    if missing_scores:
        raise ValueError(
            "No valid leave-one-out calibration scores for: "
            + ", ".join(missing_scores)
        )
    taxon_counts = {
        taxon_id: sum(reference_taxon.get(query_id) == taxon_id for query_id in query_ids)
        for taxon_id in taxa
    }
    calibration = TaxonCalibration(
        panel_sha256=panel_sha256,
        database_signature=database_signature,
        taxon_counts=taxon_counts,
        scores={
            channel: {
                taxon_id: tuple(values)
                for taxon_id, values in by_taxon.items()
            }
            for channel, by_taxon in scores.items()
        },
        k=k,
        alpha=alpha,
        snp_weight=snp_weight,
        repeat_weight=repeat_weight,
        minimum_loci=minimum_loci,
    )
    return calibration, score_rows


def run_taxon_calibration(
    *,
    reference_distances_path: str | Path,
    reference_metadata_path: str | Path,
    sequence_index_path: str | Path,
    outdir: str | Path,
    k: int = 3,
    alpha: float = 0.05,
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    minimum_loci: int = 3,
) -> dict[str, Path]:
    """Write a calibration artifact from audited leave-one-out distance rows."""
    distances_path = Path(reference_distances_path)
    metadata_path = Path(reference_metadata_path)
    index_path = Path(sequence_index_path)
    for label, path in (
        ("reference distances", distances_path),
        ("reference metadata", metadata_path),
        ("reference sequence index", index_path),
    ):
        if not path.is_file():
            raise ValueError(f"{label} path does not exist: {path}")

    with distances_path.open(newline="") as handle:
        distance_rows = list(csv.DictReader(handle, delimiter="\t"))
    required_distance_fields = {
        "query_reference_id",
        "reference_id",
        "locus_id",
        "normalized_repeat_distance",
        "normalized_snp_distance",
    }
    distance_fields = set(distance_rows[0]) if distance_rows else set()
    if not required_distance_fields.issubset(distance_fields):
        raise ValueError(
            "Reference distances require columns: "
            + ", ".join(sorted(required_distance_fields))
        )

    with metadata_path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if "\t" in (sample.splitlines() or [""])[0] else ","
        metadata_rows = list(csv.DictReader(handle, delimiter=delimiter))
    if not metadata_rows or not {"reference_id", "taxon_id"}.issubset(
        metadata_rows[0]
    ):
        raise ValueError("Reference metadata requires reference_id and taxon_id columns")
    metadata = {
        str(row["reference_id"]): {
            str(key): str(value or "") for key, value in row.items()
        }
        for row in metadata_rows
        if str(row.get("reference_id", "")).strip()
    }

    with index_path.open(newline="") as handle:
        index_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not index_rows:
        raise ValueError("Reference sequence index contains no rows")
    panel_sha256 = str(index_rows[0].get("panel_sha256", ""))
    database_signature = str(index_rows[0].get("database_signature", ""))
    if not panel_sha256 or not database_signature:
        raise ValueError(
            "Reference sequence index lacks panel_sha256 or database_signature"
        )
    if any(
        str(row.get("panel_sha256", "")) != panel_sha256
        or str(row.get("database_signature", "")) != database_signature
        for row in index_rows
    ):
        raise ValueError("Reference sequence index contains inconsistent signatures")

    calibration, score_rows = build_taxon_calibration(
        reference_locus_rows=distance_rows,
        reference_metadata=metadata,
        panel_sha256=panel_sha256,
        database_signature=database_signature,
        k=k,
        alpha=alpha,
        snp_weight=snp_weight,
        repeat_weight=repeat_weight,
        minimum_loci=minimum_loci,
    )
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    calibration_path = calibration.write(output / "taxon_calibration.json")
    scores_path = output / "taxon_calibration_scores.tsv"
    with scores_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=CALIBRATION_SCORE_FIELDS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(score_rows)
    return {"calibration": calibration_path, "scores": scores_path}


def assign_target_taxon(
    *,
    sample_id: str,
    target_taxon_id: str,
    locus_marker_rows: Iterable[Mapping],
    placement_rows: Iterable[Mapping],
    reference_metadata: Mapping[str, Mapping[str, str]],
    calibration: TaxonCalibration,
    alpha: float | None = None,
    min_loci: int | None = None,
    min_locus_fraction: float = 0.8,
    bootstrap_replicates: int = 2000,
    min_bootstrap_support: float = 0.95,
    seed: int = 12345,
    max_mean_placement_entropy: float | None = None,
    min_median_placement_lwr: float | None = None,
    expected_loci: int | None = None,
) -> TaxonAssignment:
    """Assign a requested target using only MLVA repeat and marker-SNP evidence."""
    if not 0 < min_locus_fraction <= 1:
        raise ValueError("min_locus_fraction must be greater than 0 and at most 1")
    if bootstrap_replicates < 0:
        raise ValueError("bootstrap_replicates cannot be negative")
    if not 0 <= min_bootstrap_support <= 1:
        raise ValueError("min_bootstrap_support must be between 0 and 1")
    effective_alpha = calibration.alpha if alpha is None else alpha
    if not 0 < effective_alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    effective_min_loci = calibration.minimum_loci if min_loci is None else min_loci
    if effective_min_loci < 1:
        raise ValueError("min_loci must be at least 1")

    reference_taxon, taxon_names = _reference_taxa(reference_metadata)
    if target_taxon_id not in set(reference_taxon.values()):
        raise ValueError(
            f"Target taxon {target_taxon_id!r} is absent from reference metadata"
        )
    taxa = sorted(set(reference_taxon.values()))
    if len(taxa) < 2:
        raise ValueError("Taxon assignment requires at least two labeled taxa")
    missing_calibration = [
        f"{channel}/{taxon_id}"
        for channel in CHANNELS
        for taxon_id in taxa
        if not calibration.scores.get(channel, {}).get(taxon_id)
    ]
    if missing_calibration:
        raise ValueError(
            "Taxon calibration lacks scores for: " + ", ".join(missing_calibration)
        )

    marker_rows = [dict(row) for row in locus_marker_rows]
    placements = [dict(row) for row in placement_rows]
    locus_values = _locus_reference_values(
        marker_rows,
        reference_taxon,
        calibration.snp_weight,
        calibration.repeat_weight,
    )
    callable_loci = sorted(
        locus_id
        for locus_id, by_taxon in locus_values.items()
        if all(
            bool(by_taxon.get(taxon_id, {}).get("joint"))
            for taxon_id in taxa
        )
    )
    expected = expected_loci if expected_loci is not None else len(locus_values)
    required_by_fraction = math.ceil(expected * min_locus_fraction) if expected else 0
    required_loci = max(effective_min_loci, required_by_fraction)

    distances, nearest = _aggregate_taxon_distances(
        locus_values, callable_loci, taxa, calibration.k
    )
    candidate_rows: list[dict] = []
    for taxon_id in taxa:
        nonconformity = {
            channel: _nonconformity(distances, taxon_id, channel)
            for channel in CHANNELS
        }
        p_values = {
            channel: (
                conformal_p_value(
                    float(nonconformity[channel]),
                    calibration.scores[channel][taxon_id],
                )
                if nonconformity[channel] is not None
                else None
            )
            for channel in CHANNELS
        }
        candidate_rows.append(
            {
                "sample_id": sample_id,
                "taxon_id": taxon_id,
                "taxon_name": taxon_names.get(taxon_id, ""),
                "reference_count": sum(
                    assigned_taxon == taxon_id
                    for assigned_taxon in reference_taxon.values()
                ),
                "repeat_distance": distances[taxon_id]["repeat"],
                "snp_distance": distances[taxon_id]["snp"],
                "joint_distance": distances[taxon_id]["joint"],
                "repeat_nonconformity": nonconformity["repeat"],
                "snp_nonconformity": nonconformity["snp"],
                "joint_nonconformity": nonconformity["joint"],
                "repeat_p_value": p_values["repeat"],
                "snp_p_value": p_values["snp"],
                "joint_p_value": p_values["joint"],
                "repeat_accepted": (
                    p_values["repeat"] is not None
                    and p_values["repeat"] > effective_alpha
                ),
                "snp_accepted": (
                    p_values["snp"] is not None
                    and p_values["snp"] > effective_alpha
                ),
                "joint_accepted": (
                    p_values["joint"] is not None
                    and p_values["joint"] > effective_alpha
                ),
                "nearest_references": ",".join(sorted(nearest[taxon_id])),
            }
        )
    candidate_rows.sort(
        key=lambda row: (
            float(row["joint_distance"])
            if row["joint_distance"] is not None
            else float("inf"),
            str(row["taxon_id"]),
        )
    )
    for rank, row in enumerate(candidate_rows, start=1):
        row["rank"] = rank

    target = next(row for row in candidate_rows if row["taxon_id"] == target_taxon_id)
    alternatives = [row for row in candidate_rows if row["taxon_id"] != target_taxon_id]
    best_alternative = alternatives[0]
    prediction_set = tuple(
        str(row["taxon_id"]) for row in candidate_rows if row["joint_accepted"]
    )

    locus_rows: list[dict] = []
    target_favoring = alternative_favoring = unresolved = 0
    entropy_by_locus: dict[str, list[float]] = {}
    lwr_by_locus: dict[str, list[float]] = {}
    for row in placements:
        locus_id = str(row.get("locus_id", ""))
        entropy = _finite(row.get("placement_entropy"))
        lwr = _finite(row.get("like_weight_ratio"))
        if entropy is not None:
            entropy_by_locus.setdefault(locus_id, []).append(entropy)
        if lwr is not None:
            lwr_by_locus.setdefault(locus_id, []).append(lwr)
    for locus_id in callable_loci:
        target_values = locus_values[locus_id][target_taxon_id]
        alternative_values = locus_values[locus_id][str(best_alternative["taxon_id"])]
        target_repeat = _mean_k_smallest(
            target_values["repeat"].values(), calibration.k
        )
        alternative_repeat = _mean_k_smallest(
            alternative_values["repeat"].values(), calibration.k
        )
        target_snp = _mean_k_smallest(target_values["snp"].values(), calibration.k)
        alternative_snp = _mean_k_smallest(
            alternative_values["snp"].values(), calibration.k
        )
        target_joint = float(
            _mean_k_smallest(target_values["joint"].values(), calibration.k)
        )
        alternative_joint = float(
            _mean_k_smallest(alternative_values["joint"].values(), calibration.k)
        )
        margin = alternative_joint - target_joint
        if math.isclose(margin, 0.0, rel_tol=1e-12, abs_tol=1e-12):
            interpretation = "UNRESOLVED"
            unresolved += 1
        elif margin > 0:
            interpretation = "TARGET_FAVORED"
            target_favoring += 1
        else:
            interpretation = "ALTERNATIVE_FAVORED"
            alternative_favoring += 1
        locus_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "target_taxon_id": target_taxon_id,
                "alternative_taxon_id": best_alternative["taxon_id"],
                "target_repeat_distance": target_repeat,
                "alternative_repeat_distance": alternative_repeat,
                "target_snp_distance": target_snp,
                "alternative_snp_distance": alternative_snp,
                "target_joint_distance": target_joint,
                "alternative_joint_distance": alternative_joint,
                "joint_distance_margin": margin,
                "placement_entropy": (
                    statistics.mean(entropy_by_locus.get(locus_id, []))
                    if entropy_by_locus.get(locus_id)
                    else None
                ),
                "best_placement_lwr": (
                    max(lwr_by_locus.get(locus_id, []))
                    if lwr_by_locus.get(locus_id)
                    else None
                ),
                "interpretation": interpretation,
            }
        )

    rng = random.Random(seed)
    target_wins = alternative_wins = ties = 0
    bootstrap_margins: list[float] = []
    if callable_loci and bootstrap_replicates:
        for _replicate in range(bootstrap_replicates):
            selected = [rng.choice(callable_loci) for _ in callable_loci]
            replicate_distances = _aggregate_bootstrap_joint_distances(
                locus_values, selected, taxa, calibration.k
            )
            ordered = sorted(
                (
                    (float(distance), taxon_id)
                    for taxon_id, distance in replicate_distances.items()
                    if distance is not None
                ),
                key=lambda item: (item[0], item[1]),
            )
            if len(ordered) < 2:
                continue
            target_distance = float(replicate_distances[target_taxon_id])
            alternative_distance = min(
                distance for distance, taxon_id in ordered if taxon_id != target_taxon_id
            )
            margin = alternative_distance - target_distance
            bootstrap_margins.append(margin)
            if math.isclose(margin, 0.0, rel_tol=1e-12, abs_tol=1e-12):
                ties += 1
            elif margin > 0:
                target_wins += 1
            else:
                alternative_wins += 1
    valid_bootstraps = target_wins + alternative_wins + ties
    bootstrap_support = (
        target_wins / valid_bootstraps if valid_bootstraps else None
    )
    bootstrap_tie_fraction = ties / valid_bootstraps if valid_bootstraps else None

    all_entropies = [
        statistics.mean(values) for values in entropy_by_locus.values() if values
    ]
    all_lwrs = [max(values) for values in lwr_by_locus.values() if values]
    mean_entropy = statistics.mean(all_entropies) if all_entropies else None
    median_lwr = statistics.median(all_lwrs) if all_lwrs else None

    qc_flags: list[str] = []
    if len(callable_loci) < required_loci:
        qc_flags.append("INSUFFICIENT_CALLABLE_LOCI")
    if (
        max_mean_placement_entropy is not None
        and mean_entropy is not None
        and mean_entropy > max_mean_placement_entropy
    ):
        qc_flags.append("HIGH_PLACEMENT_ENTROPY")
    if (
        min_median_placement_lwr is not None
        and median_lwr is not None
        and median_lwr < min_median_placement_lwr
    ):
        qc_flags.append("LOW_PLACEMENT_LWR")
    target_count = calibration.taxon_counts.get(target_taxon_id, 0)
    if target_count < max(2, calibration.k + 1):
        qc_flags.append("SMALL_TARGET_CALIBRATION_COHORT")

    target_channels_accepted = bool(target["repeat_accepted"]) and bool(
        target["snp_accepted"]
    )
    channels_concordant = bool(target["repeat_accepted"]) == bool(
        target["snp_accepted"]
    )
    bootstrap_passes = (
        bootstrap_support is not None
        and bootstrap_support >= min_bootstrap_support
    )
    if qc_flags:
        decision = "INDETERMINATE"
        reason = qc_flags[0]
    elif not channels_concordant:
        decision = "INDETERMINATE"
        reason = "REPEAT_SNP_EVIDENCE_DISAGREES"
    elif (
        prediction_set == (target_taxon_id,)
        and target_channels_accepted
        and bootstrap_passes
    ):
        decision = "POSITIVE"
        reason = "TARGET_UNIQUELY_SUPPORTED"
    elif target_taxon_id not in prediction_set and prediction_set:
        decision = "NEGATIVE"
        reason = "TARGET_EXCLUDED_ALTERNATIVE_SUPPORTED"
    elif not prediction_set:
        decision = "INDETERMINATE"
        reason = "OUTSIDE_CALIBRATED_REFERENCE_SPACE"
    elif len(prediction_set) > 1:
        decision = "INDETERMINATE"
        reason = "MULTIPLE_TAXA_COMPATIBLE"
    else:
        decision = "INDETERMINATE"
        reason = "BOOTSTRAP_SUPPORT_BELOW_THRESHOLD"

    target_distance = _finite(target["joint_distance"])
    alternative_distance = _finite(best_alternative["joint_distance"])
    distance_margin = (
        alternative_distance - target_distance
        if target_distance is not None and alternative_distance is not None
        else None
    )
    summary = {
        "sample_id": sample_id,
        "target_taxon_id": target_taxon_id,
        "target_taxon_name": taxon_names.get(target_taxon_id, ""),
        "decision": decision,
        "decision_reason": reason,
        "prediction_set": ",".join(prediction_set),
        "target_joint_p_value": target["joint_p_value"],
        "target_snp_p_value": target["snp_p_value"],
        "target_repeat_p_value": target["repeat_p_value"],
        "target_bootstrap_support": bootstrap_support,
        "bootstrap_tie_fraction": bootstrap_tie_fraction,
        "target_distance": target_distance,
        "best_alternative_taxon_id": best_alternative["taxon_id"],
        "best_alternative_taxon_name": best_alternative["taxon_name"],
        "best_alternative_joint_p_value": best_alternative["joint_p_value"],
        "best_alternative_distance": alternative_distance,
        "distance_margin": distance_margin,
        "bootstrap_margin_lower": _percentile(bootstrap_margins, 0.025),
        "bootstrap_margin_upper": _percentile(bootstrap_margins, 0.975),
        "callable_loci": len(callable_loci),
        "placed_loci": len(entropy_by_locus),
        "target_favoring_loci": target_favoring,
        "alternative_favoring_loci": alternative_favoring,
        "unresolved_loci": unresolved,
        "mean_placement_entropy": mean_entropy,
        "median_best_lwr": median_lwr,
        "qc_status": "PASS" if not qc_flags else "WARN",
        "qc_flags": ",".join(qc_flags),
        "calibration_version": calibration.schema_version,
    }
    return TaxonAssignment(
        summary=summary,
        candidates=tuple(candidate_rows),
        loci=tuple(locus_rows),
    )
