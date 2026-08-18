from __future__ import annotations

import csv
import json
import math
import random
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


CALIBRATION_SCHEMA_VERSION = "1.0"
CHANNELS = ("repeat", "snp", "joint")

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
        taxon_id = str(row.get("taxon_id", "")).strip()
        if not taxon_id:
            continue
        reference_taxon[str(reference_id)] = taxon_id
        name = str(row.get("taxon_name") or row.get("organism_name") or "").strip()
        if name:
            previous = taxon_names.get(taxon_id)
            if previous is not None and previous != name:
                raise ValueError(
                    f"Taxon {taxon_id!r} has conflicting names {previous!r} and {name!r}"
                )
            taxon_names[taxon_id] = name
    return reference_taxon, taxon_names


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
) -> tuple[dict[str, dict[str, float | None]], dict[str, set[str]]]:
    """Aggregate each reference across loci, then average the k nearest."""
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
                    (sum(values) / len(values), reference_id)
                    for reference_id, values in by_reference.items()
                    if len(values) == len(selected_loci)
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


def _aggregate_bootstrap_joint_distances(
    locus_values: Mapping[
        str, Mapping[str, Mapping[str, Mapping[str, float]]]
    ],
    selected_loci: Sequence[str],
    taxa: Sequence[str],
    k: int,
) -> dict[str, float | None]:
    """Aggregate only joint distances for a bootstrap replicate.

    Bootstrap decisions do not use channel-specific distances or nearest-reference
    identities. Keeping this hot path separate avoids rebuilding those unused
    structures for every replicate while retaining reference-across-loci averaging.
    """
    distances: dict[str, float | None] = {}
    locus_count = len(selected_loci)
    for taxon_id in taxa:
        by_reference: dict[str, list[float]] = {}
        for locus_id in selected_loci:
            reference_values = (
                locus_values.get(locus_id, {})
                .get(taxon_id, {})
                .get("joint", {})
            )
            for reference_id, value in reference_values.items():
                by_reference.setdefault(reference_id, []).append(value)
        complete = sorted(
            sum(values) / locus_count
            for values in by_reference.values()
            if len(values) == locus_count
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