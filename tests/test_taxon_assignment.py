from __future__ import annotations

import csv

import pytest

from mlvamaps.phylogeny import _add_taxon_assignment_outputs
from mlvamaps.taxon_assignment import (
    TaxonCalibration,
    _aggregate_bootstrap_joint_distances,
    _aggregate_taxon_distances,
    assign_target_taxon,
    build_taxon_calibration,
    conformal_p_value,
    run_taxon_calibration,
)


def test_bootstrap_joint_aggregation_matches_full_aggregation():
    taxa = ["target", "neighbor"]
    locus_values = {
        "L1": {
            "target": {
                "repeat": {"T1": 0.1, "T2": 0.2},
                "snp": {"T1": 0.3, "T2": 0.4},
                "joint": {"T1": 0.4, "T2": 0.6},
            },
            "neighbor": {
                "repeat": {"N1": 0.8, "N2": 0.9},
                "snp": {"N1": 1.0, "N2": 1.1},
                "joint": {"N1": 1.8, "N2": 2.0},
            },
        },
        "L2": {
            "target": {
                "repeat": {"T1": 0.2, "T2": 0.3},
                "snp": {"T1": 0.4, "T2": 0.5},
                "joint": {"T1": 0.6, "T2": 0.8},
            },
            "neighbor": {
                "repeat": {"N1": 0.7},
                "snp": {"N1": 0.9},
                "joint": {"N1": 1.6},
            },
        },
    }
    selected_loci = ["L2", "L1", "L2"]

    full, _nearest = _aggregate_taxon_distances(
        locus_values, selected_loci, taxa, k=2
    )
    bootstrap = _aggregate_bootstrap_joint_distances(
        locus_values, selected_loci, taxa, k=2
    )

    assert bootstrap == {
        taxon_id: full[taxon_id]["joint"] for taxon_id in taxa
    }


def _calibration(minimum_loci: int = 3) -> TaxonCalibration:
    compatible_scores = tuple([1.0] * 20)
    return TaxonCalibration(
        panel_sha256="panel",
        database_signature="database",
        taxon_counts={"target": 10, "neighbor": 10},
        scores={
            channel: {
                "target": compatible_scores,
                "neighbor": compatible_scores,
            }
            for channel in ("repeat", "snp", "joint")
        },
        k=1,
        alpha=0.05,
        minimum_loci=minimum_loci,
    )


def _metadata() -> dict[str, dict[str, str]]:
    return {
        "T1": {"taxon_id": "target", "taxon_name": "Target taxon"},
        "T2": {"taxon_id": "target", "taxon_name": "Target taxon"},
        "N1": {"taxon_id": "neighbor", "taxon_name": "Neighbor taxon"},
        "N2": {"taxon_id": "neighbor", "taxon_name": "Neighbor taxon"},
    }


def _rows(
    target_repeat: float,
    target_snp: float,
    neighbor_repeat: float,
    neighbor_snp: float,
    loci: int = 3,
) -> list[dict]:
    rows = []
    for locus_number in range(1, loci + 1):
        for reference_id, repeat, snp in (
            ("T1", target_repeat, target_snp),
            ("T2", target_repeat * 1.1, target_snp * 1.1),
            ("N1", neighbor_repeat, neighbor_snp),
            ("N2", neighbor_repeat * 1.1, neighbor_snp * 1.1),
        ):
            rows.append(
                {
                    "locus_id": f"L{locus_number}",
                    "reference_id": reference_id,
                    "normalized_repeat_distance": repeat,
                    "normalized_snp_distance": snp,
                }
            )
    return rows


def _placements(
    loci: int = 3, entropy: float = 0.05, lwr: float = 0.95
) -> list[dict]:
    return [
        {
            "locus_id": f"L{locus_number}",
            "reference_id": "T1",
            "placement_entropy": entropy,
            "like_weight_ratio": lwr,
        }
        for locus_number in range(1, loci + 1)
    ]


def _assign(rows: list[dict], **kwargs):
    return assign_target_taxon(
        sample_id="sample",
        target_taxon_id="target",
        locus_marker_rows=rows,
        placement_rows=_placements(
            len({str(row["locus_id"]) for row in rows})
        ),
        reference_metadata=_metadata(),
        calibration=_calibration(),
        min_locus_fraction=1.0,
        bootstrap_replicates=200,
        seed=7,
        **kwargs,
    )


def test_conformal_p_value_has_finite_sample_correction():
    assert conformal_p_value(0.5, [0.1, 0.2, 0.6, 0.7]) == pytest.approx(0.6)


def test_unique_target_support_is_positive_and_reproducible():
    first = _assign(_rows(0.1, 0.1, 1.0, 1.0))
    second = _assign(_rows(0.1, 0.1, 1.0, 1.0))

    assert first.summary["decision"] == "POSITIVE"
    assert first.summary["decision_reason"] == "TARGET_UNIQUELY_SUPPORTED"
    assert first.summary["prediction_set"] == "target"
    assert first.summary["target_bootstrap_support"] == 1.0
    assert first.summary == second.summary
    assert [row["interpretation"] for row in first.loci] == [
        "TARGET_FAVORED"
    ] * 3


def test_supported_neighbor_excludes_target():
    result = _assign(_rows(1.0, 1.0, 0.1, 0.1))

    assert result.summary["decision"] == "NEGATIVE"
    assert result.summary["decision_reason"] == "TARGET_EXCLUDED_ALTERNATIVE_SUPPORTED"
    assert result.summary["prediction_set"] == "neighbor"
    assert result.summary["target_bootstrap_support"] == 0.0


def test_shared_marker_evidence_is_indeterminate():
    result = _assign(_rows(0.1, 0.1, 0.1, 0.1))

    assert result.summary["decision"] == "INDETERMINATE"
    assert result.summary["decision_reason"] == "MULTIPLE_TAXA_COMPATIBLE"
    assert set(result.summary["prediction_set"].split(",")) == {"target", "neighbor"}
    assert result.summary["bootstrap_tie_fraction"] == 1.0
    assert result.summary["unresolved_loci"] == 3


def test_conflicting_repeat_and_snp_channels_are_indeterminate():
    result = _assign(_rows(0.1, 1.0, 1.0, 0.1))

    assert result.summary["decision"] == "INDETERMINATE"
    assert result.summary["decision_reason"] == "REPEAT_SNP_EVIDENCE_DISAGREES"


def test_too_few_callable_loci_are_indeterminate():
    rows = _rows(0.1, 0.1, 1.0, 1.0, loci=2)
    result = assign_target_taxon(
        sample_id="sample",
        target_taxon_id="target",
        locus_marker_rows=rows,
        placement_rows=_placements(2),
        reference_metadata=_metadata(),
        calibration=_calibration(),
        min_locus_fraction=1.0,
        bootstrap_replicates=20,
    )

    assert result.summary["decision"] == "INDETERMINATE"
    assert result.summary["decision_reason"] == "INSUFFICIENT_CALLABLE_LOCI"
    assert "INSUFFICIENT_CALLABLE_LOCI" in result.summary["qc_flags"]


def test_placement_qc_can_force_indeterminate():
    result = assign_target_taxon(
        sample_id="sample",
        target_taxon_id="target",
        locus_marker_rows=_rows(0.1, 0.1, 1.0, 1.0),
        placement_rows=_placements(entropy=2.0, lwr=0.1),
        reference_metadata=_metadata(),
        calibration=_calibration(),
        min_locus_fraction=1.0,
        bootstrap_replicates=20,
        max_mean_placement_entropy=1.0,
        min_median_placement_lwr=0.5,
    )

    assert result.summary["decision"] == "INDETERMINATE"
    assert "HIGH_PLACEMENT_ENTROPY" in result.summary["qc_flags"]
    assert "LOW_PLACEMENT_LWR" in result.summary["qc_flags"]


def test_calibration_round_trip(tmp_path):
    calibration = _calibration()
    path = calibration.write(tmp_path / "taxon_calibration.json")

    assert TaxonCalibration.read(path) == calibration


def test_leave_one_out_calibration_excludes_self_reference():
    metadata = _metadata()
    rows = []
    profiles = {
        "T1": (0.0, 0.0),
        "T2": (0.2, 0.2),
        "N1": (2.0, 2.0),
        "N2": (2.2, 2.2),
    }
    for query_id, (query_repeat, query_snp) in profiles.items():
        for locus_number in range(1, 4):
            for reference_id, (reference_repeat, reference_snp) in profiles.items():
                rows.append(
                    {
                        "query_reference_id": query_id,
                        "reference_id": reference_id,
                        "locus_id": f"L{locus_number}",
                        "normalized_repeat_distance": abs(
                            query_repeat - reference_repeat
                        ),
                        "normalized_snp_distance": abs(query_snp - reference_snp),
                    }
                )

    calibration, score_rows = build_taxon_calibration(
        reference_locus_rows=rows,
        reference_metadata=metadata,
        panel_sha256="panel",
        database_signature="database",
        k=1,
        minimum_loci=3,
    )

    assert calibration.taxon_counts == {"neighbor": 2, "target": 2}
    assert len(score_rows) == 12
    assert all(float(row["within_distance"]) > 0 for row in score_rows)


def test_calibration_command_writes_signed_artifacts(tmp_path):
    metadata_path = tmp_path / "reference_metadata.tsv"
    metadata_path.write_text(
        "reference_id\ttaxon_id\ttaxon_name\n"
        "T1\ttarget\tTarget taxon\n"
        "T2\ttarget\tTarget taxon\n"
        "N1\tneighbor\tNeighbor taxon\n"
        "N2\tneighbor\tNeighbor taxon\n"
    )
    index_path = tmp_path / "reference_sequence_index.tsv"
    index_path.write_text(
        "panel_sha256\tdatabase_signature\n"
        "panel\tdatabase\n"
        "panel\tdatabase\n"
    )
    profiles = {
        "T1": (0.0, 0.0),
        "T2": (0.2, 0.2),
        "N1": (2.0, 2.0),
        "N2": (2.2, 2.2),
    }
    distances_path = tmp_path / "reference_distances.tsv"
    lines = [
        "query_reference_id\treference_id\tlocus_id\t"
        "normalized_repeat_distance\tnormalized_snp_distance"
    ]
    for query_id, (query_repeat, query_snp) in profiles.items():
        for locus_number in range(1, 4):
            for reference_id, (reference_repeat, reference_snp) in profiles.items():
                lines.append(
                    f"{query_id}\t{reference_id}\tL{locus_number}\t"
                    f"{abs(query_repeat - reference_repeat)}\t"
                    f"{abs(query_snp - reference_snp)}"
                )
    distances_path.write_text("\n".join(lines) + "\n")

    result = run_taxon_calibration(
        reference_distances_path=distances_path,
        reference_metadata_path=metadata_path,
        sequence_index_path=index_path,
        outdir=tmp_path / "calibration",
        k=1,
        minimum_loci=3,
    )

    calibration = TaxonCalibration.read(result["calibration"])
    assert calibration.panel_sha256 == "panel"
    assert calibration.database_signature == "database"
    assert result["scores"].read_text().count("\n") == 13


def test_phylogeny_integration_writes_assignment_tables_and_checks_signatures(
    tmp_path,
):
    phylogeny = tmp_path / "phylogeny"
    phylogeny.mkdir()
    marker_path = phylogeny / "locus_marker_distances.tsv"
    marker_rows = _rows(0.1, 0.1, 1.0, 1.0)
    with marker_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "locus_id",
                "reference_id",
                "normalized_repeat_distance",
                "normalized_snp_distance",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(marker_rows)
    placement_path = phylogeny / "locus_phylogenetic_distances.tsv"
    placement_rows = _placements()
    with placement_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "locus_id",
                "reference_id",
                "placement_entropy",
                "like_weight_ratio",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(placement_rows)
    calibration_path = _calibration().write(tmp_path / "calibration.json")
    paths = {
        "phylogeny": phylogeny,
        "locus_marker_distances": marker_path,
        "phylogenetic_distances": placement_path,
    }

    result = _add_taxon_assignment_outputs(
        paths,
        sample_id="sample",
        target_taxon_id="target",
        calibration_path=calibration_path,
        reference_metadata=_metadata(),
        panel_sha256="panel",
        database_signature="database",
        expected_loci=3,
        alpha=None,
        min_loci=None,
        min_locus_fraction=1.0,
        bootstrap_replicates=20,
        min_bootstrap_support=0.95,
        max_mean_placement_entropy=None,
        min_median_placement_lwr=None,
    )

    assert result["taxon_assignment"].is_file()
    assert result["taxon_assignment_candidates"].is_file()
    assert result["taxon_assignment_loci"].is_file()
    assert "\tPOSITIVE\t" in result["taxon_assignment"].read_text()

    with pytest.raises(ValueError, match="panel signature"):
        _add_taxon_assignment_outputs(
            paths,
            sample_id="sample",
            target_taxon_id="target",
            calibration_path=calibration_path,
            reference_metadata=_metadata(),
            panel_sha256="wrong-panel",
            database_signature="database",
            expected_loci=3,
            alpha=None,
            min_loci=None,
            min_locus_fraction=1.0,
            bootstrap_replicates=20,
            min_bootstrap_support=0.95,
            max_mean_placement_entropy=None,
            min_median_placement_lwr=None,
        )


def test_target_must_be_present_in_labeled_metadata():
    with pytest.raises(ValueError, match="absent from reference metadata"):
        assign_target_taxon(
            sample_id="sample",
            target_taxon_id="missing",
            locus_marker_rows=_rows(0.1, 0.1, 1.0, 1.0),
            placement_rows=_placements(),
            reference_metadata=_metadata(),
            calibration=_calibration(),
        )