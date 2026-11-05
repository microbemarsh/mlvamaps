import csv
from dataclasses import replace

import pytest

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles
from mlvamaps.evidence_reporting import write_evidence_reports
from mlvamaps.io import write_tsv
from mlvamaps.models import Locus


def test_exact_requires_independent_complete_unambiguous_evidence(tmp_path):
    locus = Locus("L", repeat_motif="ATGC", repeat_unit_length_bp=4)
    row = CandidateEvidence("L", 5, "m", measured_repeat_count=5, raw_repeat_count=5,
                            repeat_sigma=.08, direct_product_measurement=True, full_repeat_span=True)
    limits = InferenceThresholds(minimum_molecules=1)
    for count, expected in ((1, "ESTIMATED"), (3, "EXACT")):
        evidence = [replace(row, molecule_id=str(i)) for i in range(count)] * 5
        calls, _ = infer_alleles(evidence, [locus], [], "s", "hifi", limits, require_measurement=True)
        paths = write_evidence_reports(calls, evidence, [locus], tmp_path, limits)
        assert calls[0]["genotype_status"] == expected
        with paths["genotype_evidence"].open() as handle:
            assert int(next(csv.DictReader(handle, delimiter="\t"))["independent_support"]) == count
    evidence = [replace(row, molecule_id=str(i), quality_flags=("boundary_ambiguous",)) for i in range(3)]
    calls, _ = infer_alleles(evidence, [locus], [], "s", "ont", limits, require_measurement=True)
    write_evidence_reports(calls, evidence, [locus], tmp_path, limits)
    assert calls[0]["genotype_status"] == "ESTIMATED"


@pytest.mark.parametrize("details, expected", [
    ({}, "UNRESOLVED"),
    ({"prior_only": True}, "PRIOR_ONLY"),
    ({"repeat_lower_bound": 3}, "LOWER_BOUND"),
    ({"repeat_lower_bound": 3, "finite_length_support": True}, "UNRESOLVED"),
    ({"possible_secondary": True}, "POSSIBLE_MIXTURE"),
])
def test_withheld_counts_preserve_evidence_category_and_rejection_reason(tmp_path, details, expected):
    locus = Locus("L")
    limits = InferenceThresholds()
    calls, _ = infer_alleles([], [locus], [], "s", "hifi", limits, require_measurement=True)
    calls[0].update(repeat_count="", repeat_count_raw="", predicted_repeat_count=5,
        predicted_repeat_count_raw=5.2, call_accepted=False,
        call_rejection_reason="Insufficient independent support.", estimation_reason="Uncertain evidence.", **details)
    paths = write_evidence_reports(calls, [], [locus], tmp_path, limits)
    with paths["genotype_evidence"].open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["genotype_status"] == expected and row["repeat_count"] == ""
    assert "Call withheld: Insufficient independent support." in row["reason"]
    assert row["call_accepted"] == "False" and row["predicted_repeat_count"] == "5"
    assert row["predicted_repeat_count_raw"] == "5.2" and row["repeat_count_raw"] == ""


def test_equal_repeat_counts_keep_distinct_observed_sequence_alleles(tmp_path):
    locus = Locus("L", left_flank_sequence="GACTGACCTGATCGACGTTAGCTA",
                  right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC", repeat_motif="ATGC", repeat_unit_length_bp=4)
    products = [locus.left_flank_sequence + repeat + locus.right_flank_sequence
                for repeat in ("ATGC" * 5, "ATGC" * 2 + "ACGC" + "ATGC" * 2)]
    evidence = [CandidateEvidence("L", 5, f"{j}_{i}", measured_repeat_count=5, raw_repeat_count=5,
        repeat_sigma=.08, direct_product_measurement=True, metadata={"product_sequence": product})
        for j, product in enumerate(products) for i in range(3)]
    calls, _ = infer_alleles(evidence, [locus], [], "s", "hifi", require_measurement=True)
    paths = write_evidence_reports(calls, evidence, [locus], tmp_path, InferenceThresholds())
    with paths["sequence_alleles"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2 and len({row["sequence_id"] for row in rows}) == 2
    assert {row["repeat_count"] for row in rows} == {"5"}
    assert {row["motif_mismatches"] for row in rows} == {"0", "1"}


def test_ambiguous_bases_cannot_supply_exact_or_sequence_alleles(tmp_path):
    locus = Locus("L", repeat_motif="ATGC", repeat_unit_length_bp=4)
    evidence = [CandidateEvidence("L", 5, str(i), measured_repeat_count=5,
        raw_repeat_count=5, repeat_sigma=.08, direct_product_measurement=True,
        metadata={"product_sequence": "ATGC" * 4 + "ATGN"}) for i in range(3)]
    calls, _ = infer_alleles(evidence, [locus], [], "s", "hifi", require_measurement=True)
    paths = write_evidence_reports(calls, evidence, [locus], tmp_path, InferenceThresholds())
    assert calls[0]["genotype_status"] == "ESTIMATED"
    with paths["sequence_alleles"].open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []


def test_atomic_tsv_failure_preserves_previous_result_and_removes_temporary(tmp_path):
    path = tmp_path / "calls.tsv"
    path.write_text("previous complete result\n")
    def broken_rows():
        yield {"locus": "L"}
        raise RuntimeError("interrupted evidence generator")
    with pytest.raises(RuntimeError, match="interrupted"):
        write_tsv(broken_rows(), path, ["locus"])
    assert path.read_text() == "previous complete result\n"
    assert list(tmp_path.iterdir()) == [path]


def test_cli_exposes_positive_resource_limits():
    from mlvamaps.cli import build_parser, _fastq_limits
    parser = build_parser()
    args = parser.parse_args(["call", "--fastq-max-candidates", "32"])
    assert _fastq_limits(args).maximum_candidates == 32
    with pytest.raises(SystemExit):
        parser.parse_args(["call", "--fastq-max-likelihood-cells", "0"])
