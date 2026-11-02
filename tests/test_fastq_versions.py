"""Comparison must preserve status changes, missing calls, and both directions."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from benchmark_fastq_inference import challenge_metrics, file_sha256, version_disagreements


def test_version_audit_retains_uncertainty_and_does_not_pick_a_winner(tmp_path):
    before = {("s", "a"): {"repeat_count": 5, "status": "PASS"},
              ("s", "b"): {"repeat_count": None, "status": "UNRESOLVED"}}
    after = {("s", "a"): {"repeat_count": 5, "status": "ESTIMATED"},
             ("s", "c"): {"repeat_count": 7, "status": "PASS"}}
    result = version_disagreements(before, after)
    assert [row["locus"] for row in result] == ["a", "b", "c"]
    assert result[0]["baseline"] == before[("s", "a")]
    assert result[0]["updated"] == after[("s", "a")]
    assert version_disagreements(before, before) == []
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"abc")
    assert file_sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_challenge_metrics_keep_missing_calls_in_denominator(tmp_path):
    import json
    output = tmp_path / "updated_0" / "calls"
    output.mkdir(parents=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "read_truth.json").write_text(json.dumps({"a": "correct", "b": "wrong", "c": "missing"}))
    (output / "calls.tsv").write_text(
        "sample_id\tlocus_id\trepeat_count\tallele_confidence\n"
        "s\tcorrect\t5\t.99\ns\twrong\t4\t.6\ns\tmissing\t\t0\ns\tabsent\t\t0\n")
    (output / "molecule_candidate_evidence.tsv").write_text(
        "molecule\tlocus\tquality_flags\na\tcorrect\t\na\tcorrect\t\n"
        "b\tcorrect\t\nc\twrong\tambiguous_locus\n")
    truth = {name: dict(repeat_count=value, alleles=[value], coverage=3,
                       motif_length=4, repeat_bases=20, error_rate=0)
             for name, value in (("correct", 5), ("wrong", 9), ("missing", 7), ("absent", None))}
    result = challenge_metrics(output, truth)
    assert result["call_rate"] == 2 / 3
    assert result["accuracy_among_calls"] == .5
    assert result["mae"] == 2.5
    assert result["false_positive_genotype_rate"] == 0
    assert result["locus_assignment_accuracy"] == .5
    assert result["locus_assignment_rate"] == 2 / 3
    assert result["accuracy_call_rate_curve"][-1]["accuracy"] == 1
    assert result["accuracy_call_rate_curve"][-1]["call_rate"] == 1 / 3
