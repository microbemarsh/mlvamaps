"""Native inference of an expansion from a calibrated, selected fragment pool."""

import csv
from dataclasses import replace
import shutil
from pathlib import Path
import runpy
from statistics import NormalDist

import pytest

from mlvamaps.candidate_contexts import generate_candidate_contexts, write_candidate_contexts
from mlvamaps.cli import main
from mlvamaps.io import write_fastq, write_tsv
from mlvamaps.minimap_mapping import build_competitive_indexes
from mlvamaps.models import ReadRecord
from mlvamaps.sequence import revcomp
from test_fragment_length_evidence import partial_pair


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"),
                   reason="Native tools are needed for the full length-inference workflow")
def test_calibrated_partial_pairs_recover_novel_length_without_measured_votes(tmp_path, monkeypatch):
    from mlvamaps import mapping_classification

    locus, context, _, product = partial_pair()
    locus = replace(locus, expected_min_repeats=0, expected_max_repeats=12)
    panel = tmp_path / "panel.tsv"
    write_tsv([vars(locus)], panel, list(vars(locus)))
    database = tmp_path / "database"
    database.mkdir()
    (database / "L.fasta").write_text(f">reference\n{context.sequence}\n")
    assets = write_candidate_contexts(generate_candidate_contexts([locus], database),
                                     database / "competitive_mapping")
    build_competitive_indexes(assets["fasta"], database / "competitive_mapping")
    (database / "manifest.json").write_text("{}")

    first, second = [], []
    distribution = NormalDist(260, 30)
    for index in range(2048):
        length = round(distribution.inv_cdf((index + .5) / 2048))
        right_anchor = length - 160 - 30
        # Retain exactly the modeled union: repeat-only opposite mates or
        # eligible right-boundary mates. All fragments have a 30 bp left anchor.
        if not (130 <= length <= 190 or 12 <= right_anchor <= 70):
            continue
        fragment = product[40:40 + length]
        identifier = f"fragment{index}"
        first.append(ReadRecord(identifier + "/1", fragment[:100], "I" * 100))
        second.append(ReadRecord(identifier + "/2", revcomp(fragment[-100:]), "I" * 100))
    first_path, second_path = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    write_fastq(first, first_path)
    write_fastq(second, second_path)

    measurements = []
    classify = mapping_classification.run_mapping_classification

    def record(**kwargs):
        measurements.append(dict(kwargs.get("measured_repeats") or {}))
        return classify(**kwargs)

    monkeypatch.setattr(mapping_classification, "run_mapping_classification", record)
    calls = {}
    for calibrated in (False, True):
        output = tmp_path / ("calibrated" if calibrated else "uncalibrated")
        argv = ["call", "-p", str(panel), "--database", str(database), "-o", str(output),
                "--sample-id", "sample", "--sample-mode", "isolate", "-t", "2", "--quiet",
                "-i", "sr", "--fq1", str(first_path), "--fq2", str(second_path)]
        if calibrated:
            argv.extend(["--short-insert-mean", "260", "--short-insert-sd", "30"])
        assert main(argv) == 0
        with (output / "calls.tsv").open() as handle:
            calls[calibrated], = csv.DictReader(handle, delimiter="\t")
    assert calls[True]["predicted_repeat_count"] == "40"
    assert float(calls[False]["predicted_repeat_count"]) < 40
    assert all(call["repeat_count"] == "" and call["call_accepted"] == "False" for call in calls.values())
    assert calls[True]["inference_method"] == "joint_short_read_length"
    assert calls[True]["genotype_status"] == "UNRESOLVED"
    assert calls[True]["finite_length_support"] == "True"
    assert all(not measured for measured in measurements)


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"),
                   reason="Native tools are needed for the original read-mapping challenge")
@pytest.mark.parametrize("existing_lower_bound", [0, 45])
def test_original_expansion_is_rescued_without_assembly_or_measured_votes(tmp_path, monkeypatch, existing_lower_bound):
    from mlvamaps import mapping_classification, short_read_fragments, local_assembly, unified_fastq
    from mlvamaps.io import read_fasta

    fixture = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/fastq_challenge_fixture.py"))
    argv, _, truth = fixture["prepare_challenge"](tmp_path / "inputs", "paired", copies=2)
    argv.append("--short-no-assembly")
    output = Path(argv[argv.index("-o") + 1])
    measurements = []
    classify = mapping_classification.run_mapping_classification

    def record(**kwargs):
        measurements.append(dict(kwargs.get("measured_repeats") or {}))
        return classify(**kwargs)

    def unexpected_assembly(*args, **kwargs):
        raise AssertionError("The mapping-only workflow must not invoke assembly or consensus")

    monkeypatch.setattr(mapping_classification, "run_mapping_classification", record)
    monkeypatch.setattr(short_read_fragments, "assemble_fragments", unexpected_assembly)
    monkeypatch.setattr(local_assembly, "sequence_consensus", unexpected_assembly)
    if existing_lower_bound:
        estimate = unified_fastq.estimate_unresolved_calls

        def with_stronger_bound(calls, *args, **kwargs):
            estimate(calls, *args, **kwargs)
            # A bound from another evidence path must survive a conflicting
            # conditional pair model; it cannot be silently overwritten.
            next(call for call in calls if call["locus"] == "expansion").update(
                repeat_count=45, repeat_count_raw=45, repeat_lower_bound=45,
                repeat_interval_low=45)

        monkeypatch.setattr(unified_fastq, "estimate_unresolved_calls", with_stronger_bound)
    assert main(argv) == 0
    with (output / "calls.tsv").open() as handle:
        calls = {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    call = calls["expansion"]
    assert float(call["predicted_repeat_count"]) == (existing_lower_bound or truth["expansion"]["repeat_count"])
    assert call["repeat_count"] == "" and call["call_accepted"] == "False"
    assert call["genotype_status"] == "LOWER_BOUND"
    assert call["upper_unbounded"] == "True" and call["repeat_interval_high"] == ""
    assert call["finite_length_support"] == "False"
    assert all(locus != "expansion" for locus, molecule in measurements[0])
    assert "expansion" not in dict(read_fasta(output / "query_loci.fasta"))
    assert sum(float(calls[name]["predicted_repeat_count"]) in values["alleles"]
               for name, values in truth.items() if values["repeat_count"] is not None) >= (9 if existing_lower_bound else 10)
    if existing_lower_bound:
        with (output / "joint_pair_diagnostics.tsv").open() as handle:
            diagnostic = next(row for row in csv.DictReader(handle, delimiter="\t") if row["locus"] == "expansion")
        assert diagnostic["status"] == "conflicts_with_existing_lower_bound"
        assert float(call["repeat_interval_low"]) == existing_lower_bound
    assert calls["zero_coverage"]["prior_only"] == "True"
    assert calls["long_array"]["upper_unbounded"] == "True"
    assert all(calls[name]["genotype_status"] == "AMBIGUOUS_LOCUS"
               for name in ("ambiguous_a", "ambiguous_b"))
