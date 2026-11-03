"""Read pairs can identify a conditional mode without proving a finite array."""
from dataclasses import replace
import math
from pathlib import Path
import random
import runpy

from mlvamaps.allele_inference import InferenceThresholds
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.fragment_calibration import OverlapFragmentCollector
from mlvamaps.io import read_fastq_pairs
from mlvamaps.mapping_classification import observed_reference_contexts
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.primers import read_loci_or_primers
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_pair_rescue import (
    WeakPairCollector, _model_score, _substitution_score, infer_pair_repeat_rescue,
)
from mlvamaps.short_read_likelihood import _log_normal_interval


def example():
    rng = random.Random(731)
    left, right = ("".join(rng.choices("ACGT", k=70)) for _ in range(2))
    left, right = left[:-1] + "A", "C" + right[1:]
    locus = Locus("L", forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
        left_flank_sequence=left[20:], right_flank_sequence=right[:-20], repeat_motif="ATGC",
        repeat_unit_length_bp=4, nominal_repeat_units=4, expected_product_size_bp=156)
    context = CandidateContext("ref", "L", left + "ATGC" * 4 + right, 4, 4, 70, 86,
                               left, right, repeat_motif="ATGC", observed_or_synthetic="observed")
    product = left + "ATGC" * 40 + right
    pairs = []
    for copy in range(2):
        for start in (15, 33, 51, 69, 87, 105):
            fragment = product[start:start + 180]
            identifier = f"m{copy}_{start}"
            pairs.append(ReadPair(identifier, ReadRecord(identifier + "/1", fragment[:100], "I" * 100),
                                  ReadRecord(identifier + "/2", revcomp(fragment[-100:]), "I" * 100)))
    return locus, context, pairs


def infer(loci, contexts, pairs, *, limits=None, maximum=2000, **kwargs):
    collector = WeakPairCollector(contexts, maximum_per_locus=maximum)
    for pair in pairs:
        collector.add(pair)
    calibration = dict(mean_bp=180., sd_bp=1., source="exact_overlap_component", whole_library=False)
    updates, diagnostics = infer_pair_repeat_rescue(loci, contexts, collector, calibration,
                                                   limits or InferenceThresholds(), **kwargs)
    return updates, diagnostics, collector


def test_original_challenge_uses_all_pairs_without_truth_or_assembly(tmp_path, monkeypatch):
    from mlvamaps import minimap_mapping
    monkeypatch.setattr(minimap_mapping, "build_competitive_indexes", lambda *args: None)
    fixture = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/fastq_challenge_fixture.py"))
    fixture["prepare_challenge"](tmp_path / "inputs", "paired", copies=2)
    source = tmp_path / "inputs"
    loci = read_loci_or_primers(str(source / "panel.tsv"), None)
    _, contexts, _ = observed_reference_contexts(source / "database", loci, .25)
    collector, calibrator = WeakPairCollector(contexts), OverlapFragmentCollector()
    for pair in read_fastq_pairs(source / "r1.fastq", source / "r2.fastq"):
        collector.add(pair)
        calibrator.add(pair)
    calibration, _ = calibrator.finish()
    updates, _ = infer_pair_repeat_rescue(loci, contexts, collector, calibration, InferenceThresholds())
    assert set(updates) == {"expansion"}
    call = updates["expansion"]
    assert call["repeat_count"] == 40 and call["pair_rescue_support"] == 12
    assert call["pair_rescue_boundary_support"] == 4
    assert call["status"] == "estimated" and call["upper_unbounded"]
    assert not call["finite_length_support"] and call["repeat_interval_high"] == ""
    assert "pure periodic array" in call["estimation_reason"]
    assert "measured_repeat_count" not in call


def test_pair_ids_order_mate_order_and_duplicate_references_do_not_change_fit():
    locus, context, pairs = example()
    first, _, _ = infer([locus], [context], pairs)
    changed = []
    for index, pair in enumerate(reversed(pairs)):
        identifier = f"opaque{index}"
        changed.append(ReadPair(identifier, replace(pair.read2, read_id=identifier + "/1"),
                                replace(pair.read1, read_id=identifier + "/2")))
    second, _, _ = infer([locus], [context, replace(context, candidate_id="copy")] * 3, changed + changed)
    assert first["L"]["repeat_count"] == second["L"]["repeat_count"] == 40
    for key in ("best_probability", "candidate_distribution", "pair_rescue_support", "upper_unbounded"):
        assert first["L"][key] == second["L"][key]


def test_competing_locus_remains_a_competitor_when_it_is_not_a_rescue_target():
    locus, context, pairs = example()
    other = replace(locus, locus_id="other")
    competing = replace(context, candidate_id="other", locus_id="other")
    updates, diagnostics, _ = infer([locus, other], [context, competing], pairs, target_loci={"L"})
    assert updates == {}
    assert next(row for row in diagnostics if row["locus"] == "L")["competing_pairs"] == 12


def test_pure_mates_and_one_sided_junctions_cannot_make_a_rescue():
    locus, context, pairs = example()
    for selected in ([pair for pair in pairs if pair.molecule_id.split("_")[-1] in {"15", "33", "87", "105"}],
                     [pair for pair in pairs if pair.molecule_id.split("_")[-1] != "69"]):
        updates, diagnostics, _ = infer([locus], [context], selected)
        assert updates == {}
        assert diagnostics[0]["status"] == "insufficient_reciprocal_junctions"


def test_budget_and_pool_overflow_fail_closed():
    locus, context, pairs = example()
    updates, diagnostics, collector = infer([locus], [context], pairs, maximum=3)
    assert updates == {} and diagnostics[0]["status"] == "candidate_pool_overflow"
    assert collector.diagnostics["overflow_loci"] == ["L"] and not collector.fragments["L"]
    updates, diagnostics, _ = infer([locus], [context], pairs,
                                   limits=InferenceThresholds(maximum_likelihood_cells=10))
    assert updates == {} and diagnostics[0]["status"] == "pair_profile_budget_exceeded"


def test_single_base_can_be_an_error_and_placements_are_marginalized():
    error = _substitution_score("ACGT", "ACGA", "IIII")
    assert math.isfinite(error) and error < 0
    assert _substitution_score("ACGT", "ACGA", "5555") > error
    assert _substitution_score("ACGT", "AATA", "IIII") == -math.inf
    weights = {length: _log_normal_interval((length - .5 - 180), (length + .5 - 180))
               for length in range(172, 189)}
    model = ("left", 19, 81, "ATGC" * 25, "I" * 100, "C" * 70, "ATGC", -1, 4)
    score, boundary = _model_score(model, 200, 180, 1, weights)
    matching = [weight for length, weight in weights.items() if (length - 119) % 4 == 0]
    assert math.exp(score) > math.exp(max(matching))
    assert not boundary


def test_interrupted_longer_array_counterexample_never_becomes_exact():
    locus, context, pairs = example()
    left, right = context.sequence[:70], context.sequence[86:]
    repeat = list("ATGC" * 100)
    repeat[160], repeat[239] = right[0], left[-1]
    long_product = left + "".join(repeat) + right
    for pair in pairs:
        assert any(long_product[start:start + 100] == pair.read1.sequence
                   and long_product[start + 80:start + 180] == revcomp(pair.read2.sequence)
                   for start in range(len(long_product) - 179))
    updates, _, _ = infer([locus], [context], pairs)
    assert updates["L"]["repeat_count"] == 40
    assert updates["L"]["upper_unbounded"] and not updates["L"]["finite_length_support"]
    assert updates["L"]["status"] == "estimated"


def test_pure_mate_crossing_motif_like_flank_is_not_a_junction():
    weights = {180: 0.}
    model = ("left", 19, 81, "TGCA" * 25, "I" * 100, "ATGC" * 20, "ATGC", -1, 4)
    score, boundary = _model_score(model, 160, 180, 1, weights)
    assert math.isfinite(score) and not boundary


def test_broad_external_fragment_model_uses_cached_sequence_profiles(monkeypatch):
    from mlvamaps import short_read_pair_rescue as module
    locus, context, pairs = example()
    collector = WeakPairCollector([context])
    for pair in pairs:
        collector.add(pair)
    comparisons = 0
    original = module._substitution_score

    def counted(*args):
        nonlocal comparisons
        comparisons += 1
        return original(*args)

    monkeypatch.setattr(module, "_substitution_score", counted)
    _, diagnostics = infer_pair_repeat_rescue([locus], [context], collector,
        dict(mean_bp=350., sd_bp=50., source="user", whole_library=True), InferenceThresholds())
    assert diagnostics[0]["status"] != "pair_likelihood_budget_exceeded"
    assert comparisons < 5000  # Fixed mate profiles, not 751 states x 800 inserts.


def test_conflicting_molecule_ids_are_excluded_in_either_order():
    locus, context, pairs = example()
    changed = pairs[2].read2.sequence[:-1] + ("A" if pairs[2].read2.sequence[-1] != "A" else "C")
    collision = replace(pairs[2], read2=replace(pairs[2].read2, sequence=changed))
    for ordered in (pairs + [collision], [collision] + pairs):
        updates, _, collector = infer([locus], [context], ordered)
        assert collector.diagnostics["conflicting_molecules"] == 1
        assert not updates or pairs[2].molecule_id not in updates["L"]["_molecule_ids"]


def test_nonjunction_seed_decoys_do_not_fill_the_pool():
    locus, context, pairs = example()
    collector = WeakPairCollector([context], maximum_per_locus=12)
    for index in range(100):
        read = ReadRecord(str(index), context.sequence[:12] + "A" * 88, "I" * 100)
        collector.add(ReadPair(str(index), read, read))
    for pair in pairs:
        collector.add(pair)
    assert not collector.overflow_loci
    assert len(collector.fragments["L"]) == 12


def test_missing_calibration_and_mixed_sample_do_not_get_pair_rescues():
    locus, context, pairs = example()
    collector = WeakPairCollector([context])
    for pair in pairs:
        collector.add(pair)
    updates, diagnostics = infer_pair_repeat_rescue([locus], [context], collector, None, InferenceThresholds())
    assert updates == {} and diagnostics[0]["status"] == "no_fragment_calibration"
    updates, diagnostics, _ = infer([locus], [context], pairs, sample_mode="metagenome")
    assert updates == {} and diagnostics[0]["status"] == "unsupported_sample_mode"


def test_invalid_derived_calibration_and_profile_budget_guard(monkeypatch):
    from mlvamaps import short_read_pair_rescue as module
    locus, context, pairs = example()
    collector = WeakPairCollector([context])
    for pair in pairs:
        collector.add(pair)
    updates, diagnostics = infer_pair_repeat_rescue([locus], [context], collector,
        dict(mean_bp=1e308, sd_bp=1e308), InferenceThresholds())
    assert updates == {} and diagnostics[0]["status"] == "invalid_fragment_calibration"

    def unexpected_profile(*args):
        raise AssertionError("The work limit must be checked before computing profiles")

    monkeypatch.setattr(module, "_profile", unexpected_profile)
    updates, diagnostics = infer_pair_repeat_rescue([locus], [context], collector,
        dict(mean_bp=180, sd_bp=1), InferenceThresholds(maximum_likelihood_cells=10))
    assert updates == {} and diagnostics[0]["status"] == "pair_profile_budget_exceeded"
