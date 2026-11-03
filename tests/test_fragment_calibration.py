"""Exact overlaps describe a selected fragment component, never an unseen tail."""

from dataclasses import replace
import random
import shutil

import pytest

from mlvamaps.fragment_calibration import OverlapFragmentCollector
from mlvamaps.models import ReadPair, ReadRecord
from mlvamaps.sequence import revcomp


def pairs_for_lengths(lengths):
    rng = random.Random(195)
    result = []
    for index, length in enumerate(lengths):
        sequence = "".join(rng.choices("ACGT", k=length))
        name = f"molecule{index}"
        result.append(ReadPair(name, ReadRecord(name + "/1", sequence[:100]),
                               ReadRecord(name + "/2", revcomp(sequence[-100:]))))
    return result


def calibrate(pairs, max_pairs=5000):
    collector = OverlapFragmentCollector(max_pairs)
    for pair in pairs:
        collector.add(pair)
    return collector.finish()


def test_calibration_is_bounded_deterministic_and_deduplicates_molecules():
    pairs = pairs_for_lengths([179, 180, 181] * 100)
    expected, diagnostic = calibrate(pairs, max_pairs=151)
    assert expected is not None and expected["sd_bp"] == 1
    assert expected["retained_fragments"] == diagnostic["sampled_fragments"] == 151
    assert diagnostic["sample_limited"]
    shuffled = pairs * 4
    random.Random(47).shuffle(shuffled)
    repeated, repeated_diagnostic = calibrate(shuffled, max_pairs=151)
    assert repeated == expected
    assert repeated_diagnostic["sampled_fragments"] == 151
    assert not repeated["whole_library"] and not repeated["unobserved_tail_excluded"]
    assert calibrate(pairs[:100] * 4)[0] is None


@pytest.mark.parametrize("lengths,reason", [
    ([180] * 150 + [160] * 50, "heterogeneous_overlap_lengths"),
    (list(range(110, 187)) * 3, "broad_overlap_lengths"),
    ([173] * 60 + [179] * 60, "multimodal_overlap_lengths"),
    ([185] * 120, "overlap_edge_censoring"),
    ([102] * 120, "overlap_edge_censoring"),
    (list(range(174, 221)) * 8, "overlap_edge_censoring"),
])
def test_observed_mixtures_broad_and_censored_samples_are_rejected(lengths, reason):
    component, diagnostics = calibrate(pairs_for_lengths(lengths))
    assert component is None
    assert diagnostics["calibration_status"] == reason


def test_unseen_long_component_remains_explicitly_unresolved():
    component, diagnostic = calibrate(pairs_for_lengths([180] * 150 + [260] * 150))
    assert component is not None and component["mean_bp"] == 180
    assert diagnostic["eligible_fragments"] == 150
    assert not component["whole_library"] and not component["unobserved_tail_excluded"]


def test_short_low_complexity_periodic_and_conflicting_overlaps_are_not_calibrators():
    short = pairs_for_lengths([190] * 120)
    assert calibrate(short)[0] is None  # Only 10 bp overlap.
    rng = random.Random(76)
    for overlap in ("A" * 18 + "CG", "ATGC" * 5, "NACG" * 5):
        pairs = []
        for index in range(120):
            left = "".join(rng.choices("ACGT", k=80)) + overlap
            right = overlap + "".join(rng.choices("ACGT", k=80))
            pairs.append(ReadPair(str(index), ReadRecord(str(index), left),
                                  ReadRecord(str(index), revcomp(right))))
        assert calibrate(pairs)[0] is None
    pairs = pairs_for_lengths([180] * 101)
    bad = replace(pairs[0], read1=replace(pairs[0].read1, sequence="A" * 100))
    assert calibrate(pairs + [bad]) == calibrate([bad] + pairs)
    assert calibrate(pairs + [bad])[1]["conflicting_molecule_ids"] == 1


@pytest.mark.skipif(not shutil.which("minimap2"), reason="Challenge generator builds a native index")
def test_unchanged_challenge_calibrates_from_fastq_without_truth(tmp_path):
    from mlvamaps.io import read_fastq_pairs
    from scripts.fastq_challenge_fixture import prepare_challenge

    _, inputs, _ = prepare_challenge(tmp_path / "challenge", "paired", copies=2)
    pairs = list(read_fastq_pairs(*inputs))
    component, diagnostic = calibrate(pairs)
    assert diagnostic["eligible_fragments"] > 200
    assert component["retained_fragments"] > 200
    assert component["mean_bp"] == pytest.approx(180, abs=.1)
    assert component["sd_bp"] == 1
    # Names, input ordering, and generator truth never enter calibration.
    renamed = [replace(pair, molecule_id=f"opaque{index}") for index, pair in enumerate(pairs)]
    random.Random(58).shuffle(renamed)
    assert calibrate(renamed)[0] == component
