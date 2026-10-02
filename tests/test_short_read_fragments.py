from dataclasses import replace
import csv
import random
import shutil

import pytest

from mlvamaps.allele_inference import infer_alleles
from mlvamaps.candidate_contexts import generate_candidate_contexts
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.io import write_fastq, write_tsv
from mlvamaps.repeat_calibration import assembly_equivalent_product_allele
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_evidence import extract_short_read_evidence
from mlvamaps.short_read_fragments import FlankRecruiter, RecruitedFragment, assemble_fragments


def locus_and_product(count=8, imperfect=False):
    rng = random.Random(14)
    forward, left, right, reverse = ["".join(rng.choices("ACGT", k=length)) for length in (20, 40, 40, 20)]
    motif = "AGCTACGTTCGA"
    locus = Locus("L", forward_primer=forward, reverse_primer=reverse,
                  left_flank_sequence=left, right_flank_sequence=right,
                  repeat_motif=motif, repeat_unit_length_bp=len(motif),
                  nominal_repeat_units=5, expected_min_repeats=3, expected_max_repeats=6)
    repeat = list(motif * count)
    if imperfect:
        for i in range(0, len(repeat), 6):
            repeat[i] = rng.choice("ACGT")
    return locus, forward + left + "".join(repeat) + right + revcomp(reverse)


def pair(name, first, second):
    return ReadPair(name, ReadRecord(name + "/1", first, "I" * len(first)),
                    ReadRecord(name + "/2", second, "I" * len(second)))


def recruited_tiling(locus, product, copies=3):
    recruiter = FlankRecruiter(generate_candidate_contexts([locus]))
    for length in (150, 190, 230):
        for start in sorted(set(range(0, len(product) - length + 1, 5)) | ({len(product) - length} if len(product) >= length else set())):
            for copy in range(copies):
                recruiter.add(pair(f"{length}_{start}_{copy}", product[start:start + 100],
                                   revcomp(product[start + length - 100:start + length])))
    return recruiter.fragments.get(locus.locus_id, {})


def test_flank_recruitment_rescues_mates_and_deduplicates():
    locus, product = locus_and_product()
    recruiter = FlankRecruiter(generate_candidate_contexts([locus]))
    reads = [pair("left", product[:70], locus.repeat_motif * 8),
             pair("right", "N" * 100, revcomp(product[-70:])),
             pair("both", product[:70], revcomp(product[-70:])),
             pair("repeat", locus.repeat_motif * 8, revcomp(locus.repeat_motif * 8))]
    for _ in range(2):
        for fragment in reads:
            recruiter.add(fragment)
    pool = recruiter.fragments["L"]
    assert set(pool) == {"left", "right", "both"}
    assert pool["left"].left and pool["left"].rescued
    assert pool["right"].right and pool["right"].rescued
    assert pool["right"].pair.read1.sequence == "N" * 100
    assert pool["both"].left and pool["both"].right and not pool["both"].rescued


def test_shared_flanks_do_not_recruit_and_discordant_pairs_are_rejected():
    locus, product = locus_and_product()
    contexts = generate_candidate_contexts([locus, replace(locus, locus_id="other")])
    recruiter = FlankRecruiter(contexts)
    recruiter.add(pair("shared", product[:100], revcomp(product[-100:])))
    assert not recruiter.fragments
    contexts = generate_candidate_contexts([locus])
    other = [replace(context, locus_id="other", sequence=context.sequence.translate(str.maketrans("ACGT", "CATG"))) for context in contexts]
    recruiter = FlankRecruiter(contexts + other)
    recruiter.add(pair("discordant", product[:70], revcomp(other[0].sequence[-70:])))
    assert recruiter.ambiguous_fragments == 1
    recruiter.add(pair("repeat", locus.repeat_motif * 8, locus.repeat_motif * 8))
    assert not recruiter.fragments


def test_unmapped_fragments_measure_beyond_candidate_range_without_double_counting():
    locus, product = locus_and_product(12)
    contexts = generate_candidate_contexts([locus], maximum=6)
    recruiter = FlankRecruiter(contexts)
    for i in range(3):
        recruiter.add(pair(str(i), product, revcomp(product)))
    evidence = extract_short_read_evidence([], contexts, [locus], recruited_fragments=recruiter.fragments)
    calls, _ = infer_alleles(evidence + evidence, [locus], contexts, "s", "illumina")
    assert calls[0]["repeat_count"] == 12
    assert calls[0]["molecule_support"] == 3


@pytest.mark.parametrize("count", [5, 8, 12, 16, 20, 24])
def test_imperfect_repeat_assembly_is_not_capped_by_candidates(count):
    locus, product = locus_and_product(count, imperfect=True)
    fragments = recruited_tiling(locus, product)
    sequence, support, reason = assemble_fragments(fragments, generate_candidate_contexts([locus]))
    assert reason == "assembled"
    assert sequence == product
    assert len(support) >= 3
    assert assembly_equivalent_product_allele(locus, len(sequence))[1] == count


def test_graph_rejects_unresolved_repeat_cycles_no_path_and_insufficient_coverage():
    locus, product = locus_and_product(24)
    contexts = generate_candidate_contexts([locus])
    fragments = recruited_tiling(locus, product)
    assert assemble_fragments(fragments, contexts)[2] in {"repeat_cycle", "competing_paths"}
    assert assemble_fragments({}, contexts)[2] == "insufficient_coverage"
    broken = {str(i): RecruitedFragment(pair(str(i), product[:100], "N" * 100), True, False, True)
              for i in range(3)}
    assert not assemble_fragments(broken, contexts)[0]


def test_competing_paths_are_not_collapsed_to_shortest_allele():
    locus, product = locus_and_product(8, imperfect=True)
    altered = product[:100] + ("C" if product[100] != "C" else "A") + product[101:]
    fragments = recruited_tiling(locus, product)
    fragments.update({"alt_" + name: fragment for name, fragment in recruited_tiling(locus, altered).items()})
    assert assemble_fragments(fragments, generate_candidate_contexts([locus]))[2] == "competing_paths"


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native tools unavailable")
@pytest.mark.parametrize("count", [5, 8, 12, 16, 20, 24])
def test_paired_pipeline_recovers_monotonic_counts_through_shared_parser(tmp_path, count):
    from mlvamaps.short_reads import run_short_read_call
    locus, product = locus_and_product(count, imperfect=True)
    fragments = recruited_tiling(locus, product)
    first, second = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    write_fastq((item.pair.read1 for item in fragments.values()), first)
    write_fastq((item.pair.read2 for item in fragments.values()), second)
    panel = tmp_path / "panel.tsv"
    write_tsv([vars(locus)], panel, list(vars(locus)))
    paths = run_short_read_call(str(first), str(second), str(panel), str(tmp_path / "out"), "s",
                                threads=2, short_max_candidate_repeat_count=6, show_progress=False)
    with paths["calls"].open() as handle:
        call = next(csv.DictReader(handle, delimiter="\t"))
    assert float(call["repeat_count"]) == assembly_equivalent_product_allele(locus, len(product))[1] == count
    with paths["short_read_diagnostics"].open() as handle:
        diagnostic = next(csv.DictReader(handle, delimiter="\t"))
    assert diagnostic["microassembly_attempted"] == ("no" if count == 5 else "yes")


def test_repeat_only_competitors_cannot_override_an_anchored_fragment():
    from mlvamaps.alignment_evidence import CandidateAlignment
    locus, product = locus_and_product()
    contexts = generate_candidate_contexts([locus])
    context = contexts[0]
    row = CandidateAlignment("repeat", "repeat/1", 1, "L", context.candidate_id,
                             context.repeat_count, "", 100, 60, 1, 1, 1,
                             0, 30, context.repeat_start, context.repeat_end,
                             "30M", "", True, False, False,
                             query_sequence=locus.repeat_motif * 8)
    assert extract_short_read_evidence([row], contexts, [locus], recruited_fragments={}) == []


def test_unique_merging_rejects_periodic_overlaps_with_a_sequencing_error():
    from mlvamaps.short_reads import merge_read_pair
    repeat = "ACGT" * 30
    first = "G" * 20 + repeat[:80]
    # A single error does not supply positional information within a repeat.
    first = first[:45] + "T" + first[46:]
    second = revcomp(repeat[:80] + "C" * 20)
    assert merge_read_pair(pair("periodic", first, second), require_unique=True) is None


def test_unique_merging_keeps_an_overlap_with_a_mismatch():
    from mlvamaps.short_reads import merge_read_pair
    rng = random.Random(217)
    product = "".join(rng.choices("ACGT", k=150))
    first = product[:100]
    second = list(revcomp(product[-100:]))
    second[70] = "A" if second[70] != "A" else "C"
    merged = merge_read_pair(pair("overlap", first, "".join(second)), require_unique=True)
    assert merged is not None and len(merged.sequence) == 150


def test_graph_coverage_counts_fragments_not_repeated_mates():
    locus, product = locus_and_product(8, imperfect=True)
    fragments = recruited_tiling(locus, product, copies=1)
    isolated = {"one": RecruitedFragment(pair("one", product[:100], revcomp(product[:100])), True, False, False)}
    # Unrelated fragments meet the pool-size threshold but cannot add coverage
    # to the only read covering the left anchor.
    isolated.update({name: fragment for name, fragment in fragments.items() if product[:21] not in fragment.pair.read1.sequence})
    assert not assemble_fragments(isolated, generate_candidate_contexts([locus]))[0]


def test_inconclusive_graph_reports_presence_instead_of_candidate_edge():
    from mlvamaps.alignment_evidence import CandidateEvidence
    from mlvamaps.allele_inference import InferenceThresholds
    from mlvamaps.short_read_fragments import resolve_fragment_calls
    locus, product = locus_and_product(24)
    contexts = generate_candidate_contexts([locus])
    pool = recruited_tiling(locus, product)
    evidence = [CandidateEvidence("L", 6, name, full_repeat_span=True, technology="illumina") for name in pool]
    calls, molecules = infer_alleles(evidence, [locus], contexts, "s", "illumina")
    resolved, _, _, diagnostics, _ = resolve_fragment_calls(
        calls, evidence, molecules, [], contexts, [locus], {"L": pool}, "s", InferenceThresholds(), .25)
    assert resolved[0]["repeat_count"] == ""
    assert resolved[0]["status"] == "detected_unresolved"
    assert diagnostics[0]["microassembly_success"] == "no"
