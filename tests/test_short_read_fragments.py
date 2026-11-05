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
from mlvamaps.short_read_fragments import FlankRecruiter, RecruitedFragment, _graph_fragments, assemble_fragments


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


def test_panel_primer_recruitment_rejects_shared_tiny_and_repeat_only_anchors():
    locus, product = locus_and_product()
    primer_only = replace(locus, left_flank_sequence="", right_flank_sequence="")
    recruiter = FlankRecruiter([], loci=[primer_only, replace(primer_only, locus_id="other")])
    recruiter.add(pair("shared", product, revcomp(product)))
    assert not recruiter.fragments
    tiny = replace(primer_only, forward_primer="ACGT", reverse_primer="TGCA")
    repetitive = replace(primer_only, forward_primer=locus.repeat_motif * 2,
                         reverse_primer=revcomp(locus.repeat_motif * 2))
    for panel_locus in (tiny, repetitive):
        recruiter = FlankRecruiter([], loci=[panel_locus])
        recruiter.add(pair("repeat", locus.repeat_motif * 8, revcomp(locus.repeat_motif * 8)))
        assert not recruiter.fragments and not recruiter.seeds


def test_short_panel_primer_inside_another_locus_flank_is_not_unique():
    locus, product = locus_and_product()
    short_primer = locus.forward_primer[:15]
    missing = Locus("missing", forward_primer=short_primer)
    recruiter = FlankRecruiter(generate_candidate_contexts([locus]), loci=[locus, missing])
    assert short_primer not in recruiter.seeds and revcomp(short_primer) not in recruiter.seeds
    assert all(len(seed) == 21 for seed in recruiter.seeds)
    recruiter.add(ReadPair("partial", ReadRecord("partial", short_primer)))
    assert not recruiter.fragments
    recruiter.add(pair("complete", product, revcomp(product)))
    assert set(recruiter.fragments) == {locus.locus_id}


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native tools unavailable")
@pytest.mark.parametrize("primer_only", [False, True])
@pytest.mark.parametrize("no_assembly", [False, True])
def test_panel_locus_without_reference_is_measured_and_empty_locus_retained(tmp_path, primer_only, no_assembly):
    from mlvamaps.cli import main
    from mlvamaps.candidate_contexts import write_candidate_contexts
    from mlvamaps.io import read_loci
    from mlvamaps.mapping_classification import observed_reference_contexts
    from mlvamaps.minimap_mapping import build_competitive_indexes
    from test_fastq_assembly_concordance import _sample

    records, panel, _assembly, expected = _sample(tmp_path)
    loci = read_loci(panel)[:3]
    if primer_only:
        loci = [replace(locus, left_flank_sequence="", right_flank_sequence="") for locus in loci]
    write_tsv([vars(locus) for locus in loci], panel, list(vars(loci[0])))
    database = tmp_path / "database"
    database.mkdir()
    (database / f"{loci[0].locus_id}.fasta").write_text(f">reference\n{records[0][1]}\n")
    contexts = observed_reference_contexts(database, loci)[1]
    assert {context.locus_id for context in contexts} == {loci[0].locus_id}
    assets = write_candidate_contexts(contexts, database / "competitive_mapping")
    build_competitive_indexes(assets["fasta"], database / "competitive_mapping")
    (database / "manifest.json").write_text("{}")
    reads = tmp_path / "reads.fastq"
    write_fastq((ReadRecord(f"{name}_{copy}", sequence, "I" * len(sequence))
                 for name, sequence in records[:2] for copy in range(4)), reads)
    output = tmp_path / "output"
    args = ["call", "-i", "sr", "--fq1", str(reads), "-p", str(panel),
            "--database", str(database), "--sample-id", "s", "--sample-mode", "isolate",
            "-o", str(output), "-t", "2", "--quiet"]
    assert main(args + (["--short-no-assembly"] if no_assembly else [])) == 0
    for filename in ("calls.tsv", "locus_repeat_counts.tsv", "genotype_evidence.tsv"):
        with (output / filename).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        calls = {row.get("locus") or row["locus_id"]: row for row in rows}
        assert len(rows) == len(loci) and set(calls) == {locus.locus_id for locus in loci}
        for locus in loci[:2]:
            assert calls[locus.locus_id]["repeat_count"] == expected[locus.locus_id]
            assert calls[locus.locus_id]["call_accepted"] == "True"
        assert calls[loci[2].locus_id]["repeat_count"] == ""
        assert calls[loci[2].locus_id]["call_accepted"] == "False"
    with (output / "common_locus_calls.tsv").open() as handle:
        calls = {row["locus"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert calls[loci[1].locus_id]["molecule_support"] == calls[loci[1].locus_id]["presence_support"] == "4"
    assert calls[loci[2].locus_id]["presence_support"] == "0"


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


def test_deep_tiling_retains_complete_support_and_is_order_independent():
    locus, product = locus_and_product(12, imperfect=True)
    fragments = recruited_tiling(locus, product, copies=100)
    original_ids = set(fragments)
    assert len(fragments) > 2000
    contexts = generate_candidate_contexts([locus])
    diagnostic = {}
    result = assemble_fragments(fragments, contexts, diagnostics=diagnostic)
    assert result == (product, original_ids, "assembled")
    assert set(fragments) == original_ids
    assert 0 < diagnostic["microassembly_selected_fragment_count"] < 2000
    assert diagnostic["microassembly_selection"] == "exact_sequence_representatives"
    assert assemble_fragments(dict(reversed(list(fragments.items()))), contexts) == result


def test_deep_tiling_cannot_hide_a_rare_supported_competing_branch():
    locus, product = locus_and_product(8, imperfect=True)
    altered = product[:100] + ("C" if product[100] != "C" else "A") + product[101:]
    fragments = recruited_tiling(locus, product, copies=100)
    # Only two independent copies support the minority branch. A random or
    # depth-ranked subset could omit it and falsely accept the major product.
    fragments.update({"alt_" + name: fragment for name, fragment in
                      recruited_tiling(locus, altered, copies=2).items()})
    assert len(fragments) > 2000
    assert assemble_fragments(fragments, generate_candidate_contexts([locus]))[2] == "competing_paths"


def test_graph_representatives_keep_quality_ranked_original_pairs_and_diversity():
    locus, product = locus_and_product(8, imperfect=True)
    fragments = {}
    for name, quality in (("low", "5"), ("best", "J"), ("second", "I"), ("tied", "I")):
        reads = pair(name, product[:100], revcomp(product[-100:]))
        reads = replace(reads, read1=replace(reads.read1, quality=quality * 100),
                        read2=replace(reads.read2, quality=quality * 100))
        fragments[name] = RecruitedFragment(reads, True, True, False)
    fragments["rare_left"] = RecruitedFragment(pair("rare_left", product[:90], "N" * 100), True, False, True)
    fragments["rare_right"] = RecruitedFragment(pair("rare_right", "N" * 100, revcomp(product[-90:])), False, True, True)
    selected = _graph_fragments(fragments, maximum=4)
    assert set(selected) == {"best", "second", "rare_left", "rare_right"}
    assert all(selected[name] is fragments[name] for name in selected)
    assert _graph_fragments(dict(reversed(list(fragments.items()))), maximum=4) == selected


def test_graph_reports_diversity_limit_without_sampling_away_distinct_fragments():
    locus, product = locus_and_product(8, imperfect=True)
    rng = random.Random(171)
    fragments = {str(i): RecruitedFragment(
        pair(str(i), "".join(rng.choices("ACGT", k=100)), revcomp(product[-100:])),
        False, True, True) for i in range(2001)}
    diagnostic = {}
    assert assemble_fragments(fragments, generate_candidate_contexts([locus]), diagnostics=diagnostic) == (
        "", set(), "fragment_diversity_limit")
    assert diagnostic == {"microassembly_selected_fragment_count": 0,
                          "microassembly_selection": "fragment_diversity_limit"}


@pytest.mark.parametrize("count,read_length", [(5, 100), (8, 150)])
def test_single_read_overlaps_resolve_repeats_longer_than_fixed_graph_words(count, read_length):
    locus, product = locus_and_product(count)
    contexts = generate_candidate_contexts([locus])
    recruiter = FlankRecruiter(contexts)
    for start in range(len(product) - read_length + 1):
        sequence = product[start:start + read_length]
        for copy in range(3):
            name = f"{start}_{copy}"
            recruiter.add(ReadPair(name, ReadRecord(name, sequence, "I" * len(sequence))))
    sequence, support, reason = assemble_fragments(recruiter.fragments["L"], contexts)
    assert reason == "assembled"
    assert sequence == product
    assert len(support) >= 3
    assert not any(fragment.rescued for fragment in recruiter.fragments["L"].values())


def test_reconstruction_preserves_sample_flank_variants():
    locus, product = locus_and_product(12, imperfect=True)
    product = product[:35] + ("A" if product[35] != "A" else "C") + product[36:]
    sequence, support, reason = assemble_fragments(
        recruited_tiling(locus, product), generate_candidate_contexts([locus]))
    assert reason == "assembled"
    assert sequence == product
    assert len(support) >= 3


def test_microassembly_preserves_existing_molecule_measurements(monkeypatch):
    from mlvamaps.alignment_evidence import CandidateEvidence
    from mlvamaps.allele_inference import InferenceThresholds
    from mlvamaps.short_read_fragments import resolve_fragment_calls
    locus, product = locus_and_product(5)
    direct = CandidateEvidence("L", 5, "direct", measured_repeat_count=5, raw_repeat_count=5,
                               repeat_sigma=.08, direct_product_measurement=True, full_repeat_span=True,
                               metadata={"product_sequence": product})
    pool = {name: RecruitedFragment(pair(name, product, revcomp(product)), True, True, False)
            for name in ("direct", "partial_a", "partial_b")}
    calls, molecule_calls = infer_alleles([direct], [locus], [], "s", "illumina", require_measurement=True)
    monkeypatch.setattr("mlvamaps.short_read_fragments.assemble_fragments",
                        lambda *args: (product, set(pool), "assembled"))
    monkeypatch.setattr("mlvamaps.mapped_reconstruction.call_query_products", lambda *args:
                        {"L": ({"status": "PASS", "repeat_count": 8, "repeat_count_raw": 8}, product)})
    _calls, evidence, molecule_calls, _diagnostics, _products = resolve_fragment_calls(
        calls, [direct], molecule_calls, [], [], [locus], {"L": pool}, "s", InferenceThresholds(), .25,
        require_measurement=True, sample_mode="isolate")
    assert next(row for row in evidence if row.molecule_id == "direct") is direct
    assert molecule_calls == {("L", "direct"): 5}
    assert {row.molecule_id for row in evidence if row.evidence_tier == "microassembly"} == {"partial_a", "partial_b"}


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native tools unavailable")
def test_short_read_outputs_retain_uninformative_depth_and_report_mixed_loci(tmp_path):
    from mlvamaps.short_reads import run_short_read_call
    locus, product = locus_and_product(5)
    _locus, alternate = locus_and_product(8)
    panel = tmp_path / "panel.tsv"
    write_tsv([vars(locus)], panel, list(vars(locus)))
    database = tmp_path / "database"
    database.mkdir()
    (database / "L.fasta").write_text(f">reference\n{product}\n")
    reads = tmp_path / "reads.fastq"
    sequences = [product] * 6 + [alternate] * 4 + [product[:70]] * 5
    write_fastq((ReadRecord(f"m{i}", sequence, "I" * len(sequence)) for i, sequence in enumerate(sequences)), reads)
    paths = run_short_read_call(str(reads), None, str(panel), str(tmp_path / "out"), "s",
        database_path=str(database), threads=2, show_progress=False, sample_mode="metagenome")
    with paths["calls"].open() as handle:
        call = next(csv.DictReader(handle, delimiter="\t"))
    assert call["status"] == "MULTIPLE_VARIANTS"
    assert int(call["read_depth"]) == 15
    assert int(call["informative_molecule_count"]) == 10
    assert int(call["uninformative_locus_reads"]) == 5
    with paths["sample_summary"].open() as handle:
        assert next(csv.DictReader(handle, delimiter="\t"))["mixed_loci"] == "1"


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native tools unavailable")
@pytest.mark.parametrize("layout", ["single", "orphan_r1", "orphan_r2", "paired_flank_snp", "mapping_only",
                                   "orphan_r1_mapping_only", "orphan_r2_mapping_only"])
def test_reconstructed_short_read_fingerprint_matches_assembly(tmp_path, layout, monkeypatch):
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.short_reads import run_short_read_call

    count = 12 if layout == "paired_flank_snp" else 5
    locus, product = locus_and_product(count, imperfect=layout == "paired_flank_snp")
    if "mapping_only" in layout:
        # A sample primer SNP changes the graph's terminal seed. Mapping must
        # locate the observed end without supplying the reference nucleotide.
        product = product[:10] + ("A" if product[10] != "A" else "C") + product[11:]
        monkeypatch.setattr(FlankRecruiter, "add", lambda *args: None)
    if layout == "paired_flank_snp":
        product = product[:35] + ("A" if product[35] != "A" else "C") + product[36:]
        fragments = [item.pair for item in recruited_tiling(locus, product).values()]
    else:
        fragments = []
        for start in range(len(product) - 99):
            for copy in range(3):
                name = f"{start}_{copy}"
                fragment = pair(name, product[start:start + 100], "N" * 100)
                if layout.startswith("orphan_r2"):
                    fragment = pair(name, "N" * 100, revcomp(product[start:start + 100]))
                    fragment = replace(fragment, read1=replace(fragment.read1, quality="!" * 100))
                else:
                    fragment = replace(fragment, read2=replace(fragment.read2, quality="!" * 100))
                fragments.append(fragment)
    first, second, panel = (tmp_path / name for name in ("r1.fastq", "r2.fastq", "panel.tsv"))
    write_fastq((item.read1 for item in fragments), first)
    write_fastq((item.read2 for item in fragments), second)
    write_tsv([vars(locus)], panel, list(vars(locus)))
    database = tmp_path / "database"
    database.mkdir()
    # Deliberately different reference length: reconstruction must recover the
    # sample allele rather than improve agreement by copying a database allele.
    _, reference = locus_and_product(4)
    (database / "L.fasta").write_text(f">reference\n{reference}\n")
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(f">contig\n{product}\n")
    assembled = run_assembly_call(
        assembly_path=str(assembly), loci_path=str(panel), sample_id="s",
        outdir=str(tmp_path / "assembly"), threads=2,
    )
    paths = run_short_read_call(
        str(first), None if layout in {"single", "mapping_only"} else str(second), str(panel),
        str(tmp_path / "short"), "s", database_path=str(database), threads=2,
        keep_intermediates=True, show_progress=False, sample_mode="isolate",
        short_min_pair_retention=0,
    )
    assert paths["fingerprint"].read_text() == assembled["fingerprint"].read_text()
    with paths["calls"].open() as handle:
        call = next(csv.DictReader(handle, delimiter="\t"))
    assert float(call["repeat_count"]) == count
    assert call["status"] == "PASS"
    with paths["short_read_diagnostics"].open() as handle:
        diagnostic = next(csv.DictReader(handle, delimiter="\t"))
    assert diagnostic["microassembly_success"] == "yes"
    # Two tiled endpoints (three molecules each) now resolve the repeat from
    # calibrated partial flanks; none individually spans the primer product.
    assert diagnostic["spanning_fragment_count"] == ("0" if layout == "paired_flank_snp" else "6")
    assert product == "".join(paths["short_read_microassembly"].read_text().splitlines()[1:])
    assert product == "".join(paths["query_loci"].read_text().splitlines()[1:])
    if "mapping_only" in layout:
        assert int(diagnostic["mapping_recruited_fragment_count"]) > 0
        assert diagnostic["flank_recruited_fragment_count"] == "0"
    if layout.startswith("orphan"):
        with paths["short_read_qc"].open() as handle:
            qc = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        assert int(qc["orphan_reads"]) == len(fragments)


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
    # Flank-only measurements no longer bypass full-product reconstruction.
    assert diagnostic["microassembly_attempted"] == "yes"


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


def test_inconclusive_graph_retains_candidate_estimate_with_low_confidence():
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
    assert resolved[0]["repeat_count"] == 6
    assert resolved[0]["status"] == "estimated"
    assert resolved[0]["confidence"] == "low"
    assert resolved[0]["best_probability"] == calls[0]["best_probability"]
    assert resolved[0]["candidate_distribution"] == calls[0]["candidate_distribution"]
    assert "candidate range boundary" in resolved[0]["estimation_reason"]
    assert diagnostics[0]["microassembly_success"] == "no"
    assert diagnostics[0]["final_repeat_count"] == 6
    assert diagnostics[0]["call_method"] == "candidate_likelihood"


def test_ambiguous_candidate_estimate_is_retained_without_inventing_missing_calls():
    from mlvamaps.alignment_evidence import CandidateEvidence
    from mlvamaps.allele_inference import InferenceThresholds
    from mlvamaps.short_read_fragments import resolve_fragment_calls
    from mlvamaps.unified_fastq import common_calls_to_compatibility
    locus, _ = locus_and_product()
    missing = replace(locus, locus_id="missing")
    contexts = generate_candidate_contexts([locus, missing])
    evidence = [CandidateEvidence("L", count, "m", left_boundary_span=True, technology="illumina")
                for count in (5, 6)]
    calls, molecules = infer_alleles(evidence, [locus, missing], contexts, "s", "illumina")
    assert calls[0]["status"] == "ambiguous"
    resolved, _, _, _, _ = resolve_fragment_calls(
        calls, evidence, molecules, [], contexts, [locus, missing], {}, "s", InferenceThresholds(), .25)
    estimate, absent = resolved
    assert estimate["repeat_count"] == 5
    assert estimate["status"] == "estimated"
    assert estimate["best_probability"] == .5
    assert absent["repeat_count"] == "" and absent["status"] == "not_found"
    row = common_calls_to_compatibility(resolved)[0]
    assert row["repeat_count"] == 5 and row["status"] == "ESTIMATED"
    assert row["allele_confidence"] == .5
    assert "not a calibrated probability" in row["evidence"]


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native tools unavailable")
def test_failed_assembly_estimate_remains_diagnostic_in_report_and_summary(tmp_path):
    from mlvamaps.short_reads import run_short_read_call
    locus, product = locus_and_product(24)
    fragments = recruited_tiling(locus, product)
    first, second, panel = (tmp_path / name for name in ("r1.fastq", "r2.fastq", "panel.tsv"))
    write_fastq((item.pair.read1 for item in fragments.values()), first)
    write_fastq((item.pair.read2 for item in fragments.values()), second)
    write_tsv([vars(locus)], panel, list(vars(locus)))
    paths = run_short_read_call(str(first), str(second), str(panel), str(tmp_path / "out"), "s",
                                threads=2, show_progress=False)
    with paths["calls"].open() as handle:
        call = next(csv.DictReader(handle, delimiter="\t"))
    assert call["repeat_count"] == "" and call["status"] == "ESTIMATED"
    assert call["predicted_repeat_count"] != "" and call["call_accepted"] == "False"
    assert call["allele_confidence"] != "" and call["allele_distribution"] != ""
    with paths["fingerprint"].open() as handle:
        assert next(csv.DictReader(handle, delimiter="\t"))["L"] == call["repeat_count"]
    with paths["sample_summary"].open() as handle:
        summary = next(csv.DictReader(handle, delimiter="\t"))
    assert summary["partial_loci"] == "1" and summary["complete_loci"] == "0"
    report = paths["report"].read_text()
    assert "Withheld" in report and "Diagnostic prediction (not called)" in report
    assert call["call_rejection_reason"] in report
