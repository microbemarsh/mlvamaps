from dataclasses import replace
import gzip
import shutil

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import infer_alleles, estimate_unresolved_calls
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.io import read_fastq
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.short_read_evidence import estimate_insert_distribution, extract_short_read_evidence
from mlvamaps.short_read_fragments import RecruitedFragment


def test_conflicting_mates_keep_their_alleles_and_independent_molecule_count():
    from mlvamaps.sequence import revcomp
    locus = Locus("L", forward_primer="ACGTTAGCTAGTCGTACGAT", reverse_primer="TGCATCAGTCGATGCTAGCA",
        left_flank_sequence="GACTGACCTGATCGACGTTAGCTA", right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC",
        repeat_motif="ATGC", repeat_unit_length_bp=4)
    fragments = {}
    for name, alleles in (("a", (4, 8)), ("b", (6, 8))):
        reads = []
        for mate, allele in enumerate(alleles, 1):
            sequence = (locus.forward_primer + locus.left_flank_sequence + "ATGC" * allele
                        + locus.right_flank_sequence + revcomp(locus.reverse_primer))
            sequence = sequence if mate == 1 else revcomp(sequence)
            reads.append(ReadRecord(f"{name}/{mate}", sequence, "I" * len(sequence)))
        fragments[name] = RecruitedFragment(ReadPair(name, *reads), True, True, False)
    evidence = extract_short_read_evidence([], [], [locus], recruited_fragments={"L": fragments})
    assert {(row.molecule_id, row.measured_repeat_count) for row in evidence if row.measured_repeat_count is not None} == {
        ("a", 4), ("a", 8), ("b", 6), ("b", 8)}
    calls, molecule_calls = infer_alleles(evidence, [locus], [], "s", "illumina", require_measurement=True)
    estimate_unresolved_calls(calls, evidence, [])
    assert not molecule_calls
    assert calls[0]["repeat_count"] == 8
    assert calls[0]["presence_support"] == calls[0]["estimate_support"] == 2
    assert calls[0]["status"] == "estimated"


def test_partial_bound_uses_observed_nonrepeat_length_for_short_boundary_anchors():
    from mlvamaps.short_read_evidence import partial_repeat_lower_bound
    locus = Locus("L", forward_primer="A" * 20, reverse_primer="T" * 20,
                  left_flank_sequence="GACTGACCTGATCGACGTTAGCTA", right_flank_sequence="C" * 24,
                  repeat_motif="ATGC", repeat_unit_length_bp=4,
                  expected_product_size_bp=140, nominal_repeat_units=7)
    partial = locus.left_flank_sequence + "ATGC" * 12
    # Observed nonrepeat span is 120 bp; 48 bp visible repeat gives allele 14.
    assert partial_repeat_lower_bound(partial, locus, reference_nonrepeat_bp=(120,)) == 14
    assert partial_repeat_lower_bound(partial, locus, reference_nonrepeat_bp=(120, 124)) is None


@pytest.mark.parametrize("phase", range(4))
def test_right_anchored_partial_bounds_handle_motif_phase_and_calibration(phase):
    from mlvamaps.short_read_evidence import partial_repeat_lower_bound
    locus = Locus("L", forward_primer="A" * 20, reverse_primer="T" * 20,
                  left_flank_sequence="GACTGACCTGATCGACGTTAGCTA",
                  right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC", repeat_motif="ATGC",
                  repeat_unit_length_bp=4, expected_product_size_bp=140, nominal_repeat_units=7)
    repeat = ("ATGC" * 14)[phase:phase + 49]
    sequence = repeat + locus.right_flank_sequence
    assert partial_repeat_lower_bound(sequence, locus, (120,)) == 14.25
    assert partial_repeat_lower_bound(sequence, locus, (120, 124)) is None
    assert partial_repeat_lower_bound("N" + sequence[1:], locus, (120,)) is None
    assert partial_repeat_lower_bound(sequence + locus.right_flank_sequence, locus, (120,)) is None


def test_right_anchored_reverse_mate_supplies_one_bound_without_a_measured_count():
    from mlvamaps.sequence import revcomp
    locus = Locus("L", right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC",
                  repeat_motif="ATGC", repeat_unit_length_bp=4)
    sequence = revcomp("ATGC" * 12 + locus.right_flank_sequence)
    read = ReadRecord("m/2", sequence, "I" * len(sequence))
    fragments = {"L": {"m": RecruitedFragment(ReadPair("m", read), False, True, True)}}
    evidence = extract_short_read_evidence([], [], [locus], recruited_fragments=fragments)
    assert len(evidence) == 1 and evidence[0].metadata["repeat_lower_bound"] == 12
    assert evidence[0].measured_repeat_count is None
    calls, molecules = infer_alleles(evidence * 3, [locus], [], "s", "illumina", require_measurement=True)
    assert not molecules and calls[0]["presence_support"] == 1
    assert calls[0]["repeat_lower_bound"] == 12


@pytest.mark.parametrize("cigar,delta", [("50M1I50M", .25), ("50M3D47M", -.75),
    ("39M4I61M", 0), ("61M4I39M", 0), ("35M4D61M", 0), ("59M3D38M", -.25)])
def test_repeat_indels_preserve_partial_units_and_exclude_flank_bases(cigar, delta):
    from mlvamaps.short_read_evidence import _repeat_indel_delta
    context = CandidateContext("c", "L", "A" * 100, 5, 4, 40, 60, "", "")
    row = CandidateAlignment("m", "m", None, "L", "c", 5, "r", 100, 60, 1, 1, 1,
                             0, 100, 0, 100, cigar, "", True, False, False)
    assert _repeat_indel_delta(row, context) == delta


def test_same_count_reference_backgrounds_keep_competing_pair_lengths():
    locus = Locus("L", expected_product_size_bp=200, nominal_repeat_units=5, repeat_unit_length_bp=4)
    contexts = [CandidateContext(name, "L", "A" * 200, 5, 4, 90, 110, "", "") for name in ("a", "b")]
    first = CandidateAlignment("m", "m/1", 1, "L", "a", 5, "ref", 100, 60, 1, 1, .2,
                               0, 40, 20, 60, "40M", "", True, False, False)
    rows = [first, replace(first, mate=2, reverse=True, reference_start=120, reference_end=160),
            replace(first, candidate_id="b"),
            replace(first, candidate_id="b", mate=2, reverse=True, reference_start=140, reference_end=180)]
    evidence = extract_short_read_evidence(rows, contexts, [locus], insert_mean=160, insert_sd=2)
    assert {row.metadata["pair_repeat_count"] for row in evidence} == {5, 10}
    assert {row.molecule_id for row in evidence} == {"m"}


@pytest.mark.skipif(not shutil.which("sassy"), reason="native primer indel calibration required")
@pytest.mark.parametrize("reference_index", [7, 8, 9, 10])
def test_pair_geometry_excludes_reference_primer_indels_outside_the_fragment(tmp_path, reference_index):
    from test_fastq_assembly_concordance import _sample
    from mlvamaps.io import read_loci
    from mlvamaps.sequence import revcomp

    records, panel, _assembly, expected = _sample(tmp_path)
    locus = read_loci(panel)[reference_index]
    sequence = records[reference_index][1]
    if reference_index % 2:
        sequence = revcomp(sequence)
    start = sequence.index(locus.left_flank_sequence) + len(locus.left_flank_sequence)
    end = sequence.index(locus.right_flank_sequence)
    count = float(expected[locus.locus_id])
    context = CandidateContext("c", locus.locus_id, sequence, count, 4, start, end, "", "")
    first = CandidateAlignment("m", "m/1", 1, locus.locus_id, "c", count, "R", 42, 60, 1, 1, .2,
        0, 21, start - 21, start, "21M", "", True, False, False,
        query_sequence=sequence[start - 21:start], query_quality="I" * 21)
    second = replace(first, read_id="m/2", mate=2, reverse=True, reference_start=end, reference_end=end + 21,
                     query_sequence=revcomp(sequence[end:end + 21]))
    evidence = extract_short_read_evidence([first, second], [context], [locus],
        insert_mean=end + 21 - (start - 21), insert_sd=2)
    assert evidence[0].measured_repeat_count is None
    assert evidence[0].metadata["pair_repeat_count"] == count == 5


def test_insert_calibration_uses_independent_flank_pairs_and_rejects_repeat_geometry():
    context = CandidateContext("c", "L", "A" * 2000, 5, 4, 1000, 1020, "", "")
    first = CandidateAlignment("", "", 1, "L", "c", 5, "ref", 100, 60, 1, 1, 1,
                               0, 100, 200, 300, "100M", "", True, False, False)
    rows = [replace(first, molecule_id=str(i)) for i in range(120)]
    rows += [replace(first, molecule_id=str(i), mate=2, reverse=True,
                     reference_start=400 + i % 20, reference_end=500 + i % 20) for i in range(120)]
    diagnostics = {}
    fit = estimate_insert_distribution(rows, [context], diagnostics=diagnostics)
    assert fit[0] == 309.5 and 5 < fit[1] < 7 and fit[2] == 120
    assert diagnostics == dict(calibration_status="calibrated", eligible_fragments=120,
                               retained_fragments=120, required_fragments=100)
    assert estimate_insert_distribution(rows[:10], [context]) is None
    assert estimate_insert_distribution(rows, [replace(context, repeat_start=350)]) is None
    assert estimate_insert_distribution([replace(row, mapping_quality=0) for row in rows], [context]) is None
    assert estimate_insert_distribution([], [context], diagnostics=diagnostics) is None
    assert diagnostics["calibration_status"] == "no_alignments" and diagnostics["eligible_fragments"] == 0
    assert estimate_insert_distribution(rows[:10], [context], diagnostics=diagnostics) is None
    assert diagnostics["calibration_status"] == "insufficient_unique_flank_pairs"
    constant = [replace(row, reference_start=400, reference_end=500) if row.mate == 2 else row for row in rows]
    assert estimate_insert_distribution(constant, [context], diagnostics=diagnostics) is None
    assert diagnostics["calibration_status"] == "zero_dispersion"
    outliers = [replace(row, reference_start=row.reference_start + 300, reference_end=row.reference_end + 300)
                if row.mate == 2 and int(row.molecule_id) >= 90 else row for row in constant]
    assert estimate_insert_distribution(outliers, [context], diagnostics=diagnostics) is None
    assert diagnostics["calibration_status"] == "insufficient_pairs_after_trimming"
    assert diagnostics["eligible_fragments"] == 120 and diagnostics["retained_fragments"] == 90
    edge = [replace(row, reference_start=row.reference_start - 200, reference_end=row.reference_end - 200) for row in rows]
    assert estimate_insert_distribution(edge, [context], diagnostics=diagnostics) is None
    assert diagnostics["calibration_status"] == "target_edge_censoring"


@pytest.mark.parametrize("insert", [None, (350, 30)])
def test_library_calibration_export_distinguishes_unavailable_and_supplied(tmp_path, monkeypatch, insert):
    import csv
    from mlvamaps.unified_fastq import run_unified_fastq_inference
    monkeypatch.setattr("mlvamaps.unified_fastq.map_reads_to_candidates_bam", lambda *args, **kwargs: [])
    locus = Locus("L", repeat_unit_length_bp=4, expected_min_repeats=4, expected_max_repeats=5)
    context = CandidateContext("c", "L", "ACGT" * 30, 4, 4, 50, 66, "", "")
    calls, _, _, paths = run_unified_fastq_inference(
        reads1=tmp_path / "reads.fastq", reads2=None, loci=[locus], contexts=[context],
        database_path=None, outdir=tmp_path, sample_id="s", technology="illumina",
        minimap2_bin="minimap2", threads=1, minimum_molecules=3, minimum_probability=.8,
        insert_mean=insert[0] if insert else None, insert_sd=insert[1] if insert else None)
    with paths["library_insert_distribution"].open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["calibration_status"] == ("user_supplied" if insert else "no_alignments")
    assert row["source"] == ("user_supplied" if insert else "unavailable")
    assert row["mean_bp"] == ("350" if insert else "")
    assert row["eligible_fragments"] == ("" if insert else "0")
    assert calls[0]["prior_only"] and calls[0]["independent_support"] == 0


def test_unspanned_long_array_reports_a_bound_and_never_copies_the_reference_count():
    locus = Locus("L", left_flank_sequence="GACTGACCTGATCGACGTTAGCTA", repeat_motif="ATGC", repeat_unit_length_bp=4)
    read = ReadRecord("m", locus.left_flank_sequence + "ATGC" * 12, "I" * 72)
    fragments = {"L": {"m": RecruitedFragment(ReadPair("m", read), True, False, False)}}
    evidence = extract_short_read_evidence([], [], [locus], recruited_fragments=fragments)
    calls, molecules = infer_alleles(evidence, [locus], [], "s", "illumina", require_measurement=True)
    estimate_unresolved_calls(calls, evidence, [])
    assert not molecules and calls[0]["repeat_count"] == ""
    assert calls[0]["repeat_lower_bound"] == 12
    assert evidence[0].repeat_count is None


def test_reconstructed_product_does_not_multiply_likelihood_or_classification_votes():
    def run(n):
        evidence = [CandidateEvidence("L", 5, str(i), measured_repeat_count=5, raw_repeat_count=5,
                                      repeat_sigma=.3, evidence_tier="microassembly") for i in range(n)]
        return infer_alleles(evidence, [Locus("L", repeat_unit_length_bp=4)], [], "s", "illumina", require_measurement=True)
    three, molecules = run(3)
    thirty, more = run(30)
    assert not molecules and not more
    assert three[0]["best_probability"] == thirty[0]["best_probability"]
    assert thirty[0]["derived_molecule_support"] == 30


def test_ambiguous_boundaries_survive_fragment_extraction():
    from mlvamaps.sequence import revcomp
    locus = Locus("L", forward_primer="ACGTTAGCTAGTCGTACGAT", reverse_primer="TGCATCAGTCGATGCTAGCA",
        left_flank_sequence="GACTGACCTGATCGACGTTAGCTA", right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC",
        repeat_motif="ATGC", repeat_unit_length_bp=4)
    sequence = (locus.forward_primer + locus.left_flank_sequence + "ATGC" * 2
                + locus.left_flank_sequence + "ATGC" * 5 + locus.right_flank_sequence + revcomp(locus.reverse_primer))
    read = ReadRecord("m", sequence, "I" * len(sequence))
    fragments = {"L": {"m": RecruitedFragment(ReadPair("m", read), True, True, False)}}
    evidence = extract_short_read_evidence([], [], [locus], recruited_fragments=fragments)
    assert evidence[0].measured_repeat_count is not None
    assert "boundary_ambiguous" in evidence[0].quality_flags
    assert evidence[0].repeat_sigma >= 1


@pytest.mark.parametrize("content", [b"@\nACGT\n+\nIIII\n", b"@x\n\n+\n\n", b"@x\nACGT\n+\nIII\n"])
def test_malformed_records_fail_instead_of_creating_missing_genotypes(tmp_path, content):
    path = tmp_path / "reads.fastq"
    path.write_bytes(content)
    with pytest.raises(ValueError, match="FASTQ"):
        list(read_fastq(path))


def test_truncated_gzip_is_an_error_and_empty_fastq_has_zero_records(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(b"@x\nACGT\n+\nIIII\n")[:-8])
    with pytest.raises(ValueError, match="FASTQ"):
        list(read_fastq(path))
    path.write_bytes(gzip.compress(b""))
    assert list(read_fastq(path)) == []


def test_duplicate_fastq_ids_are_rejected_at_ingestion(tmp_path):
    from mlvamaps.io import read_fastq_pairs
    path = tmp_path / "reads.fastq"
    path.write_text("@m\nACGT\n+\nIIII\n@m\nTGCA\n+\nIIII\n")
    with pytest.raises(ValueError, match="Duplicate FASTQ molecule"):
        list(read_fastq_pairs(path))
