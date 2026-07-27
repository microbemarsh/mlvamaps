from mlvamaps.models import Locus, ReadRecord
from mlvamaps.recruitment import (
    build_recruitment_references,
    local_product_records,
    parse_recruitment_sam,
    recruitment_fallback_evidence,
    recruitment_summary_rows,
)


def _locus() -> Locus:
    return Locus(
        locus_id="VNTR",
        forward_primer="AAAACCCC",
        reverse_primer="GGGGTTTT",
        left_flank_sequence="GCTA",
        right_flank_sequence="TAGC",
        repeat_motif="ATGC",
        repeat_unit_length_bp=4,
        expected_min_repeats=5,
        expected_max_repeats=6,
        nominal_repeat_units=5,
        expected_product_size_bp=44,
    )


def test_recruitment_reference_bank_contains_expected_alleles():
    references = build_recruitment_references([_locus()])

    assert [row["candidate_allele"] for row in references] == [
        5,
        5.5,
        5.5,
        5.5,
        6,
    ]
    assert [len(row["sequence"]) for row in references] == [44, 45, 46, 47, 48]
    assert {row["reference_source"] for row in references} == {
        "synthetic_panel_product"
    }


def test_database_product_is_preferred_and_modal_length_selected(tmp_path):
    database = tmp_path / "database"
    database.mkdir()
    canonical = "AAAACCCCGCTA" + ("ATGC" * 5) + "TAGCAAAACCCC"
    (database / "VNTR.fasta").write_text(
        f">short\n{canonical[:-4]}\n>modal1\n{canonical}\n>modal2\n{canonical}\n"
    )

    references = build_recruitment_references([_locus()], database)

    assert references[0]["reference_source"] == "database_product"
    assert len(references[0]["sequence"]) == len(canonical)


def test_competitive_mapping_separates_presence_from_full_product(tmp_path):
    locus = _locus()
    references = build_recruitment_references([locus])
    reference = next(row for row in references if row["candidate_allele"] == 5)
    sequence = reference["sequence"]
    partial = sequence[4:20]
    repeat_spanning = sequence[8:36]
    sam = tmp_path / "recruitment.sam"
    sam.write_text(
        f"@SQ\tSN:{reference['reference_name']}\tLN:{len(sequence)}\n"
        f"full\t0\t{reference['reference_name']}\t1\t60\t{len(sequence)}M\t*\t0\t0\t"
        f"{sequence}\t{'I' * len(sequence)}\tNM:i:0\n"
        f"partial\t0\t{reference['reference_name']}\t5\t40\t{len(partial)}M\t*\t0\t0\t"
        f"{partial}\t{'I' * len(partial)}\tNM:i:0\n"
        f"repeat_partial\t0\t{reference['reference_name']}\t9\t50\t{len(repeat_spanning)}M\t*\t0\t0\t"
        f"{repeat_spanning}\t{'I' * len(repeat_spanning)}\tNM:i:0\n"
    )
    reads = [
        ReadRecord("full", sequence, "I" * len(sequence)),
        ReadRecord("partial", partial, "I" * len(partial)),
        ReadRecord(
            "repeat_partial",
            repeat_spanning,
            "I" * len(repeat_spanning),
        ),
    ]

    rows, assignments = parse_recruitment_sam(
        sam,
        references,
        reads,
        [locus],
        min_aligned_bp=10,
    )
    summaries = recruitment_summary_rows("sample", [locus], rows)

    by_read = {row["read_id"]: row for row in rows}
    assert by_read["full"]["evidence_class"] == "FULL_PRODUCT"
    assert by_read["partial"]["evidence_class"] == "PRESENCE_ONLY"
    assert by_read["repeat_partial"]["evidence_class"] == "REPEAT_INFORMATIVE"
    assert len(assignments) == 1
    assert assignments[0].product_size_bp == len(sequence)
    assert summaries[0]["presence_status"] == "PRESENT_PROVISIONAL"
    assert summaries[0]["mapped_reads"] == 3
    assert summaries[0]["full_product_reads"] == 1
    assert summaries[0]["genotype_informative_reads"] == 2


def test_local_product_uses_modal_length_and_read_consensus(tmp_path):
    locus = _locus()
    references = build_recruitment_references([locus])
    reference = next(row for row in references if row["candidate_allele"] == 5)
    sequence = reference["sequence"]
    altered = sequence[:15] + "T" + sequence[16:]
    sam = tmp_path / "consensus.sam"
    sam.write_text(
        f"@SQ\tSN:{reference['reference_name']}\tLN:{len(sequence)}\n"
        f"r1\t0\t{reference['reference_name']}\t1\t60\t{len(sequence)}M\t*\t0\t0\t"
        f"{sequence}\t{'I' * len(sequence)}\tNM:i:0\n"
        f"r2\t0\t{reference['reference_name']}\t1\t60\t{len(sequence)}M\t*\t0\t0\t"
        f"{sequence}\t{'I' * len(sequence)}\tNM:i:0\n"
        f"r3\t0\t{reference['reference_name']}\t1\t60\t{len(sequence)}M\t*\t0\t0\t"
        f"{altered}\t{'I' * len(sequence)}\tNM:i:1\n"
    )
    reads = [
        ReadRecord("r1", sequence, "I" * len(sequence)),
        ReadRecord("r2", sequence, "I" * len(sequence)),
        ReadRecord("r3", altered, "I" * len(sequence)),
    ]

    _rows, assignments = parse_recruitment_sam(
        sam,
        references,
        reads,
        [locus],
        min_aligned_bp=10,
    )
    products = local_product_records(assignments)

    assert products == [("VNTR_local_primary", sequence)]


def test_repeat_spanning_partial_read_produces_provisional_fallback():
    locus = _locus()
    rows = [
        {
            "read_id": "partial",
            "locus_id": "VNTR",
            "candidate_allele": 5,
            "alignment_identity": 0.99,
            "genotype_informative": "yes",
        }
    ]

    asvs, predictions = recruitment_fallback_evidence(
        rows,
        [locus],
        feature_loci=set(),
        sample_id="sample",
    )

    assert asvs[0]["repeat_count"] == 5
    assert asvs[0]["support_reads"] == 1
    assert predictions[0].predicted_repeat_count == 5
    assert predictions[0].evidence_weight == 0.99


def test_recruitment_requires_a_score_lead_over_other_loci(tmp_path):
    locus = _locus()
    other = Locus(**{**locus.__dict__, "locus_id": "OTHER"})
    reference = next(
        row
        for row in build_recruitment_references([locus])
        if row["candidate_allele"] == 5
    )
    competing = {
        **reference,
        "reference_name": "other_reference",
        "locus_id": "OTHER",
    }
    sequence = reference["sequence"]
    sam = tmp_path / "competitive.sam"
    sam.write_text(
        f"@SQ\tSN:{reference['reference_name']}\tLN:{len(sequence)}\n"
        f"@SQ\tSN:{competing['reference_name']}\tLN:{len(sequence)}\n"
        f"read\t0\t{reference['reference_name']}\t1\t0\t{len(sequence)}M\t*\t0\t0\t"
        f"{sequence}\t{'I' * len(sequence)}\tNM:i:0\tAS:i:40\n"
        f"read\t256\t{competing['reference_name']}\t1\t0\t{len(sequence)}M\t*\t0\t0\t"
        f"{sequence}\t{'I' * len(sequence)}\tNM:i:0\tAS:i:35\n"
    )
    reads = [ReadRecord("read", sequence, "I" * len(sequence))]

    rejected, _assignments = parse_recruitment_sam(
        sam,
        [reference, competing],
        reads,
        [locus, other],
        min_aligned_bp=10,
        min_locus_score_margin=10,
    )
    accepted, _assignments = parse_recruitment_sam(
        sam,
        [reference, competing],
        reads,
        [locus, other],
        min_aligned_bp=10,
        min_locus_score_margin=5,
    )

    assert rejected == []
    assert accepted[0]["locus_id"] == "VNTR"
    assert accepted[0]["locus_score_margin"] == 5
