import pysam

from mlvamaps.locus_measurement import (
    extract_reference_interval_from_original_read,
    measure_locus_product,
    original_read_in_locus_orientation,
    reference_interval_to_query_interval,
)
from mlvamaps.models import Locus
from mlvamaps.sequence import revcomp


def _alignment(cigar, query, reference_start=0, reverse=False, supplementary=False):
    alignment = pysam.AlignedSegment()
    alignment.query_name = "read"
    alignment.query_sequence = query
    alignment.flag = (16 if reverse else 0) | (2048 if supplementary else 0)
    alignment.reference_id = 0
    alignment.reference_start = reference_start
    alignment.mapping_quality = 60
    alignment.cigarstring = cigar
    return alignment


def _locus():
    return Locus(
        locus_id="L1",
        forward_primer="AAAACCCC",
        reverse_primer="GGGGTTTT",
        left_flank_sequence="GCTA",
        right_flank_sequence="TAGC",
        repeat_motif="ATGC",
        repeat_unit_length_bp=4,
        expected_min_repeats=3,
        expected_max_repeats=8,
        expected_product_size_bp=44,
        nominal_repeat_units=5,
    )


def _product(repeats=5):
    return "AAAACCCCGCTA" + ("ATGC" * repeats) + "TAGCAAAACCCC"


def test_reference_projection_perfect_and_soft_clipped():
    perfect = _alignment("10M", "A" * 10, reference_start=5)
    clipped = _alignment("3S10M2S", "CCC" + "A" * 10 + "GG", reference_start=5)
    assert reference_interval_to_query_interval(perfect, 5, 15) == (0, 10)
    assert reference_interval_to_query_interval(clipped, 5, 15) == (3, 13)


def test_reference_projection_includes_insertions_and_rejects_deleted_boundary():
    insertion = _alignment("5M2I5M", "A" * 12)
    deletion = _alignment("4M2D6M", "A" * 10)
    assert reference_interval_to_query_interval(insertion, 0, 10) == (0, 12)
    assert reference_interval_to_query_interval(deletion, 4, 10) is None


def test_measurement_uses_anchors_in_original_soft_clipped_read():
    product = _product(5)
    sequence = "GGG" + product + "TT"
    measurement = measure_locus_product(sequence, _locus(), "I" * len(sequence), source="fastq_read")
    assert measurement.status == "FULL_PRODUCT"
    assert measurement.product_start == 3
    assert measurement.product_sequence == product
    assert measurement.repeat_length_bp == 20
    assert measurement.raw_repeat_count == 5
    assert measurement.called_allele == 5


def test_repeat_unit_insertion_changes_measured_not_recruitment_allele():
    measurement = measure_locus_product(_product(6), _locus(), source="fastq_read")
    assert measurement.raw_repeat_count == 6
    assert measurement.called_allele == 6
    assert 6 in measurement.allele_likelihoods


def test_partial_flank_spanning_read_is_repeat_informative():
    sequence = "GCTA" + ("ATGC" * 5) + "TAGC"
    measurement = measure_locus_product(sequence, _locus(), source="fastq_read")
    assert measurement.status == "REPEAT_INFORMATIVE"
    assert measurement.called_allele == 5


def test_reverse_read_and_quality_are_synchronized():
    product = _product(5)
    original = revcomp(product)
    quality = "".join(chr(33 + index) for index in range(len(product)))
    alignment = _alignment(f"{len(product)}M", product, reverse=True)
    oriented, oriented_quality = original_read_in_locus_orientation(alignment, original, quality)
    assert oriented == product
    assert oriented_quality == quality[::-1]
    measurement = measure_locus_product(oriented, _locus(), oriented_quality, source="fastq_read")
    assert measurement.status == "FULL_PRODUCT"
    assert measurement.product_quality == quality[::-1]


def test_incomplete_anchor_never_becomes_full_product():
    measurement = measure_locus_product(_product(5)[:-8], _locus(), source="fastq_read")
    assert measurement.status == "REPEAT_INFORMATIVE"
    assert measurement.called_allele == 5


def test_iupac_primer_anchor_is_supported():
    locus = Locus(**{**_locus().__dict__, "forward_primer": "AAAACCCN"})
    measurement = measure_locus_product(_product(5), locus, source="fastq_read")
    assert measurement.status == "FULL_PRODUCT"
    assert measurement.forward_anchor.edit_distance == 0


def test_product_extraction_keeps_softclips_and_quality_in_sync():
    sequence = "CCC" + ("A" * 10) + "GG"
    quality = "0123456789ABCDE"
    alignment = _alignment("3S10M2S", sequence, reference_start=5)
    extracted = extract_reference_interval_from_original_read(
        alignment, sequence, quality, 5, 15, padding=2
    )
    assert (extracted.start, extracted.end) == (1, 15)
    assert extracted.sequence == sequence[1:15]
    assert extracted.quality == quality[1:15]
