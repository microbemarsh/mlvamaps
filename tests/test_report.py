from mlvamaps.models import Locus
from mlvamaps.report import _assembly_gel_svg, _gel_svg


def test_gels_prefer_exact_phylogenetic_reference_amplicon_sizes():
    loci = [
        Locus(
            locus_id="VNTR_01",
            forward_primer="ACGTTGCAAC",
            reverse_primer="TGCATGCAAA",
            left_flank_sequence="GGTA",
            right_flank_sequence="CCAT",
            repeat_motif="ATG",
        )
    ]
    reference_bands = [
        {
            "reference_id": "PHYLO_R1",
            "locus_id": "VNTR_01",
            "product_size_bp": "47",
            "repeat_count": "6",
        }
    ]
    profile = {"profile_id": "PROFILE_R1", "VNTR_01": "5"}

    fastq_gel = _gel_svg(
        "QUERY",
        loci,
        [
            {
                "locus_id": "VNTR_01",
                "called_repeat_count": "5",
                "read_depth": "12",
            }
        ],
        profile,
        [],
        reference_bands,
    )
    assert "PHYLO_R1" in fastq_gel
    assert "VNTR_01 reference: 47 bp" in fastq_gel
    assert ">PROFILE_R1<" not in fastq_gel

    assembly_gel = _assembly_gel_svg(
        "QUERY",
        [
            {
                "locus_id": "VNTR_01",
                "present": "yes",
                "product_size_bp": "44",
                "read_depth": "",
                "repeat_count": "5",
            }
        ],
        profile,
        loci,
        reference_bands,
    )
    assert "PHYLO_R1" in assembly_gel
    assert "VNTR_01 (6U) reference: 47 bp" in assembly_gel
    assert ">PROFILE_R1<" not in assembly_gel
