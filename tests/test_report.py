from mlvamaps.models import Locus
from mlvamaps.report import (
    _assembly_gel_svg,
    _gel_svg,
    _phylogenetic_warning_html,
    write_assembly_report,
    write_report,
)


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
                "primary_product_size_bp": "44",
                "primary_read_depth": "12",
                "read_depth": "12",
            }
        ],
        profile,
        [
            {
                "locus_id": "VNTR_01",
                "repeat_count": "8",
                "support_reads": "12",
                "frequency": "1.0",
            }
        ],
        reference_bands,
    )
    assert "PHYLO_R1" in fastq_gel
    assert "VNTR_01 reference: 47 bp" in fastq_gel
    assert "VNTR_01: 44 bp" in fastq_gel
    assert "VNTR_01: 53 bp" not in fastq_gel
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


def test_phylogenetic_warning_explains_exact_match_override():
    warning = _phylogenetic_warning_html(
        [
            {
                "reference_id": "R1",
                "ranking_warning": "EXACT_MATCH_OVERRIDES_PLACEMENT",
            }
        ]
    )
    assert "Exact-match placement warning" in warning
    assert "R1" in warning
    assert "likelihood-weighted placement distance" in warning


def test_reports_prioritize_sample_findings_and_remove_novelty(tmp_path):
    locus = Locus(locus_id="L1", repeat_motif="AT")
    phylogenetic_rows = [
        {
            "rank": "1",
            "reference_id": "R1",
            "match_status": "EXACT_AMPLICON_MATCH",
            "combined_marker_distance": "0.00000000",
            "whole_genome_exact_match": "yes",
            "whole_genome_snps": "0",
            "whole_genome_indel_bases": "0",
            "whole_genome_align_fraction_ref": "100.00000000",
            "whole_genome_align_fraction_query": "100.00000000",
            "tie_break_status": "APPLIED",
        }
    ]
    write_report(
        tmp_path / "reads",
        "sample",
        [
            {
                "locus_id": "L1",
                "called_repeat_count": "5",
                "posterior_probability": "0.99",
                "read_depth": "20",
                "call_status": "PASS",
            }
        ],
        [locus],
        phylogenetic_rows=phylogenetic_rows,
        local_assembly_rows=[
            {
                "locus_id": "L1",
                "dominant_variant_id": "L1_ASV1",
                "input_reads": "20",
                "unique_sequences": "2",
                "observed_min_product_bp": "99",
                "observed_modal_product_bp": "100",
                "observed_max_product_bp": "101",
                "poa_consensus_bp": "100",
                "pcr_product_size_bp": "100",
                "raw_repeat_count": "5",
                "called_repeat_count": "5",
                "measurement_source": "dominant_cluster_poa_assembly",
                "pcr_status": "PASS",
            }
        ],
    )
    write_assembly_report(
        tmp_path / "assembly",
        "sample",
        [
            {
                "locus_id": "L1",
                "present": "yes",
                "repeat_count": "5",
                "product_size_bp": "100",
                "status": "PASS",
            }
        ],
        [{"locus_id": "L1", "product_size_bp": "100"}],
        loci=[locus],
        phylogenetic_rows=phylogenetic_rows,
    )

    for report_path in (
        tmp_path / "reads" / "report.html",
        tmp_path / "assembly" / "report.html",
    ):
        report = report_path.read_text()
        assert "Sample Overview" in report
        assert "Closest Reference Genomes" in report
        assert "Exact whole-genome match" in report
        assert "Technical marker-distance components" in report
        assert "Novelty" not in report
    fastq_report = (tmp_path / "reads" / "report.html").read_text()
    assert "FASTQ Local Assembly Concordance" in fastq_report
    assert "Raw bp min / mode / max" in fastq_report
    assert "dominant_cluster_poa_assembly" in fastq_report
    assert "POA assembly calls" in fastq_report
    assert "Generated " in fastq_report
    assert "Cache-Control" in fastq_report
