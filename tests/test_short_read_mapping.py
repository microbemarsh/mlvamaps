from __future__ import annotations


from mlvamaps.models import Locus
from mlvamaps.minimap_mapping import minimap2_competitive_command


def _locus(unit: int = 4) -> Locus:
    return Locus(
        "L1", forward_primer="ACGT", reverse_primer="TGCA",
        left_flank_sequence="AACCGGTT", right_flank_sequence="GGTTAACC",
        repeat_motif="ATGC" if unit == 4 else "A" * unit,
        repeat_unit_length_bp=unit, expected_min_repeats=2,
        expected_max_repeats=8, nominal_repeat_units=4,
    )








def test_minimap2_command_is_single_competitive_paired_mapping(tmp_path):
    command = minimap2_competitive_command(
        tmp_path / "contexts.fasta", tmp_path / "r1.fq",
        tmp_path / "r2.fq", 8, "illumina", "/opt/minimap2",
    )
    assert command[:4] == ["/opt/minimap2", "-a", "-x", "sr"]
    assert "--cs=long" in command and "--secondary=yes" in command
    assert command[-2].endswith("r1.fq") and command[-1].endswith("r2.fq")


def test_candidate_contexts_are_structured_and_collapsed(tmp_path):
    from mlvamaps.candidate_contexts import generate_candidate_contexts, write_candidate_contexts
    contexts = generate_candidate_contexts([_locus()], maximum=6)
    first = write_candidate_contexts(contexts, tmp_path)
    assert first["metadata"].is_file()
    assert first["fasta"].is_file()
    assert [context.repeat_count for context in contexts] == [2, 3, 4, 5, 6]
