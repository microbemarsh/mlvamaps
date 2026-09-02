import csv
from pathlib import Path

from mlvamaps.assembly_call import run_assembly_call
from mlvamaps.in_silico_pcr import (
    _legacy_first_primer_orientations,
    _new_searcher,
)
from mlvamaps.locus_measurement import find_anchor


ORACLE_DIR = Path(__file__).parent / "data" / "mlva_finder_oracle"


def test_mlva_finder_oracle_matches_with_parallel_pcr(tmp_path):
    with (ORACLE_DIR / "expected.csv").open(newline="") as handle:
        expected = {
            row["locus_id"]: row["repeat_count"]
            for row in csv.DictReader(handle)
        }
    result = run_assembly_call(
        assembly_path=str(ORACLE_DIR / "assembly.fasta"),
        loci_path=None,
        primers_path=str(ORACLE_DIR / "primers.tsv"),
        outdir=str(tmp_path / "result"),
        sample_id="oracle",
        threads=4,
    )
    with result["calls"].open(newline="") as handle:
        observed = {
            row["locus_id"]: row["repeat_count"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert observed == expected


def test_degenerate_primer_preserves_per_expansion_reverse_fallback():
    # AAR expands to AAA and AAG. AAA occurs forward while the reverse
    # complement of AAG (CTT) occurs on the input strand. MLVA_finder retains
    # both because reverse fallback is evaluated separately per expansion.
    sequence = "GGGAAACCCCTTGGG"
    orientations = _legacy_first_primer_orientations(
        _new_searcher(),
        "AAR",
        sequence,
        sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1],
        0,
        True,
    )
    assert {strand for strand, _sequence, _matches in orientations} == {"+", "-"}


def test_native_anchor_search_reports_substitution_and_indels():
    substitution = find_anchor("ACGTACGA", "GGACGTTCGATT", 3)
    insertion = find_anchor("ACGTACGA", "GGACGTTACGATT", 3)
    deletion = find_anchor("ACGTACGA", "GGACGACGATT", 3)
    assert (
        substitution.mismatches,
        substitution.insertions,
        substitution.deletions,
    ) == (1, 0, 0)
    assert (insertion.mismatches, insertion.insertions, insertion.deletions) == (
        0,
        1,
        0,
    )
    assert (deletion.mismatches, deletion.insertions, deletion.deletions) == (
        0,
        0,
        1,
    )
