import pytest

from mlvamaps.reference_database import read_sequence_database, canonical_assembly_digest

def test_sequence_database_directory_uses_locus_filenames(tmp_path):
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(">R1\nAAAA\n>R2\nAAAT\n")
    (database / "unrelated.fasta").write_text(">R1\nCCCC\n")
    assert read_sequence_database(database, {"L1"}) == {
        "L1": [("R1", "AAAA"), ("R2", "AAAT")]
    }


def test_empty_sequence_database_remains_an_error(tmp_path):
    database = tmp_path / "database"
    database.mkdir()
    with pytest.raises(ValueError, match="contains no FASTA files"):
        read_sequence_database(database, {"L1"})


def test_canonical_assembly_digest_ignores_headers_order_and_orientation(tmp_path):
    first = tmp_path / "first.fa"
    second = tmp_path / "second.fa"
    first.write_text(">alpha\nAAGC\n>beta\nTTAA\n")
    second.write_text(">renamed_beta\nTTAA\n>renamed_alpha\nGCTT\n")

    assert canonical_assembly_digest(
        first
    ) == canonical_assembly_digest(second)
