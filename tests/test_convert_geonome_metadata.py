from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "convert_geonome_metadata.py"
SPEC = importlib.util.spec_from_file_location("convert_geonome_metadata", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_convert_geonome_reference_directory(tmp_path):
    reference = tmp_path / "geonome-reference"
    reference.mkdir()
    (reference / "reference_manifest.json").write_text(
        json.dumps({"artifacts": {"normalized_metadata": "tables/metadata.tsv"}})
    )
    metadata = reference / "tables" / "metadata.tsv"
    metadata.parent.mkdir()
    metadata.write_text(
        "genome_id\tcollection_date\tlatitude\tlongitude\tcountry\tstate\tcounty\t"
        "location_raw\tisolation_source\tsample_type_normalized\n"
        "GCF_000000001.1\t2025-01-02\t40.1\t-75.2\tUSA\tPA\tBerks\tReading\tsoil\tenvironment\n"
        "GCF_000000002.1\t2025\t\t\tUSA\tNY\t\t\t\twater\n"
    )
    output = tmp_path / "reference_metadata.tsv"

    assert MODULE.convert_metadata(reference, output) == 2

    rows = _read_tsv(output)
    assert list(rows[0]) == MODULE.OUTPUT_FIELDS
    assert rows[0] == {
        "reference_id": "GCF_000000001.1",
        "collection_date": "2025-01-02",
        "latitude": "40.1",
        "longitude": "-75.2",
        "location": "Reading",
        "source": "soil",
    }
    assert rows[1]["location"] == "NY, USA"
    assert rows[1]["source"] == "water"


def test_convert_rejects_duplicate_identifiers(tmp_path):
    metadata = tmp_path / "normalized_metadata.tsv"
    metadata.write_text("genome_id\tlatitude\tlongitude\nR1\t\t\nR1\t\t\n")

    try:
        MODULE.convert_metadata(metadata, tmp_path / "output.tsv")
    except ValueError as exc:
        assert "duplicate reference identifier" in str(exc)
    else:
        raise AssertionError("duplicate identifiers should fail conversion")
