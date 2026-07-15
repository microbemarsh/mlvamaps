from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _short_locus_name(locus_id: str) -> str:
    parts = locus_id.split("_")
    if len(parts) >= 4 and parts[-1].lower().endswith("u") and parts[-2].lower().endswith("bp"):
        return "_".join(parts[:-3])
    return locus_id


def _clean_value(value: str) -> str:
    value = value.strip()
    if not value:
        return ""
    try:
        numeric = float(value)
    except ValueError:
        return value
    if numeric.is_integer():
        return str(int(numeric))
    return str(numeric).rstrip("0").rstrip(".")


def read_primer_locus_map(primers_path: Path) -> dict[str, str]:
    with primers_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        mapping = {_short_locus_name(row["locus_id"]): row["locus_id"] for row in reader}
    return mapping


def convert_profiles(
    input_path: Path,
    primers_path: Path,
    output_path: Path,
    keep_empty: bool = False,
) -> tuple[int, int, list[str]]:
    locus_map = read_primer_locus_map(primers_path)
    with input_path.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        if not reader.fieldnames or "Access_number" not in reader.fieldnames:
            raise ValueError("Expected an Access_number column in the UF VNTR table.")
        input_loci = [field for field in reader.fieldnames if field != "Access_number"]
        unmatched = [field for field in input_loci if field not in locus_map]
        locus_pairs = [(field, locus_map[field]) for field in input_loci if field in locus_map]
        output_loci = [output_locus for _input_locus, output_locus in locus_pairs]
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="") as out_handle:
            writer = csv.DictWriter(
                out_handle,
                fieldnames=["profile_id", "strain_id", "metadata"] + output_loci,
                delimiter="\t",
            )
            writer.writeheader()
            rows_written = 0
            for row in reader:
                accession = row.get("Access_number", "").strip()
                if not accession:
                    continue
                converted_values = {
                    output_locus: _clean_value(row.get(input_locus, ""))
                    for input_locus, output_locus in locus_pairs
                }
                if not keep_empty and not any(converted_values.values()):
                    continue
                out_row = {
                    "profile_id": accession,
                    "strain_id": accession,
                    "metadata": "UF Ba VNTR database",
                }
                out_row.update(converted_values)
                writer.writerow(out_row)
                rows_written += 1
    return rows_written, len(input_loci), unmatched


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert the UF B. anthracis VNTR table into mlvamaps profiles.")
    parser.add_argument("input", type=Path)
    parser.add_argument("--primers", type=Path, default=Path("examples/seer_lab_Ba/mlvamaps_primers.example.tsv"))
    parser.add_argument("--output", type=Path, default=Path("data/uf_ba_mlva_profiles.tsv"))
    parser.add_argument("--keep-empty", action="store_true", help="Keep accessions with no VNTR values")
    args = parser.parse_args()

    rows_written, input_loci, unmatched = convert_profiles(args.input, args.primers, args.output, args.keep_empty)
    print(f"Wrote {rows_written} profiles across {input_loci - len(unmatched)} mapped loci to {args.output}")
    if unmatched:
        print("Unmatched input loci: " + ", ".join(unmatched))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
