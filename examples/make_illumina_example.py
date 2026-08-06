from __future__ import annotations

import gzip
import sys
from pathlib import Path


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def main() -> None:
    output = Path(sys.argv[1] if len(sys.argv) > 1 else "illumina_demo")
    output.mkdir(parents=True, exist_ok=True)
    forward = "ACGTTGCAACGTTGCAACGT"
    reverse = "AGTCAGTCAGTCAGTCAGTC"
    left = "CCGGAATTCGACCTGA"
    right = "TTAACCGGCTAGGTCA"
    product = forward + left + "ATGC" * 4 + right + reverse_complement(reverse)
    (output / "panel.tsv").write_text(
        "locus_id\tforward_primer\treverse_primer\tleft_flank_sequence\tright_flank_sequence\trepeat_motif\trepeat_unit_length_bp\texpected_min_repeats\texpected_max_repeats\tnominal_repeat_units\n"
        f"demo\t{forward}\t{reverse}\t{left}\t{right}\tATGC\t4\t2\t10\t4\n"
    )
    with gzip.open(output / "truth.fasta.gz", "wt") as handle:
        handle.write(f">SRR_DEMO\nTTTT{product}CCCC\n")
    (output / "metadata.tsv").write_text(
        "run_accession\tbiosample\tcollection_date\tlatitude\tlongitude\tgeo_loc_name\thost\tisolation_source\n"
        "SRR_DEMO\tSAMN_DEMO\t2026-01-15\t38.9072\t-77.0369\tUSA: Washington, DC\thuman\tsynthetic example\n"
    )
    write_fastq(
        output / "SRR_DEMO_1.fastq.gz",
        [(f"SRR_DEMO.{index}/1", product[:60]) for index in range(1, 5)],
    )
    write_fastq(
        output / "SRR_DEMO_2.fastq.gz",
        [(f"SRR_DEMO.{index}/2", reverse_complement(product[-60:])) for index in range(1, 5)],
    )


if __name__ == "__main__":
    main()
