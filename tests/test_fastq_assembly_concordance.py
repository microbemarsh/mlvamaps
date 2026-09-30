"""Compare final fingerprints from identical synthetic FASTQ/FASTA molecules.

Expected alleles were checked against i2bc/MLVA_finder, -c -r {0.25,0.4,0},
on 2026-09-30. Upstream MLVA_finder.py SHA256:
61d922dbfa7006c8e10ed11b19a114df41e4b6f6f07e9676683d6ea8329f5b45.
"""

import csv
import random
import shutil

import pytest

from mlvamaps.assembly_call import run_assembly_call
from mlvamaps.cli import main
from mlvamaps.io import write_tsv
from mlvamaps.sequence import revcomp


def _sample(directory):
    rng = random.Random(31)
    records, panel, expected = [], [], {}
    # Raw quarter alleles exercise MLVA_finder's strict rounding boundaries.
    cases = [
        ("integer", 20, "", "5"),
        ("quarter", 21, "", "5.5"),
        ("half", 22, "", "5.5"),
        ("threequarter", 23, "", "5.5"),
        ("novel", 36, "", "9"),
        ("contraction", 4, "", "1"),
        ("zero", 0, "", "0"),
        ("forwardins", 20, "forward_ins", "5"),
        ("forwarddel", 20, "forward_del", "5"),
        ("reverseins", 20, "reverse_ins", "5"),
        ("reversedel", 20, "reverse_del", "5"),
    ]
    for index, (name, repeat_bp, edit, allele) in enumerate(cases):
        forward, left, right, reverse = [
            "".join(rng.choices("ACGT", k=size)) for size in (20, 25, 25, 20)
        ]
        locus = f"{name}_4bp_106bp_4U"
        panel.append(dict(
            locus_id=locus, forward_primer=forward, reverse_primer=reverse,
            left_flank_sequence=left, right_flank_sequence=right,
            repeat_motif="ATGC", repeat_unit_length_bp=4,
            expected_product_size_bp=106, nominal_repeat_units=4,
            expected_min_repeats=3, expected_max_repeats=6,
        ))
        reverse = revcomp(reverse)
        if edit == "forward_ins":
            forward = forward[:10] + "A" + forward[10:]
        elif edit == "forward_del":
            forward = forward[:10] + forward[11:]
        elif edit == "reverse_ins":
            reverse = reverse[:10] + "A" + reverse[10:]
        elif edit == "reverse_del":
            reverse = reverse[:10] + reverse[11:]
        product = forward + left + ("ATGC" * 10)[:repeat_bp] + right + reverse
        # Exercise both input strands without changing the product convention.
        records.append((locus, revcomp(product) if index % 2 else product))
        expected[locus] = allele
    panel_path = directory / "panel.tsv"
    write_tsv(panel, panel_path, list(panel[0]))
    assembly = directory / "assembly.fasta"
    assembly.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
    return records, panel_path, assembly, expected


def _fingerprint(path):
    with path.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    return {name: value for name, value in row.items() if name != "sample_id"}


@pytest.mark.skipif(
    not shutil.which("minimap2") or not shutil.which("sassy"),
    reason="minimap2 and sassy are needed for FASTQ/assembly concordance",
)
@pytest.mark.parametrize("layout", ["single", "paired", "overlap", "shotgun"])
@pytest.mark.parametrize("tolerance", [0.25, 0.4, 0])
def test_short_read_fingerprint_matches_assembly(tmp_path, layout, tolerance):
    records, panel, assembly, expected = _sample(tmp_path)
    if tolerance == 0.4:
        expected.update(quarter_4bp_106bp_4U="5", threequarter_4bp_106bp_4U="6")
    elif tolerance == 0:
        expected.update(quarter_4bp_106bp_4U="5.25", threequarter_4bp_106bp_4U="5.75")
    assembly_result = run_assembly_call(
        assembly_path=str(assembly), loci_path=str(panel),
        outdir=str(tmp_path / "assembly"), sample_id="sample", threads=2,
        assembly_round_tolerance=tolerance,
    )
    assert _fingerprint(assembly_result["fingerprint"]) == expected
    first, second = [], []
    for locus, product in records:
        if layout == "shotgun":
            rng = random.Random(locus)
            contig = ("".join(rng.choices("ACGT", k=100)) + product
                      + "".join(rng.choices("ACGT", k=100)))
            fragments = [contig[start:start + 180] for start in range(0, len(contig) - 179, 5)]
            length = 100
        else:
            fragments = [product] * 4
            length = len(product) - 21 if layout == "overlap" else len(product)
        for molecule, fragment in enumerate(fragments):
            for mate, sequence, output in (
                (1, fragment[:length], first),
                (2, revcomp(fragment[-length:]), second),
            ):
                output.append(f"@{locus}_{molecule}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    reads1, reads2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    reads1.write_text("".join(first))
    reads2.write_text("".join(second))
    output = tmp_path / "fastq"
    args = [
        "call", "-p", str(panel), "-i", "sr", "--fq1", str(reads1),
        "-o", str(output), "--sample-id", "sample", "--sample-mode", "isolate",
        "-t", "2", "--quiet", "--assembly-round-tolerance", str(tolerance),
    ]
    if layout != "single":
        args.extend(["--fq2", str(reads2)])
    assert main(args) == 0
    assert _fingerprint(output / "mlva_fingerprint.tsv") == expected
