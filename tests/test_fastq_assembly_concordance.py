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
from mlvamaps.io import read_loci
from mlvamaps.candidate_contexts import generate_candidate_contexts, write_candidate_contexts
from mlvamaps.minimap_mapping import build_competitive_indexes
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
@pytest.mark.parametrize("layout", ["single", "paired", "overlap", "shotgun", "long", "long_unmapped"])
@pytest.mark.parametrize("tolerance", [0.25, 0.4, 0])
def test_short_read_fingerprint_matches_assembly(tmp_path, layout, tolerance, monkeypatch):
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
    # Use a real observed reference allele (4), distinct from most sample
    # alleles. Old expanded candidate assets remain valid database assets but
    # must not define the FASTQ measurement or restrict novel/half alleles.
    database = tmp_path / "database"
    database.mkdir()
    loci = read_loci(panel)
    for locus in loci:
        sequence = (locus.forward_primer + locus.left_flank_sequence + "ATGC" * 4
                    + locus.right_flank_sequence + revcomp(locus.reverse_primer))
        alternate = list(sequence)
        for position in (24, 29, 34, 39):
            alternate[position] = "ACGT"[("ACGT".index(alternate[position]) + 1) % 4]
        (database / f"{locus.locus_id}.fasta").write_text(
            f">R1\n{sequence}\n>R1_copy\n{sequence}\n>R2\n{''.join(alternate)}\n"
        )
    assets = write_candidate_contexts(generate_candidate_contexts(loci, database), database / "competitive_mapping")
    build_competitive_indexes(assets["fasta"], database / "competitive_mapping")
    (database / "manifest.json").write_text("{}")

    def redundant_mapping(*args, **kwargs):
        pytest.fail("Reference classification must reuse the allele caller's observed alignments")

    monkeypatch.setattr("mlvamaps.mapping_classification.map_reads_to_candidates_bam", redundant_mapping)
    if layout == "long_unmapped":
        monkeypatch.setattr("mlvamaps.unified_fastq.map_reads_to_candidates_bam", lambda *args, **kwargs: [])
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
        "--database", str(database),
    ]
    if layout.startswith("long"):
        args[args.index("sr")] = str(reads1)
        del args[args.index("--fq1"):args.index("--fq1") + 2]
    elif layout != "single":
        args.extend(["--fq2", str(reads2)])
    assert main(args) == 0
    assert _fingerprint(output / "mlva_fingerprint.tsv") == expected
    if layout != "long_unmapped":
        with (output / "classification" / "mapping_reference_matches.tsv").open() as handle:
            match = next(csv.DictReader(handle, delimiter="\t"))
        assert match["reference_id"] == "R1"
        assert match["equivalent_references"] == "R1;R1_copy"


@pytest.mark.skipif(
    not shutil.which("minimap2") or not shutil.which("sassy"),
    reason="minimap2 and sassy are needed for workflow estimates",
)
@pytest.mark.parametrize("technology", ["short", "long"])
def test_alignment_estimate_reaches_html_tables_and_fingerprint(tmp_path, monkeypatch, technology):
    from types import SimpleNamespace
    from mlvamaps.pipeline import run_call
    from mlvamaps.short_reads import run_short_read_call

    records, panel, _assembly, _expected = _sample(tmp_path)
    loci = read_loci(panel)[:2]
    write_tsv([vars(locus) for locus in loci], panel, list(vars(loci[0])))
    database = tmp_path / "database"
    database.mkdir()
    for locus in loci:
        reference = (locus.forward_primer + locus.left_flank_sequence + "ATGC" * 4
                     + locus.right_flank_sequence + revcomp(locus.reverse_primer))
        (database / f"{locus.locus_id}.fasta").write_text(f">R1\n{reference}\n")
    # Keep genuine native mappings but simulate failure of direct anchor sizing.
    unresolved = lambda *args, **kwargs: SimpleNamespace(called_allele=None, status="PRESENCE_ONLY", confidence=0)
    monkeypatch.setattr("mlvamaps.short_read_evidence.measure_locus_product", unresolved)
    monkeypatch.setattr("mlvamaps.long_read_evidence.measure_locus_product", unresolved)
    monkeypatch.setattr("mlvamaps.pipeline.assignments_from_pcr", lambda *args, **kwargs: [])
    reads = tmp_path / "reads.fastq"
    product = records[0][1]
    reads.write_text("".join(f"@m{i}\n{product}\n+\n{'I' * len(product)}\n" for i in range(4)))
    kwargs = dict(loci_path=str(panel), database_path=str(database), outdir=str(tmp_path / "out"),
                  sample_id="s", threads=2, show_progress=False, sample_mode="isolate")
    if technology == "short":
        paths = run_short_read_call(reads1_path=str(reads), reads2_path=None, **kwargs)
    else:
        paths = run_call(reads_path=str(reads), **kwargs)
    def rows(path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    calls = rows(paths["calls"])
    assert calls[0]["repeat_count"] == "5" and calls[0]["status"] == "ESTIMATED"
    assert calls[0]["estimation_method"] == "alignment_repeat_span"
    assert calls[0]["estimate_support"] == "4"
    assert calls[1]["repeat_count"] == ""
    assert rows(paths["repeat_counts"])[0]["repeat_count"] == "5"
    assert _fingerprint(paths["fingerprint"])[loci[0].locus_id] == "5"
    report = paths["report"].read_text()
    assert "<strong>≈5</strong>" in report and "repeat-region CIGAR indels" in report
    assert loci[1].locus_id in report and "Unavailable" in report
    evidence = rows(tmp_path / "out" / "molecule_candidate_evidence.tsv")
    assert all(row["inferred_molecule_repeat"] == "" for row in evidence)
