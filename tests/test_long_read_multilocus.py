"""A genomic read can independently support several non-overlapping loci."""
import csv
import random
import shutil

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.mapped_reconstruction import recruit_mapping_fragments, resolve_query_products
from mlvamaps.models import Locus


def _mapped(locus, start, end, *, reverse=False, offset=0):
    return CandidateAlignment(
        "m", "m", None, locus, locus, 5, "reference", 140, 60, 1, .14, 1,
        start, end, 0, 140, f"{start}S140M", "", True, False, reverse,
        query_sequence="A" * (160 if offset else 1000), query_quality="I" * (160 if offset else 1000),
        query_offset=offset, query_length=1000,
    )


def _contexts(*names):
    return [CandidateContext(name, name, "A" * 140, 5, 4, 60, 80, "", "") for name in names]


@pytest.mark.parametrize("reverse,cropped", [(False, False), (True, False), (True, True)])
def test_distinct_long_read_windows_each_retain_original_molecule(reverse, cropped):
    first = _mapped("A", 10, 150)
    second = (_mapped("B", 10, 150, reverse=True, offset=250) if cropped else
              _mapped("B", 260 if reverse else 600, 400 if reverse else 740, reverse=reverse))
    contexts = _contexts("A", "B")
    # The paired-fragment policy continues to reject cross-locus competitors.
    assert not recruit_mapping_fragments([first, second], contexts)
    ambiguous = {}
    pools = recruit_mapping_fragments([first, second], contexts, allow_multi_locus=True,
                                     ambiguous_assignments=ambiguous)
    assert set(pools) == {"A", "B"}
    assert all(set(pool) == {"m"} for pool in pools.values())
    assert not ambiguous


def test_competing_locus_windows_remain_ambiguous_without_erasing_an_independent_locus():
    rows = [_mapped("A", 10, 150), _mapped("B", 10, 150), _mapped("C", 600, 740)]
    ambiguous = {}
    pools = recruit_mapping_fragments(rows, _contexts("A", "B", "C"), allow_multi_locus=True,
                                     ambiguous_assignments=ambiguous)
    assert set(pools) == {"C"}
    assert ambiguous == {"m": ("A", "B")}


def test_reconstruction_ambiguity_is_local_to_the_locus(monkeypatch):
    calls = [{"locus": "C", "repeat_count": 5, "status": "called"}]
    evidence = [CandidateEvidence("C", 5, str(i), direct_product_measurement=True,
        measured_repeat_count=5, raw_repeat_count=5, metadata={"product_sequence": "ACGT"}) for i in range(3)]
    evidence += [CandidateEvidence("A", None, str(i), quality_flags=("ambiguous_locus",)) for i in range(3)]
    monkeypatch.setattr("mlvamaps.mapped_reconstruction.call_query_products", lambda queries, *args:
        {"C": ({"status": "PASS", "repeat_count": 5, "repeat_count_raw": 5}, "ACGT")} if queries else {})
    products, audit, _updated, _changed = resolve_query_products(
        calls, evidence, [Locus("C")], [], "sample", InferenceThresholds(), .25)
    assert products == {"C": ("ACGT", 5)}
    assert audit[0]["reconstruction_molecules"] == 3


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native FASTQ tools required")
@pytest.mark.parametrize("mapped_loci", [("A", "B"), ("A",), ()])
def test_genomic_long_reads_call_both_strands_like_the_matching_assembly(tmp_path, monkeypatch, mapped_loci):
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.io import write_tsv
    from mlvamaps.pipeline import run_call
    from mlvamaps.sequence import revcomp
    from mlvamaps import unified_fastq

    native_mapping = unified_fastq.map_reads_to_candidates_bam
    def selected_mapping(*args, **kwargs):
        return [row for row in native_mapping(*args, **kwargs) if row.locus_id in mapped_loci]
    monkeypatch.setattr(unified_fastq, "map_reads_to_candidates_bam", selected_mapping)

    rng = random.Random(3102)
    loci, products = [], []
    database = tmp_path / "database"
    database.mkdir()
    for name, count in (("A", 5), ("B", 8)):
        forward, left, right, reverse = ["".join(rng.choices("ACGT", k=size)) for size in (20, 40, 40, 20)]
        locus = Locus(name, forward_primer=forward, reverse_primer=reverse,
            left_flank_sequence=left, right_flank_sequence=right, repeat_motif="ATGC",
            repeat_unit_length_bp=4, expected_product_size_bp=136, nominal_repeat_units=4)
        loci.append(locus)
        products.append(forward + left + "ATGC" * count + right + revcomp(reverse))
        reference = forward + left + "ATGC" * 4 + right + revcomp(reverse)
        (database / f"{name}.fasta").write_text(f">reference\n{reference}\n")
    panel = tmp_path / "panel.tsv"
    write_tsv([vars(locus) for locus in loci], panel, list(vars(loci[0])))
    genomic = products[0] + "".join(rng.choices("ACGT", k=200)) + revcomp(products[1])
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(f">contig\n{genomic}\n")
    reads = tmp_path / "reads.fastq"
    reads.write_text("".join(f"@m{i}\n{genomic}\n+\n{'I' * len(genomic)}\n" for i in range(4)))
    expected = run_assembly_call(str(assembly), str(panel), str(tmp_path / "assembly"), sample_id="s", threads=2)
    observed = run_call(str(reads), str(panel), str(tmp_path / "fastq"), sample_id="s",
                       database_path=str(database), sample_mode="isolate", threads=2, show_progress=False)
    assert observed["fingerprint"].read_text() == expected["fingerprint"].read_text()
    with observed["calls"].open() as handle:
        calls = {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert {name: row["repeat_count"] for name, row in calls.items()} == {"A": "5", "B": "8"}
    assert all(row["read_depth"] == "4" for row in calls.values())
    with (tmp_path / "fastq/locus_recruited_reads.tsv").open() as handle:
        recruited = list(csv.DictReader(handle, delimiter="\t"))
    assert {(row["read_id"], row["locus_id"]) for row in recruited} == {
        (f"m{i}", locus) for i in range(4) for locus in ("A", "B")}
    assert len(recruited) == 8
    assert all(row["measured_read_allele"] == calls[row["locus_id"]]["repeat_count"] for row in recruited)
    with (tmp_path / "fastq/read_locus_assignments.tsv").open() as handle:
        assignments = list(csv.DictReader(handle, delimiter="\t"))
    assert len(assignments) == 8
    assert {(row["read_id"], row["assigned_locus"]) for row in assignments} == {
        (f"m{i}", locus) for i in range(4) for locus in ("A", "B")}
    with (tmp_path / "fastq/molecule_candidate_evidence.tsv").open() as handle:
        rescued = {row["locus"] for row in csv.DictReader(handle, delimiter="\t")
                   if row["evidence_tier"] == "primer_recovery"}
    assert {"A", "B"} - set(mapped_loci) <= rescued

