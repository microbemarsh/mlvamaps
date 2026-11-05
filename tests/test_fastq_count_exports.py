"""Unresolved FASTQ calls retain diagnostic predictions and observed depth."""
import csv
import math
import shutil

import pytest

from mlvamaps.io import read_loci, write_tsv
from mlvamaps.pipeline import run_call
from mlvamaps.sequence import revcomp
from mlvamaps.short_reads import run_short_read_call
from test_fastq_assembly_concordance import _sample


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native FASTQ tools required")
@pytest.mark.parametrize("technology", ["short", "long"])
@pytest.mark.parametrize("mode", ["isolate", "metagenome"])
@pytest.mark.parametrize("evidence", ["low_depth", "conflicting", "empty"])
def test_rejected_counts_keep_predictions_and_observed_depth_in_every_export(tmp_path, technology, mode, evidence):
    records, panel, _assembly, _expected = _sample(tmp_path)
    loci = read_loci(panel)[:2]
    write_tsv([vars(locus) for locus in loci], panel, list(vars(loci[0])))
    database = tmp_path / "database"
    database.mkdir()
    for locus in loci:
        reference = (locus.forward_primer + locus.left_flank_sequence + "ATGC" * 4
                     + locus.right_flank_sequence + revcomp(locus.reverse_primer))
        (database / f"{locus.locus_id}.fasta").write_text(f">reference\n{reference}\n")
    product = records[0][1]
    sequences = [] if evidence == "empty" else [product, product]
    if evidence == "conflicting":
        sequences[1] = product.replace(loci[0].right_flank_sequence, "ATGC" * 4 + loci[0].right_flank_sequence, 1)
    reads = tmp_path / "reads.fastq"
    reads.write_text("".join(f"@m{index}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                             for index, sequence in enumerate(sequences)))
    output = tmp_path / "output"
    kwargs = dict(loci_path=str(panel), database_path=str(database), outdir=str(output),
                  sample_id="sample", sample_mode=mode, min_depth=20, threads=2, show_progress=False)
    if technology == "short":
        run_short_read_call(reads1_path=str(reads), reads2_path=None, **kwargs)
    else:
        run_call(reads_path=str(reads), **kwargs)

    def rows(name):
        with (output / name).open() as handle:
            return list(csv.DictReader(handle, delimiter="," if name.endswith(".csv") else "\t"))

    common = {row["locus"]: row for row in rows("common_locus_calls.tsv")}
    assert set(common) == {locus.locus_id for locus in loci}
    assert all(math.isfinite(float(row["predicted_repeat_count"])) and float(row["predicted_repeat_count"]) >= 0 for row in common.values())
    assert all(row["repeat_count"] == "" and row["repeat_count_raw"] == ""
               and row["call_accepted"] == "False" and row["call_rejection_reason"]
               for row in common.values())
    assert int(common[loci[0].locus_id]["presence_support"]) == len(sequences)
    assert common[loci[1].locus_id]["prior_only"] == "True"
    for filename in ("calls.tsv", "locus_repeat_counts.tsv", "genotype_evidence.tsv",
                     "mlva_fingerprint_probabilistic.tsv", "myoga_loci.csv"):
        exported = rows(filename)
        assert len(exported) == len(loci)
        for row in exported:
            locus = row.get("locus") or row["locus_id"]
            assert row["repeat_count"] == ""
            assert float(row["predicted_repeat_count"]) == float(common[locus]["predicted_repeat_count"])
            assert row["predicted_repeat_count_raw"] == common[locus]["predicted_repeat_count_raw"]
            assert row["call_accepted"] == "False" and row["call_rejection_reason"]
            if "read_depth" in row:
                assert int(row["read_depth"]) == int(common[locus]["presence_support"])
    fingerprint = rows("mlva_fingerprint.tsv")[0]
    assert all(fingerprint[name] == "" for name in common)
    if technology == "long":
        for row in rows("allele_calls.tsv"):
            assert row["called_repeat_count"] == ""
            assert float(row["predicted_repeat_count"]) == float(common[row["locus_id"]]["predicted_repeat_count"])
            assert int(row["read_depth"]) == int(common[row["locus_id"]]["presence_support"])


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native FASTQ tools required")
@pytest.mark.parametrize("technology", ["short", "long"])
@pytest.mark.parametrize("status", ["estimated", "ambiguous", "mixed"])
def test_profile_tree_excludes_rejected_predictions_regardless_status(tmp_path, monkeypatch, technology, status):
    from mlvamaps import unified_fastq
    original = unified_fastq.run_unified_fastq_inference

    def with_low_confidence_status(**kwargs):
        result = original(**kwargs)
        result[0][0]["status"] = status
        return result

    monkeypatch.setattr(unified_fastq, "run_unified_fastq_inference", with_low_confidence_status)
    test_rejected_counts_keep_predictions_and_observed_depth_in_every_export(
        tmp_path, technology, "isolate", "low_depth")
    with (tmp_path / "output/classification/mlva_profile_tree_status.tsv").open() as handle:
        tree_status = next(csv.DictReader(handle, delimiter="\t"))
    # Neither the low-depth prediction nor the prior becomes a distance input.
    assert tree_status["loci"] == ""
