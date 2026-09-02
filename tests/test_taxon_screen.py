import json
import csv
import shutil
from pathlib import Path

import pytest

from mlvamaps.pipeline import run_call
from mlvamaps.sequence import revcomp
from mlvamaps.taxon_screen import (
    build_deacon_filter_command,
    run_taxon_screen,
)


def test_deacon_command_retains_target_reads_with_native_threads(tmp_path):
    command = build_deacon_filter_command(
        "target.idx",
        "reads.fastq.gz",
        "screened.fastq.gz",
        "summary.json",
        threads=12,
        absolute_threshold=3,
        relative_threshold=0.02,
        executable="/usr/bin/deacon",
    )
    assert command == [
        "/usr/bin/deacon",
        "filter",
        "--abs-threshold",
        "3",
        "--rel-threshold",
        "0.02",
        "--threads",
        "12",
        "--output",
        "screened.fastq.gz",
        "--summary",
        "summary.json",
        "--quiet",
        "target.idx",
        "reads.fastq.gz",
    ]


def test_taxon_screen_keeps_deacon_summary_and_fastq(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@target\nACGT\n+\nIIII\n")
    index = tmp_path / "target.idx"
    index.write_bytes(b"index")
    executable = tmp_path / "deacon"
    executable.write_text(
        """#!/usr/bin/env python3
import gzip, json, pathlib, sys
args = sys.argv
output = pathlib.Path(args[args.index("--output") + 1])
summary = pathlib.Path(args[args.index("--summary") + 1])
input_path = pathlib.Path(args[-1])
with gzip.open(output, "wt") as handle:
    handle.write(input_path.read_text())
summary.write_text(json.dumps({
    "seqs_in": 2, "seqs_out": 1, "seqs_removed": 1,
    "bp_in": 8, "bp_out": 4, "bp_removed": 4
}))
"""
    )
    executable.chmod(0o755)

    output, summary_path, summary = run_taxon_screen(
        reads,
        index,
        tmp_path / "screen",
        threads=4,
        executable=str(executable),
    )

    assert output.name == "taxon_screened_reads.fastq.gz"
    assert summary_path.is_file()
    assert json.loads(summary_path.read_text()) == summary
    assert summary["seqs_out"] == 1


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 unavailable")
def test_pipeline_screens_metagenome_before_mlva_and_reports_counts(tmp_path):
    loci = tmp_path / "loci.tsv"
    loci.write_text(
        "locus_id\tchrom_or_contig\tstart\tend\tforward_primer\t"
        "reverse_primer\tleft_flank_sequence\tright_flank_sequence\t"
        "repeat_motif\texpected_min_repeats\texpected_max_repeats\t"
        "expected_amplicon_min_bp\texpected_amplicon_max_bp\tpool_id\n"
        "VNTR\tchr1\t1\t100\tACGTTGCAAC\tTGCATGCAAA\tGGTA\tCCAT\t"
        "ATG\t3\t8\t30\t80\tA\n"
    )
    product = (
        "ACGTTGCAAC"
        + "GGTA"
        + ("ATG" * 5)
        + "CCAT"
        + revcomp("TGCATGCAAA")
    )
    reads = tmp_path / "metagenome.fastq"
    reads.write_text(
        f"@target\n{product}\n+\n{'I' * len(product)}\n"
        "@background\nACACACACACACACACACACACACACACAC\n+\n"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
    )
    index = tmp_path / "target.idx"
    index.write_bytes(b"index")
    executable = tmp_path / "deacon"
    executable.write_text(
        """#!/usr/bin/env python3
import gzip, json, pathlib, sys
args = sys.argv
output = pathlib.Path(args[args.index("--output") + 1])
summary = pathlib.Path(args[args.index("--summary") + 1])
input_path = pathlib.Path(args[-1])
first = "\\n".join(input_path.read_text().splitlines()[:4]) + "\\n"
with gzip.open(output, "wt") as handle:
    handle.write(first)
summary.write_text(json.dumps({
    "seqs_in": 2, "seqs_out": 1, "seqs_removed": 1,
    "bp_in": 73, "bp_out": 43, "bp_removed": 30
}))
"""
    )
    executable.chmod(0o755)

    result = run_call(
        reads_path=str(reads),
        loci_path=str(loci),
        outdir=str(tmp_path / "result"),
        sample_id="META",
        min_read_length=20,
        min_depth=1,
        locus_mapping=False,
        threads=1,
        taxon_screen_index=str(index),
        deacon_bin=str(executable),
    )

    with (tmp_path / "result" / "qc_summary.tsv").open(newline="") as handle:
        qc = {
            row["metric"]: row["value"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    assert qc["taxon_screen_input_reads"] == "2"
    assert qc["taxon_screen_retained_reads"] == "1"
    assert qc["input_reads"] == "1"
    assert result["taxon_screened_reads"].is_file()
    assert result["taxon_screen_summary"].is_file()
    assert "Target-taxon screen" in result["report"].read_text()
