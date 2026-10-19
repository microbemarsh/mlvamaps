from __future__ import annotations

import subprocess
import shutil
from pathlib import Path

import pytest

from mlvamaps import sassy_cli
from mlvamaps.in_silico_pcr import _record_batch_product_rows, _record_product_rows
from mlvamaps.models import Locus


def test_searcher_builds_cli_command_and_parses_half_open_matches(monkeypatch):
    observed = {}

    def fake_run(command, capture_output, text, check):
        observed["command"] = command
        return subprocess.CompletedProcess(
            command,
            0,
            "pat_id\ttext_id\tcost\tstrand\tstart\tend\tmatch_region\tcigar\n"
            "pattern\ttext\t1\t+\t2\t11\tACGTTACGA\t3=1D5=\n",
            "",
        )

    monkeypatch.setattr(sassy_cli, "_resolve_executable", lambda: "/conda/bin/sassy")
    monkeypatch.setattr(sassy_cli.subprocess, "run", fake_run)
    matches = sassy_cli.Searcher("ascii", rc=False).search(b"ACGTACGA", b"GGACGTTACGATT", 3)

    assert observed["command"][0:2] == ["/conda/bin/sassy", "search"]
    assert "--no-rc" in observed["command"]
    assert observed["command"][observed["command"].index("--alphabet") + 1] == "iupac"
    assert matches == [sassy_cli.Match(2, 11, 1, "+", "3=1D5=", 0, 8)]


def test_missing_sassy_reports_conda_installation(monkeypatch):
    monkeypatch.delenv("SASSY_BIN", raising=False)
    monkeypatch.setattr(sassy_cli.shutil, "which", lambda _name: None)
    with pytest.raises(RuntimeError, match="conda install -c bioconda sassy"):
        sassy_cli._resolve_executable()


def test_searcher_batches_patterns_and_reuses_results_for_lower_thresholds(monkeypatch):
    calls = []

    def fake_run(command, capture_output, text, check):
        calls.append(command)
        return subprocess.CompletedProcess(
            command,
            0,
            "pat_id\ttext_id\tcost\tstrand\tstart\tend\tmatch_region\tcigar\n"
            "p0\ttext\t0\t+\t2\t6\tACGT\t4=\n"
            "p0\ttext\t2\t+\t8\t12\tACGA\t3=1X\n"
            "p1\ttext\t1\t+\t4\t8\tTTAA\t3=1X\n",
            "",
        )

    monkeypatch.setattr(sassy_cli, "_resolve_executable", lambda: "/conda/bin/sassy")
    monkeypatch.setattr(sassy_cli.subprocess, "run", fake_run)
    searcher = sassy_cli.Searcher("ascii", rc=False)
    text = b"GGACGTTTAAGGACGA"
    searcher.prime([b"ACGT", b"TTAA"], text, 2)

    assert len(calls) == 1
    assert len(searcher.search(b"ACGT", text, 0)) == 1
    assert len(searcher.search(b"ACGT", text, 1)) == 1
    assert len(searcher.search(b"ACGT", text, 2)) == 2
    assert [match.cost for match in searcher.search(b"TTAA", text, 1)] == [1]
    assert len(calls) == 1


def test_record_extraction_primes_both_orientations_once(monkeypatch):
    instances = []

    class FakeSearcher:
        def __init__(self):
            self.prime_calls = []
            instances.append(self)

        def prime(self, patterns, text, k):
            self.prime_calls.append((set(patterns), text, k))

        def search(self, pattern, text, k):
            return []

    monkeypatch.setattr("mlvamaps.in_silico_pcr._new_searcher", FakeSearcher)
    locus = Locus(locus_id="L1", forward_primer="AAR", reverse_primer="CCC")
    task = ("assembly.fa", "contig", "AAACCC", [locus], 2, 1, 100, True, False, 0.0)

    assert _record_product_rows(task) == []
    assert len(instances[0].prime_calls) == 2
    assert instances[0].prime_calls[0][0] == {b"AAA", b"AAG", b"GGG"}
    assert all(call[2] == 2 for call in instances[0].prime_calls)


def test_batched_search_keeps_targets_and_cached_thresholds_separate(monkeypatch):
    calls = []

    def fake_run(command, **kwargs):
        targets = Path(command[-1]).read_text()
        calls.append(targets)
        rows = ("p0\ttext\t2\t+\t2\t6\tACGT\t2=2X\n"
                "p0\ttext1\t1\t+\t4\t8\tACGT\t3=1X\n")
        return subprocess.CompletedProcess(command, 0,
            "pat_id\ttext_id\tcost\tstrand\tstart\tend\tmatch_region\tcigar\n" + rows, "")

    monkeypatch.setattr(sassy_cli, "_resolve_executable", lambda: "/conda/bin/sassy")
    monkeypatch.setattr(sassy_cli.subprocess, "run", fake_run)
    searcher = sassy_cli.Searcher(rc=False)
    first, second = b"GGACGTAA", b"TTTTACGT"
    searcher.prime_many([b"ACGT"], [first, second, first], 2)
    assert len(calls) == 1
    assert calls[0].count(">") == 2
    assert searcher.search(b"ACGT", first, 1) == []
    assert searcher.search(b"ACGT", first, 2)[0].text_start == 2
    assert searcher.search(b"ACGT", second, 1)[0].text_start == 4
    searcher.prime_many([b"ACGT"], [second, first], 1)
    assert len(calls) == 1


@pytest.mark.skipif(not shutil.which("sassy"), reason="Sassy is needed for batch equivalence")
def test_native_batch_matches_individual_pcr_without_merging_records(monkeypatch):
    from mlvamaps.sequence import revcomp

    locus = Locus("L1", forward_primer="ACGTTGCAACGTCAGTACGA",
                  reverse_primer="TGACCGTACAGTCAGGATCA")
    product = locus.forward_primer + "GATA" * 5 + revcomp(locus.reverse_primer)
    reads = [product, product, revcomp(product), product[:8] + "C" + product[8:],
             product[:8] + product[9:], "T" * 80]
    tasks = [("reads.fa", f"read{i}", sequence, [locus], 2, 1, 1000, True, False, 0.0)
             for i, sequence in enumerate(reads)]
    expected = [row for task in tasks for row in _record_product_rows(task)]
    assert expected
    real_run, calls = sassy_cli.subprocess.run, []

    def counted_run(command, **kwargs):
        calls.append(command)
        return real_run(command, **kwargs)

    monkeypatch.setattr(sassy_cli.subprocess, "run", counted_run)
    assert _record_batch_product_rows(tasks) == expected
    assert len(calls) == 1
    assert {row["reference_id"] for row in expected} >= {"read0", "read1", "read2"}
