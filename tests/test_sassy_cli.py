from __future__ import annotations

import subprocess

import pytest

from mlvamaps import sassy_cli
from mlvamaps.in_silico_pcr import _record_product_rows
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