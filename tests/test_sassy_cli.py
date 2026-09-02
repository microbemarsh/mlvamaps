from __future__ import annotations

import subprocess

import pytest

from mlvamaps import sassy_cli


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