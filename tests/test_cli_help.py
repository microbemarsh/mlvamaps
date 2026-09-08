import pytest

from mlvamaps.cli import build_parser, main


@pytest.mark.parametrize("flags", [["-h"], ["--help"], ["--advanced"], ["-h", "--advanced"], ["--advanced", "-h"]])
def test_call_help(flags, capsys):
    with pytest.raises(SystemExit) as exc:
        main(["call", *flags])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    for option in ("--input", "--panel", "--database", "--output", "--threads", "--advanced"):
        assert option in output
    for option in ("--min-posterior", "--minimap2-bin", "--short-min-mapq"):
        assert (option in output) == ("--advanced" in flags)


def test_hidden_options_still_parse():
    args = build_parser().parse_args(["call", "-i", "sample.fastq", "--min-posterior", "0.9"])
    assert args.min_posterior == 0.9
