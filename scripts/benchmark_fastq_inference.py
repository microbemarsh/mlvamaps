"""Time real FASTQ workflows against the assembly-concordance fixture.

Run from a checkout with dev dependencies, minimap2 and sassy on PATH:
python scripts/benchmark_fastq_inference.py --repetitions 3 --copies 25

Only FASTQ calling is timed; database construction and assembly truth are not.
The fixture asserts the final fingerprint, including novel/fractional alleles.
"""

import argparse
import contextlib
import io
import json
from pathlib import Path
import runpy
import statistics
import sys
import tempfile
import time

import pytest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--copies", type=int, default=25)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--layouts", nargs="+", choices=("paired", "shotgun", "long"),
                        default=["paired", "shotgun", "long"])
    parser.add_argument("--sample-mode", choices=("isolate", "metagenome"), default="isolate")
    args = parser.parse_args()
    if args.repetitions < 1 or args.copies < 1:
        parser.error("repetitions and copies must be positive")
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    fixture = runpy.run_path(str(root / "tests/test_fastq_assembly_concordance.py"))
    test = fixture["test_short_read_fingerprint_matches_assembly"]
    original = test.__globals__["main"]
    timings = {}
    for layout in args.layouts:
        elapsed = []

        def timed_call(argv):
            argv[argv.index("--sample-mode") + 1] = args.sample_mode
            first = argv[argv.index("--fq1") + 1] if "--fq1" in argv else argv[argv.index("-i") + 1]
            paths = [first]
            if "--fq2" in argv:
                paths.append(argv[argv.index("--fq2") + 1])
            for name in paths:
                path = Path(name)
                lines = path.read_text().splitlines()
                with path.open("w") as handle:
                    for copy in range(args.copies):
                        for start in range(0, len(lines), 4):
                            header, sequence, plus, quality = lines[start:start + 4]
                            handle.write(f"@copy{copy}_{header[1:]}\n{sequence}\n{plus}\n{quality}\n")
            started = time.perf_counter()
            result = original(argv)
            elapsed.append(time.perf_counter() - started)
            return result

        test.__globals__["main"] = timed_call
        for _ in range(args.repetitions):
            with tempfile.TemporaryDirectory(prefix="mlvamaps-benchmark-") as directory:
                with pytest.MonkeyPatch.context() as patch, contextlib.redirect_stdout(io.StringIO()):
                    test(Path(directory), layout, .25, patch)
        timings[layout] = {"seconds": elapsed, "median_seconds": statistics.median(elapsed)}
    report = {"copies": args.copies, "repetitions": args.repetitions,
              "sample_mode": args.sample_mode,
              "threads": 2, "fingerprints_match_assembly": True, "workflows": timings}
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
