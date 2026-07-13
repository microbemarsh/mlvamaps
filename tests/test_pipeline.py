from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from types import SimpleNamespace

from mlva_seer import sequence
from mlva_seer.assembly_call import (
    build_minimap2_command,
    extract_primer_products,
    read_alignment_depth,
    read_minimap2_depth,
    run_assembly_call,
)
from mlva_seer.cli import build_parser, main
from mlva_seer.in_silico_pcr import build_amplirust_command, expected_amplicon_bounds, write_amplirust_primers
from mlva_seer.io import read_loci
from mlva_seer.locus_assignment import assign_reads
from mlva_seer.models import Locus, ReadRecord
from mlva_seer.pipeline import run_call
from mlva_seer.primers import read_primer_pairs
from mlva_seer.simulation import simulate_reads
from scripts.convert_uf_ba_vntrs import convert_profiles


def write_panel(tmp_path):
    loci = tmp_path / "mlva_loci.tsv"
    loci.write_text(
        "\t".join(
            [
                "locus_id",
                "chrom_or_contig",
                "start",
                "end",
                "forward_primer",
                "reverse_primer",
                "left_flank_sequence",
                "right_flank_sequence",
                "repeat_motif",
                "expected_min_repeats",
                "expected_max_repeats",
                "expected_amplicon_min_bp",
                "expected_amplicon_max_bp",
                "pool_id",
            ]
        )
        + "\n"
        + "VNTR_01\tchr1\t1\t100\tACGTTGCAAC\tTGCATGCAAA\tGGTA\tCCAT\tATG\t3\t8\t30\t80\tA\n"
        + "VNTR_02\tchr1\t200\t300\tTTGACCGTAA\tAACCGGTTCA\tCCTA\tTAGG\tGATA\t2\t6\t30\t90\tA\n"
    )
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text("profile_id\tstrain_id\tVNTR_01\tVNTR_02\tmetadata\nP1\tS1\t5\t4\tknown\n")
    return loci, profiles


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_simulate_and_call_pipeline(tmp_path):
    loci, profiles = write_panel(tmp_path)
    sim = simulate_reads(
        loci_path=str(loci),
        profiles_path=str(profiles),
        profile_id="P1",
        outdir=str(tmp_path / "sim"),
        sample_id="SIM1",
        depth=25,
        error_rate=0.0,
    )
    result = run_call(
        reads_path=str(sim["reads"]),
        loci_path=str(loci),
        profiles_path=str(profiles),
        outdir=str(tmp_path / "results"),
        sample_id="SIM1",
        min_read_length=20,
        min_depth=5,
    )

    calls = {row["locus_id"]: row for row in read_tsv(result["allele_calls"])}
    assert calls["VNTR_01"]["called_repeat_count"] == "5"
    assert calls["VNTR_02"]["called_repeat_count"] == "4"
    assert calls["VNTR_01"]["call_status"] == "PASS"
    easy_calls = {row["locus_id"]: row for row in read_tsv(result["calls"])}
    assert easy_calls["VNTR_01"]["present"] == "yes"
    assert easy_calls["VNTR_01"]["repeat_count"] == "5"
    matches = read_tsv(tmp_path / "results" / "profile_matches.tsv")
    assert matches[0]["best_profile_id"] == "P1"
    assert matches[0]["distance"] == "0.0"
    report = result["report"].read_text()
    assert "Generated MLVA agarose gel comparison" in report
    assert "P1" in report
    assert "query-band" in report
    assert "reference-band" in report
    assert "query band intensity = fragment read support" in report
    assert "25 reads" in report


def test_dropout_is_reported(tmp_path):
    loci, _profiles = write_panel(tmp_path)
    reads = tmp_path / "empty.fastq"
    reads.write_text("")
    result = run_call(
        reads_path=str(reads),
        loci_path=str(loci),
        outdir=str(tmp_path / "dropout"),
        sample_id="EMPTY",
        min_read_length=20,
    )
    statuses = {row["locus_id"]: row["call_status"] for row in read_tsv(result["allele_calls"])}
    assert statuses == {"VNTR_01": "LOCUS_DROPOUT", "VNTR_02": "LOCUS_DROPOUT"}
    easy_calls = {row["locus_id"]: row for row in read_tsv(result["calls"])}
    assert easy_calls["VNTR_01"]["present"] == "no"


def test_amplirust_primer_export_and_command(tmp_path):
    loci_path, _profiles = write_panel(tmp_path)
    loci = read_loci(loci_path)
    primers = write_amplirust_primers(loci, tmp_path / "amplirust_primers.csv")
    primer_rows = list(csv.DictReader(primers.open()))
    assert primer_rows[0] == {"name": "VNTR_01", "forward": "ACGTTGCAAC", "reverse": "TGCATGCAAA"}
    assert expected_amplicon_bounds(loci) == (30, 90)

    command = build_amplirust_command(
        input_path="assembly.fasta",
        primers_path=primers,
        output_fasta=tmp_path / "products.fasta",
        stats_tsv=tmp_path / "stats.tsv",
        min_len=30,
        max_len=90,
        max_errors=3,
        circular=True,
    )
    assert command[:7] == [
        "amplirust",
        "--input",
        "assembly.fasta",
        "--primers",
        str(primers),
        "--output",
        str(tmp_path / "products.fasta"),
    ]
    assert "--search-rc" in command
    assert "--circular" in command
    assert command[command.index("--max-errors") + 1] == "3"


def test_assembly_call_from_primer_products(tmp_path):
    primers = tmp_path / "primers.tsv"
    primers.write_text(
        "locus_id\tforward_primer\treverse_primer\trepeat_unit_length_bp\texpected_product_size_bp\tnominal_repeat_units\n"
        "VNTR_01\tACGTTGCAAC\tTGCATGCAAA\t3\t43\t5\n"
        "VNTR_02\tTTGACCGTAA\tAACCGGTTCA\t4\t34\t4\n"
    )
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nNNNNACGTTGCAACGGTA" + "ATG" * 7 + "CCATTTGCATGCAANNNN\n")
    sam = tmp_path / "assembly_reads.sam"
    sam.write_text(
        "@SQ\tSN:contig1\tLN:56\n"
        "read1\t0\tcontig1\t5\t60\t48M\t*\t0\t0\tACGT\tIIII\n"
        "read2\t0\tcontig1\t5\t60\t24M\t*\t0\t0\tACGT\tIIII\n"
        "read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
    )
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text("profile_id\tstrain_id\tVNTR_01\tVNTR_02\tmetadata\nASM_MATCH\tstrain_A\t7\t\tassembly profile\n")

    result = run_assembly_call(
        assembly_path=str(assembly),
        loci_path=None,
        primers_path=str(primers),
        outdir=str(tmp_path / "assembly_results"),
        sample_id="ASM1",
        alignments_path=str(sam),
        profiles_path=str(profiles),
    )

    calls = {row["locus_id"]: row for row in read_tsv(result["calls"])}
    assert calls["VNTR_01"]["present"] == "yes"
    assert calls["VNTR_01"]["repeat_count"] == "7"
    assert calls["VNTR_01"]["product_size_bp"] == "48"
    assert calls["VNTR_01"]["read_depth"] == "2"
    assert calls["VNTR_01"]["mean_coverage"] == "1.5"
    assert calls["VNTR_02"]["present"] == "no"
    report = result["report"].read_text()
    assert "MLVA Seer Assembly Report: ASM1" in report
    assert "VNTR_01" in report
    assert "Assembly Amplicons" in report
    assert "Generated MLVA assembly gel electrophoresis image" in report
    assert "intensity = depth support" in report
    assert "reference-band" in report
    assert "Closest-reference bands are drawn in magenta" in report
    assert "Closest MLVA Profiles" in report
    assert "ASM_MATCH" in report
    matches = read_tsv(result["profile_matches"])
    assert matches[0]["best_profile_id"] == "ASM_MATCH"
    assert matches[0]["distance"] == "0.0"


def test_assembly_report_uses_default_band_intensity_without_depth(tmp_path):
    primers = tmp_path / "primers.tsv"
    primers.write_text(
        "locus_id\tforward_primer\treverse_primer\trepeat_unit_length_bp\texpected_product_size_bp\tnominal_repeat_units\n"
        "VNTR_01\tACGTTGCAAC\tTGCATGCAAA\t3\t43\t5\n"
    )
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nNNNNACGTTGCAACGGTA" + "ATG" * 7 + "CCATTTGCATGCAANNNN\n")
    result = run_assembly_call(
        assembly_path=str(assembly),
        loci_path=None,
        primers_path=str(primers),
        outdir=str(tmp_path / "assembly_no_depth"),
        sample_id="ASM_NO_DEPTH",
    )
    report = result["report"].read_text()
    assert "uniform default intensity" in report
    assert 'opacity="0.740"' in report


def test_easy_cli_accepts_primer_and_fastq_positionals(tmp_path):
    loci, profiles = write_panel(tmp_path)
    sim = simulate_reads(
        loci_path=str(loci),
        profiles_path=str(profiles),
        profile_id="P1",
        outdir=str(tmp_path / "sim"),
        sample_id="SIM1",
        depth=12,
        error_rate=0.0,
    )
    exit_code = main(
        [
            "call",
            str(loci),
            str(sim["reads"]),
            "--outdir",
            str(tmp_path / "cli_results"),
            "--sample-id",
            "SIM1",
            "--min-read-length",
            "20",
            "--min-depth",
            "5",
        ]
    )
    assert exit_code == 0
    calls = {row["locus_id"]: row for row in read_tsv(tmp_path / "cli_results" / "calls.tsv")}
    assert calls["VNTR_01"]["present"] == "yes"
    assert calls["VNTR_02"]["repeat_count"] == "4"


def test_cli_has_conventional_output_and_thread_options():
    parser = build_parser()
    call_args = parser.parse_args(["call", "primers.tsv", "sample.fastq.gz", "-o", "run", "-t", "4"])
    assert call_args.outdir == "run"
    assert call_args.threads == 4

    default_call_args = parser.parse_args(["call", "primers.tsv", "sample.fastq.gz"])
    assert default_call_args.outdir == "results"
    assert default_call_args.threads == 32

    extract_args = parser.parse_args(["extract-amplicons", "--input", "assembly.fasta", "--primers", "p.tsv"])
    assert extract_args.threads == 32


def test_minimap2_depth_parser(tmp_path):
    sam = tmp_path / "support.sam"
    sam.write_text(
        "@SQ\tSN:VNTR_01|contig1|forward|1-39\tLN:39\n"
        "read1\t0\tVNTR_01|contig1|forward|1-39\t1\t60\t39M\t*\t0\t0\tACGT\tIIII\n"
        "read2\t256\tVNTR_01|contig1|forward|1-39\t1\t60\t39M\t*\t0\t0\tACGT\tIIII\n"
        "read3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
    )
    depth = read_minimap2_depth(sam, {"VNTR_01|contig1|forward|1-39": 39})
    assert depth["VNTR_01|contig1|forward|1-39"]["mapped_reads"] == 1
    assert depth["VNTR_01|contig1|forward|1-39"]["mean_coverage"] == 1.0
    command = build_minimap2_command("amplicons.fasta", "reads.fastq.gz", threads=2)
    assert command == ["minimap2", "-a", "-t", "2", "amplicons.fasta", "reads.fastq.gz"]
    preset_command = build_minimap2_command("amplicons.fasta", "reads.fastq.gz", threads=2, preset="sr")
    assert preset_command == ["minimap2", "-a", "-t", "2", "-x", "sr", "amplicons.fasta", "reads.fastq.gz"]


def test_assembly_alignment_depth_from_sam(tmp_path):
    products = [
        {
            "product_id": "VNTR_01|contig1|forward|5-52",
            "locus_id": "VNTR_01",
            "contig": "contig1",
            "contig_start": 5,
            "contig_end": 52,
            "product_size_bp": 48,
        }
    ]
    sam = tmp_path / "assembly_reads.sam"
    sam.write_text(
        "@SQ\tSN:contig1\tLN:56\n"
        "read1\t0\tcontig1\t5\t60\t48M\t*\t0\t0\tACGT\tIIII\n"
        "read2\t0\tcontig1\t20\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
        "read3\t0\tother\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
    )
    depth = read_alignment_depth(sam, products)
    assert depth["VNTR_01|contig1|forward|5-52"]["mapped_reads"] == 2
    assert round(depth["VNTR_01|contig1|forward|5-52"]["mean_coverage"], 3) == 1.208


def test_sassy_is_preferred_for_approximate_matching(monkeypatch):
    class FakeSearcher:
        def __init__(self, alphabet, rc=False):
            self.alphabet = alphabet
            self.rc = rc

        def search(self, pattern, text, k):
            return [
                SimpleNamespace(text_start=7, cost=1),
                SimpleNamespace(text_start=3, cost=0),
            ]

    fake_sassy = SimpleNamespace(Searcher=FakeSearcher)
    monkeypatch.setattr(sequence, "sassy", fake_sassy)
    sequence._clear_sassy_searchers()
    assert sequence.find_best("ACGT", "TTTACGTTT", 1) == (3, 0)
    sequence._clear_sassy_searchers()


def test_sassy_searchers_are_thread_local(monkeypatch):
    class FakeSearcher:
        rendezvous = Barrier(2)

        def __init__(self, alphabet, rc=False):
            self.active = False

        def search(self, pattern, text, k):
            if self.active:
                raise RuntimeError("Already borrowed")
            self.active = True
            try:
                self.rendezvous.wait(timeout=2)
                return [SimpleNamespace(text_start=3, cost=0)]
            finally:
                self.active = False

    monkeypatch.setattr(sequence, "sassy", SimpleNamespace(Searcher=FakeSearcher))
    sequence._clear_sassy_searchers()
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(lambda _index: sequence.find_best("ACGT", "TTTACGTTT", 1), range(2)))
    assert results == [(3, 0), (3, 0)]


def test_sassy_batch_search_owns_requested_threads(monkeypatch):
    class FakeSearcher:
        calls = []

        def __init__(self, alphabet, rc=False):
            self.rc = rc

        def search_many(self, patterns, texts, k, threads, mode):
            self.calls.append((threads, mode, self.rc, len(patterns), len(texts)))
            matches = []
            for pattern_idx, pattern in enumerate(patterns):
                reverse_pattern = sequence.revcomp(pattern.decode()).encode()
                for text_idx, text in enumerate(texts):
                    for strand, query in (("+", pattern), ("-", reverse_pattern)):
                        position = text.find(query)
                        if position >= 0:
                            matches.append(
                                SimpleNamespace(
                                    pattern_idx=pattern_idx,
                                    text_idx=text_idx,
                                    text_start=position,
                                    text_end=position + len(query),
                                    cost=0,
                                    strand=strand,
                                )
                            )
            return matches

    monkeypatch.setattr(sequence, "sassy", SimpleNamespace(Searcher=FakeSearcher))
    sequence._clear_sassy_searchers()
    locus = Locus(
        locus_id="VNTR",
        forward_primer="ACGTAC",
        reverse_primer="AACCGT",
        expected_amplicon_min_bp=10,
        expected_amplicon_max_bp=30,
    )
    amplicon = "ACGTAC" + "GATA" * 2 + sequence.revcomp(locus.reverse_primer)
    reads = [
        ReadRecord("forward", amplicon, "I" * len(amplicon)),
        ReadRecord("reverse", sequence.revcomp(amplicon), "I" * len(amplicon)),
    ]

    assignments = assign_reads(reads, [locus], "SAMPLE", max_primer_mismatches=0, threads=7)

    assert FakeSearcher.calls == [(7, "batch_texts", True, 2, 2)]
    assert [assignment.orientation for assignment in assignments] == ["forward", "reverse"]
    assert all(assignment.assigned_locus == "VNTR" for assignment in assignments)
    assert all(assignment.passes_assignment_qc for assignment in assignments)


def test_cleaned_mlva_seer_primer_tsv_is_ingestible():
    loci = read_primer_pairs("examples/seer_lab_Ba/mlva_seer_primers.example.tsv")
    assert len(loci) == 31
    assert loci[0].locus_id == "vrrA_12bp_314bp_10U"
    assert loci[0].forward_primer == "CACAACTACCACCGATGGCACA"
    assert loci[0].reverse_primer == "GCGCGTTTCGTTTGATTCATAC"
    assert len(loci[0].repeat_motif) == 12
    assert loci[0].expected_amplicon_min_bp == 194
    assert loci[0].expected_amplicon_max_bp == 434


def test_raw_legacy_primer_file_is_ingestible():
    loci = read_primer_pairs("examples/seer_lab_Ba/insilicoMLVAprimers_all.raw.example.csv")
    assert len(loci) == 31
    assert loci[-1].locus_id == "Bavntr35_6bp_115bp_5U"
    assert len(loci[-1].repeat_motif) == 6


def test_uf_ba_profile_converter_maps_short_locus_names(tmp_path):
    source = tmp_path / "uf_ba_vntrs.tsv"
    source.write_text("Access_number\tvrrA\tBams13\nEMPTY\t\t\nUF1\t10\t72\n")
    output = tmp_path / "profiles.tsv"
    rows_written, mapped_loci, unmatched = convert_profiles(
        source,
        primers_path=Path("examples/seer_lab_Ba/mlva_seer_primers.example.tsv"),
        output_path=output,
    )
    assert rows_written == 1
    assert mapped_loci == 2
    assert unmatched == []
    rows = read_tsv(output)
    assert rows[0]["profile_id"] == "UF1"
    assert rows[0]["vrrA_12bp_314bp_10U"] == "10"
    assert rows[0]["Bams13_9bp_814bp_70U"] == "72"
