from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from types import SimpleNamespace

from mlva_seer import sequence
from mlva_seer.assembly_call import (
    build_minimap2_command,
    read_alignment_depth,
    read_minimap2_depth,
    run_assembly_call,
)
from mlva_seer.cli import build_parser, main
from mlva_seer.clustering import build_savont_command, cluster_vntr_asvs
from mlva_seer.in_silico_pcr import build_amplirust_command, expected_amplicon_bounds, write_amplirust_primers
from mlva_seer.io import read_loci
from mlva_seer.models import Locus, ReadRecord, RepeatFeature
from mlva_seer.ml_classifier import predict_read_alleles
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


def make_repeat_feature(read_id, sequence, repeat_count=4):
    return RepeatFeature(
        read_id=read_id,
        locus_id="VNTR",
        repeat_region_start=0,
        repeat_region_end=len(sequence),
        repeat_region_length_bp=len(sequence),
        repeat_motif="ATG",
        raw_repeat_count_estimate=float(repeat_count),
        nearest_integer_repeat_count=repeat_count,
        flank_quality_score=1.0,
        repeat_pattern="ATG-ATG-ATG-ATG",
        repeat_sequence=sequence,
        mean_qscore=30.0,
        mismatch_count_in_repeat_region=0,
        motif_kmer_count=repeat_count,
        left_primer_score=1.0,
        right_primer_score=1.0,
        left_flank_score=0.0,
        right_flank_score=0.0,
        amplicon_sequence=sequence,
        amplicon_quality="I" * len(sequence),
    )


def write_fake_savont(tmp_path):
    executable = tmp_path / "savont"
    executable.write_text(
        """#!/usr/bin/env python3
import collections, pathlib, sys
if '--version' in sys.argv:
    print('savont 0.6.1')
    raise SystemExit(0)
args = sys.argv[2:]
flag = args.index('--output-dir')
inputs = [pathlib.Path(value) for value in args[:flag]]
out = pathlib.Path(args[flag + 1]); out.mkdir(parents=True, exist_ok=True)
samples = [path.stem for path in inputs]
clusters = []
for path in inputs:
    lines = path.read_text().splitlines()
    reads = [(lines[i][1:].split()[0], lines[i + 1]) for i in range(0, len(lines), 4)]
    consensus = collections.Counter(sequence for _, sequence in reads).most_common(1)[0][0]
    clusters.append((reads, consensus))
with (out / 'final_asvs.fasta').open('w') as handle:
    for idx, (reads, consensus) in enumerate(clusters):
        handle.write(f'>final_consensus_{idx}_depth_{len(reads)}\\n{consensus}\\n')
with (out / 'feature-table.tsv').open('w') as handle:
    handle.write('#OTU ID\\t' + '\\t'.join(samples) + '\\n')
    for idx, (reads, _consensus) in enumerate(clusters):
        depths = ['0'] * len(samples); depths[idx] = str(len(reads))
        handle.write(f'final_consensus_{idx}_depth_{len(reads)}\\t' + '\\t'.join(depths) + '\\n')
with (out / 'final_clusters.tsv').open('w') as handle:
    for idx, (reads, _consensus) in enumerate(clusters):
        handle.write(f'final_cluster_{idx}\\tsize_{len(reads)}\\trepresentative_0\\tmembers\\n')
        for read_id, _sequence in reads:
            handle.write(f'{read_id} 100\\n')
"""
    )
    executable.chmod(0o755)
    return executable


def write_fake_amplirust(tmp_path):
    executable = tmp_path / "amplirust"
    executable.write_text(
        """#!/usr/bin/env python3
import csv, pathlib, sys
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
def rc(seq): return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
records = []
name = None
parts = []
for line in pathlib.Path(value('--input')).read_text().splitlines():
    if line.startswith('>'):
        if name is not None: records.append((name, ''.join(parts)))
        name = line[1:].split()[0]; parts = []
    else: parts.append(line.strip())
if name is not None: records.append((name, ''.join(parts)))
with pathlib.Path(value('--primers')).open() as handle:
    primers = list(csv.DictReader(handle))
products = []
for reference, original in records:
    case = 0
    for primer in primers:
        for strand, sequence in [('+', original), ('-', rc(original))]:
            fwd = sequence.find(primer['forward'])
            rev_seq = rc(primer['reverse'])
            rev = sequence.find(rev_seq, fwd + len(primer['forward'])) if fwd >= 0 else -1
            if fwd < 0 or rev < 0: continue
            end = rev + len(rev_seq); full_len = end - fwd
            min_len = int(value('--min-len')); max_len = int(value('--max-len'))
            if not min_len <= full_len <= max_len: continue
            case += 1
            start0, end0 = (fwd, end) if strand == '+' else (len(original) - end, len(original) - fwd)
            suffix = '_rc' if strand == '-' else ''
            amplicon_id = f"{reference}:{primer['name']}{suffix}:{case}"
            product = sequence[fwd:end]
            products.append((amplicon_id, reference, primer['name'], product, full_len, fwd, rev, strand, start0, end0, len(primer['forward']), len(rev_seq)))
            break
with pathlib.Path(value('--output')).open('w') as fasta:
    for row in products:
        fasta.write(f'>{row[0]}\\tpos={row[8]}-{row[9]}\\tstrand={row[7]}\\tlen={row[4]}\\n{row[3]}\\n')
header = 'amplicon_id\treference_id\tsource_file\tprimer_name\tproduct_len\tfull_len\tfwd_start\tfwd_end\tfwd_mismatches\tfwd_identity\tfwd_cigar\trev_start\trev_end\trev_mismatches\trev_identity\trev_cigar\tstrand\tis_circular_wrap\tproduct_seq\\n'
with pathlib.Path(value('--tsv')).open('w') as stats:
    stats.write(header)
    for row in products:
        stats.write(f'{row[0]}\t{row[1]}\tfake.fa\t{row[2]}\t{len(row[3])}\t{row[4]}\t{row[5]}\t{row[5] + row[10]}\t0\t1.0\t{row[10]}=\t{row[6]}\t{row[6] + row[11]}\t0\t1.0\t{row[11]}=\t{row[7]}\tfalse\t{row[3]}\\n')
"""
    )
    executable.chmod(0o755)
    return executable


def test_mlva_clustering_uses_savont_and_retains_indels(tmp_path):
    reference = "ATG" * 4
    insertion = reference[:4] + "A" + reference[4:]
    deletion = reference[:4] + reference[5:]
    features = [
        make_repeat_feature("ref1", reference),
        make_repeat_feature("ref2", reference),
        make_repeat_feature("ref3", reference),
        make_repeat_feature("insertion", insertion),
        make_repeat_feature("deletion", deletion),
    ]
    savont = write_fake_savont(tmp_path)
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=5)
    rows, fasta, memberships = cluster_vntr_asvs(
        features, [locus], tmp_path / "savont-out", threads=4, executable=str(savont)
    )

    dominant = rows[0]
    assert dominant["support_reads"] == 5
    assert dominant["unique_sequences"] == 3
    assert dominant["consensus_sequence"] == reference
    assert dominant["reads_with_indels"] == 2
    assert dominant["total_insertions"] == 1
    assert dominant["total_deletions"] == 1
    assert len(rows) == 1
    assert fasta[0] == ("VNTR_ASV1", reference)

    insertion_membership = next(row for row in memberships if row["read_id"] == "insertion")
    assert insertion_membership["variant_id"] == "VNTR_ASV1"
    assert insertion_membership["insertions_vs_consensus"] == 1
    assert insertion_membership["deletions_vs_consensus"] == 0
    assert "-" in insertion_membership["aligned_consensus_sequence"]

    deletion_membership = next(row for row in memberships if row["read_id"] == "deletion")
    assert deletion_membership["insertions_vs_consensus"] == 0
    assert deletion_membership["deletions_vs_consensus"] == 1
    assert "-" in deletion_membership["aligned_repeat_sequence"]

    predictions = predict_read_alleles(features, [locus], memberships)
    by_read = {prediction.read_id: prediction for prediction in predictions}
    assert by_read["insertion"].insertions_vs_consensus == 1
    assert by_read["deletion"].deletions_vs_consensus == 1
    assert by_read["insertion"].evidence_weight < by_read["ref1"].evidence_weight
    assert by_read["deletion"].evidence_weight < by_read["ref1"].evidence_weight


def test_savont_receives_all_loci_and_full_thread_count(tmp_path):
    inputs = [tmp_path / "locus_0000.fastq", tmp_path / "locus_0001.fastq"]
    command = build_savont_command(
        inputs,
        tmp_path / "savont",
        threads=32,
        min_read_length=60,
        max_read_length=400,
        min_cluster_size=2,
    )
    assert command[:4] == ["savont", "asv", str(inputs[0]), str(inputs[1])]
    assert command[command.index("--threads") + 1] == "32"
    assert "--pooled-samples" in command
    assert "--single-strand" in command


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
    savont = write_fake_savont(tmp_path)
    amplirust = write_fake_amplirust(tmp_path)
    result = run_call(
        reads_path=str(sim["reads"]),
        loci_path=str(loci),
        profiles_path=str(profiles),
        outdir=str(tmp_path / "results"),
        sample_id="SIM1",
        min_read_length=20,
        min_depth=5,
        savont_bin=str(savont),
        amplirust_bin=str(amplirust),
    )

    calls = {row["locus_id"]: row for row in read_tsv(result["allele_calls"])}
    assert calls["VNTR_01"]["called_repeat_count"] == "5"
    assert calls["VNTR_02"]["called_repeat_count"] == "4"
    assert calls["VNTR_01"]["call_status"] == "PASS"
    assert calls["VNTR_01"]["effective_read_depth"] == "25.0"
    asv_rows = read_tsv(result["asv_table"])
    assert asv_rows[0]["support_reads"] == "25"
    assert "total_insertions" in asv_rows[0]
    memberships = read_tsv(result["asv_memberships"])
    assert len(memberships) == 50
    assert {row["sample_id"] for row in memberships} == {"SIM1"}
    assert all("aligned_repeat_sequence" in row for row in memberships)
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
    assert command[command.index("--max-n-fraction") + 1] == "0.0"


def test_assembly_call_from_primer_products(tmp_path):
    amplirust = write_fake_amplirust(tmp_path)
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
        amplirust_bin=str(amplirust),
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
    amplirust = write_fake_amplirust(tmp_path)
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
        amplirust_bin=str(amplirust),
    )
    report = result["report"].read_text()
    assert "uniform default intensity" in report
    assert 'opacity="0.740"' in report


def test_easy_cli_accepts_primer_and_fastq_positionals(tmp_path):
    loci, profiles = write_panel(tmp_path)
    savont = write_fake_savont(tmp_path)
    amplirust = write_fake_amplirust(tmp_path)
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
            "--savont-bin",
            str(savont),
            "--amplirust-bin",
            str(amplirust),
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
    assert default_call_args.savont_min_cluster_size == 2
    assert default_call_args.savont_bin == "savont"
    assert default_call_args.amplirust_bin == "amplirust"

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
