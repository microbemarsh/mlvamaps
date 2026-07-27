from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from types import SimpleNamespace

import pytest

from mlvamaps import sequence
from mlvamaps.bayesian_caller import call_loci
from mlvamaps.assembly_call import (
    _primer_match_display,
    amplirust_rows_to_products,
    assembly_call_rows,
    build_minimap2_command,
    legacy_amplicon_bounds,
    read_alignment_depth,
    read_minimap2_depth,
    run_assembly_call,
    run_minimap2_depth,
)
from mlvamaps.calling import legacy_round_repeat_count
from mlvamaps.cli import build_parser, main
from mlvamaps.clustering import (
    _alignment_metrics,
    build_vsearch_cluster_command,
    build_vsearch_derep_command,
    cluster_vntr_asvs,
)
from mlvamaps.in_silico_pcr import build_amplirust_command, expected_amplicon_bounds, write_amplirust_primers
from mlvamaps.io import read_loci
from mlvamaps.mapping import (
    build_minimap2_map_command,
    parse_minimap2_sam,
)
from mlvamaps.mixture import estimate_variant_mixtures
from mlvamaps.models import Locus, ReadPrediction, ReadRecord, RepeatFeature
from mlvamaps.ml_classifier import predict_read_alleles
from mlvamaps.pipeline import run_call
from mlvamaps.primers import read_primer_pairs
from mlvamaps.simulation import simulate_reads
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


def write_fake_vsearch(tmp_path):
    executable = tmp_path / "vsearch"
    executable.write_text(
        """#!/usr/bin/env python3
import collections, pathlib, sys
if '--version' in sys.argv:
    print('vsearch v2.30.0_linux_x86_64')
    raise SystemExit(0)
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
def uc(kind, query, target='*', cluster=0, length='*', identity='*'):
    return '\\t'.join([kind, str(cluster), str(length), str(identity), '+', '0', '0', '*', query, target]) + '\\n'
if '--fastx_uniques' in args:
    lines = pathlib.Path(value('--fastx_uniques')).read_text().splitlines()
    records = [(lines[i][1:].split()[0], lines[i + 1]) for i in range(0, len(lines), 4)]
    groups = collections.OrderedDict()
    for read_id, sequence in records: groups.setdefault(sequence, []).append(read_id)
    ordered = sorted(groups.items(), key=lambda item: -len(item[1]))
    with pathlib.Path(value('--fastaout')).open('w') as fasta, pathlib.Path(value('--uc')).open('w') as out:
        for cluster, (sequence, read_ids) in enumerate(ordered):
            centroid = read_ids[0]
            fasta.write(f'>{centroid};size={len(read_ids)};\\n{sequence}\\n')
            out.write(uc('S', centroid, cluster=cluster, length=len(sequence)))
            for read_id in read_ids[1:]:
                out.write(uc('H', read_id, centroid, cluster, len(sequence), '100.0'))
    raise SystemExit(0)
if '--cluster_size' in args:
    lines = pathlib.Path(value('--cluster_size')).read_text().splitlines()
    records = [(lines[i][1:].split()[0], lines[i + 1]) for i in range(0, len(lines), 2)]
    centroid, centroid_sequence = records[0]
    with pathlib.Path(value('--centroids')).open('w') as fasta:
        fasta.write(f'>{centroid}\\n{centroid_sequence}\\n')
    with pathlib.Path(value('--uc')).open('w') as out:
        out.write(uc('S', centroid, cluster=0, length=len(centroid_sequence)))
        for query, sequence in records[1:]:
            out.write(uc('H', query, centroid, 0, len(sequence), '99.0'))
    raise SystemExit(0)
raise SystemExit(2)
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
            if product.count('N') / len(product) > float(value('--max-n-fraction')): continue
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


def test_vsearch_clustering_uses_observed_centroid_and_retains_indels(tmp_path):
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
    vsearch = write_fake_vsearch(tmp_path)
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=5)
    rows, fasta, memberships = cluster_vntr_asvs(
        features, [locus], tmp_path / "vsearch-out", threads=4, executable=str(vsearch)
    )

    dominant = rows[0]
    assert dominant["support_reads"] == 5
    assert dominant["unique_sequences"] == 3
    assert dominant["representative_sequence"] == reference
    assert dominant["representative_read_id"] == "ref1"
    assert dominant["reads_with_indels"] == 2
    assert dominant["total_insertions"] == 1
    assert dominant["total_deletions"] == 1
    assert len(rows) == 1
    assert fasta[0] == ("VNTR_ASV1", reference)

    insertion_membership = next(row for row in memberships if row["read_id"] == "insertion")
    assert insertion_membership["variant_id"] == "VNTR_ASV1"
    assert insertion_membership["insertions_vs_representative"] == 1
    assert insertion_membership["deletions_vs_representative"] == 0
    assert "-" in insertion_membership["aligned_representative_sequence"]

    deletion_membership = next(row for row in memberships if row["read_id"] == "deletion")
    assert deletion_membership["insertions_vs_representative"] == 0
    assert deletion_membership["deletions_vs_representative"] == 1
    assert "-" in deletion_membership["aligned_repeat_sequence"]

    predictions = predict_read_alleles(features, [locus], memberships)
    by_read = {prediction.read_id: prediction for prediction in predictions}
    assert by_read["insertion"].insertions_vs_representative == 1
    assert by_read["deletion"].deletions_vs_representative == 1
    assert by_read["insertion"].evidence_weight < by_read["ref1"].evidence_weight
    assert by_read["deletion"].evidence_weight < by_read["ref1"].evidence_weight


def test_parasail_global_alignment_retains_indels_and_substitutions():
    insertion = _alignment_metrics("AACGT", "ACGT")
    deletion = _alignment_metrics("ACGT", "AACGT")
    substitution = _alignment_metrics("ACAT", "ACGT")
    empty = _alignment_metrics("", "ACGT")

    assert insertion["insertions_vs_representative"] == 1
    assert "-" in insertion["aligned_representative_sequence"]
    assert deletion["deletions_vs_representative"] == 1
    assert "-" in deletion["aligned_repeat_sequence"]
    assert substitution["substitutions_vs_representative"] == 1
    assert substitution["edit_distance_to_representative"] == 1
    assert empty["aligned_repeat_sequence"] == "----"
    assert empty["deletions_vs_representative"] == 4


def test_emu_style_mixture_filters_trace_variants_and_refines_status():
    asv_rows = [
        {
            "sample_id": "S1",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV1",
            "repeat_count": 5,
            "support_reads": 1000,
            "representative_sequence": "ACGTACGT",
            "representative_length_bp": 8,
            "total_insertions": 2,
            "total_deletions": 2,
            "total_substitutions": 4,
        },
        {
            "sample_id": "S1",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV2",
            "repeat_count": 6,
            "support_reads": 100,
            "representative_sequence": "ACGTACGTTTTT",
            "representative_length_bp": 12,
            "total_insertions": 1,
            "total_deletions": 1,
            "total_substitutions": 2,
        },
        {
            "sample_id": "S1",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV3",
            "repeat_count": 5,
            "support_reads": 2,
            "representative_sequence": "ACGTACGA",
            "representative_length_bp": 8,
            "total_insertions": 0,
            "total_deletions": 0,
            "total_substitutions": 0,
        },
    ]
    mixture = estimate_variant_mixtures(asv_rows, min_fraction=0.01)
    by_variant = {row["variant_id"]: row for row in mixture}
    assert sum(float(row["estimated_fraction"]) for row in mixture) == 1.0
    assert by_variant["VNTR_ASV1"]["abundance_class"] == "DOMINANT"
    assert by_variant["VNTR_ASV2"]["meaningful"] == "yes"
    assert by_variant["VNTR_ASV3"]["meaningful"] == "no"
    assert by_variant["VNTR_ASV3"]["estimated_fraction"] == 0.0

    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=7)
    predictions = [
        ReadPrediction(
            f"read{index}", "VNTR", 5, 1.0, None, 0.0, "VNTR_ASV1", 0, 0, 0, 1.0
        )
        for index in range(20)
    ]
    call = call_loci(
        predictions,
        [locus],
        asv_rows,
        mixture_rows=mixture,
        sample_mode="isolate",
    )[0]
    assert call["call_status"] == "PASS"
    assert call["num_vntr_asvs"] == 3
    assert call["num_meaningful_variants"] == 2

    balanced = [dict(row) for row in mixture]
    balanced[0]["estimated_fraction"] = 0.7
    balanced[1]["estimated_fraction"] = 0.3
    call = call_loci(predictions, [locus], asv_rows, mixture_rows=balanced)[0]
    assert call["call_status"] == "MULTIPLE_VARIANTS"

    dominant_only = [dict(row, meaningful="no") for row in mixture]
    dominant_only[0]["meaningful"] = "yes"
    call = call_loci(predictions, [locus], asv_rows, mixture_rows=dominant_only)[0]
    assert call["call_status"] == "PASS"
    assert call["num_meaningful_variants"] == 1


def test_native_repeat_motif_statistics_preserve_patterns_and_partials():
    parts, mismatches, motif_kmers = sequence.repeat_motif_statistics(
        "ATGATCAT", "ATG", 3
    )
    assert parts == ["ATG", "ATC", "AT:partial"]
    assert mismatches == 1
    assert motif_kmers == 1

    parts, mismatches, motif_kmers = sequence.repeat_motif_statistics(
        "ACGT", "N", 2
    )
    assert parts == ["AC", "GT"]
    assert mismatches == 4
    assert motif_kmers == 0
    assert sequence.mean_qscore("IIII") == 40.0


def test_vsearch_commands_use_native_dereplication_and_full_thread_count(tmp_path):
    input_path = tmp_path / "locus_0000.fastq"
    derep = build_vsearch_derep_command(
        input_path,
        tmp_path / "uniques.fasta",
        tmp_path / "derep.uc",
    )
    command = build_vsearch_cluster_command(
        tmp_path / "uniques.fasta",
        tmp_path / "centroids.fasta",
        tmp_path / "clusters.uc",
        threads=32,
        min_identity=0.97,
    )
    assert derep[:3] == ["vsearch", "--fastx_uniques", str(input_path)]
    assert "--sizeout" in derep
    assert command[:2] == ["vsearch", "--cluster_size"]
    assert command[command.index("--threads") + 1] == "32"
    assert command[command.index("--id") + 1] == "0.97"
    assert command[command.index("--iddef") + 1] == "1"
    assert command[command.index("--qmask") + 1] == "none"
    assert command[command.index("--wordlength") + 1] == "3"
    assert command[command.index("--minwordmatches") + 1] == "1"
    assert command[command.index("--gapopen") + 1] == "4"
    assert "--sizein" in command


def test_vsearch_work_directory_is_reset_before_clustering(tmp_path):
    reference = "ATG" * 4
    features = [make_repeat_feature("ref1", reference), make_repeat_feature("ref2", reference)]
    vsearch = write_fake_vsearch(tmp_path)
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=5)
    retry_dir = tmp_path / "vsearch-retry"
    retry_dir.mkdir()
    (retry_dir / "stale-partial-output.tsv").write_text("incomplete\n")

    rows, _fasta, memberships = cluster_vntr_asvs(
        features, [locus], retry_dir, threads=2, executable=str(vsearch)
    )

    assert rows[0]["support_reads"] == 2
    assert len(memberships) == 2
    assert not (retry_dir / "stale-partial-output.tsv").exists()


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
    vsearch = write_fake_vsearch(tmp_path)
    amplirust = write_fake_amplirust(tmp_path)
    result = run_call(
        reads_path=str(sim["reads"]),
        loci_path=str(loci),
        profiles_path=str(profiles),
        outdir=str(tmp_path / "results"),
        sample_id="SIM1",
        min_read_length=20,
        min_depth=5,
        vsearch_bin=str(vsearch),
        amplirust_bin=str(amplirust),
        locus_mapping=False,
        fastq_strategy="primer",
    )

    calls = {row["locus_id"]: row for row in read_tsv(result["allele_calls"])}
    assert calls["VNTR_01"]["called_repeat_count"] == "5"
    assert calls["VNTR_02"]["called_repeat_count"] == "4"
    assert calls["VNTR_01"]["call_status"] == "PASS"
    assert calls["VNTR_01"]["effective_read_depth"] == "25.0"
    asv_rows = read_tsv(result["asv_table"])
    assert asv_rows[0]["support_reads"] == "25"
    assert asv_rows[0]["representative_read_id"]
    assert asv_rows[0]["representative_sequence"]
    assert "total_insertions" in asv_rows[0]
    assert result["asv_representatives"].exists()
    mixture_rows = read_tsv(result["mixture_abundance"])
    assert {row["abundance_class"] for row in mixture_rows} == {"DOMINANT"}
    assert all(row["estimated_fraction"] == "1.0" for row in mixture_rows)
    assert not (tmp_path / "results" / "vntr_asv_consensus.fasta").exists()
    memberships = read_tsv(result["asv_memberships"])
    assert len(memberships) == 50
    assert {row["sample_id"] for row in memberships} == {"SIM1"}
    assert all("aligned_repeat_sequence" in row for row in memberships)
    assert all("aligned_representative_sequence" in row for row in memberships)
    easy_calls = {row["locus_id"]: row for row in read_tsv(result["calls"])}
    assert easy_calls["VNTR_01"]["present"] == "yes"
    assert easy_calls["VNTR_01"]["repeat_count"] == "5"
    assert easy_calls["VNTR_01"]["primary_read_depth"] == "25"
    assert easy_calls["VNTR_01"]["num_candidate_variants"] == "0"
    assert easy_calls["VNTR_01"]["num_confirmed_secondary_variants"] == "0"
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
    assert "Variant Mixture Abundance" in report
    assert "EM-estimated variant abundance plot" in report
    assert "Locus Confidence" in report
    assert "<details>" in report


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


def test_legacy_amplirust_primer_export_and_removed_command(tmp_path):
    loci_path, _profiles = write_panel(tmp_path)
    loci = read_loci(loci_path)
    primers = write_amplirust_primers(loci, tmp_path / "amplirust_primers.csv")
    primer_rows = list(csv.DictReader(primers.open()))
    assert primer_rows[0] == {"name": "VNTR_01", "forward": "ACGTTGCAAC", "reverse": "TGCATGCAAA"}
    assert expected_amplicon_bounds(loci) == (30, 90)

    with pytest.raises(RuntimeError, match="replaced by MLVAMaps"):
        build_amplirust_command(
            input_path="assembly.fasta",
            primers_path=primers,
            output_fasta=tmp_path / "products.fasta",
            stats_tsv=tmp_path / "stats.tsv",
            min_len=30,
            max_len=90,
            max_errors=3,
            circular=True,
        )


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
    profiles.write_text("profile_id\tstrain_id\tVNTR_01\tVNTR_02\tmetadata\nASM_MATCH\tstrain_A\t6.5\t\tassembly profile\n")

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
    assert calls["VNTR_01"]["repeat_count"] == "6.5"
    assert calls["VNTR_01"]["product_size_bp"] == "48"
    assert calls["VNTR_01"]["read_depth"] == "2"
    assert calls["VNTR_01"]["mean_coverage"] == "1.5"
    assert calls["VNTR_02"]["present"] == "no"
    report = result["report"].read_text()
    assert "MLVAMaps Assembly Report: ASM1" in report
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
    legacy = list(csv.DictReader(result["legacy_details"].open()))
    assert legacy[0]["size"] == "48"
    assert legacy[0]["allele"] == "6.5"
    assert legacy[0]["position1"] == "4"
    assert legacy[0]["position2"] == "42"
    assert legacy[0]["predicted PCR target"].startswith("ACGTTGCAAC")
    legacy_fingerprint = list(csv.reader(result["legacy_fingerprint"].open()))
    assert legacy_fingerprint == [
        ["key", "Access_number", "VNTR_01", "VNTR_02"],
        ["001", "ASM1", "6.5", ""],
    ]
    distribution = read_tsv(result["allele_distribution"])
    assert distribution[0]["allele"] == "6.5"
    assert distribution[0]["selected"] == "yes"
    assert distribution[0]["inference_method"] == "legacy_minimum_allele"


def test_legacy_rounding_and_mismatch_round_product_selection():
    assert legacy_round_repeat_count(6.24) == 6
    assert legacy_round_repeat_count(6.25) == 6.5
    assert legacy_round_repeat_count(6.75) == 6.5
    assert legacy_round_repeat_count(6.76) == 7
    assert _primer_match_display("ACGT", "AGGT") == "AgGT"
    assert _primer_match_display("ACGT", "AGT") == "A.GT"
    assert _primer_match_display("ACGT", "ACAGT") == "ACaGT"

    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
    )
    products = [
        {
            "locus_id": "VNTR",
            "product_id": "one-error-small-allele",
            "primer_error_round": 1,
            "forward_mismatches": 1,
            "reverse_mismatches": 0,
            "size_derived_repeat_count": "2",
            "product_size_bp": 28,
        },
        {
            "locus_id": "VNTR",
            "product_id": "perfect-large-allele",
            "primer_error_round": 0,
            "forward_mismatches": 0,
            "reverse_mismatches": 0,
            "size_derived_repeat_count": "6.25",
            "product_size_bp": 45,
        },
    ]
    row = assembly_call_rows([locus], products, "S1", algorithm="novel")[0]
    assert row["evidence"] == "perfect-large-allele"
    assert row["repeat_count_raw"] == "6.25"
    assert row["repeat_count"] == 6
    assert row["status"] == "AMBIGUOUS"


def test_amplirust_products_use_mlva_finder_size_formula_for_primer_indels():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="CCCC",
        repeat_unit_length_bp=2,
        expected_product_size_bp=20,
        nominal_repeat_units=5,
    )
    row = {
        "primer_name": "VNTR",
        "reference_id": "contig1 description",
        "is_circular_wrap": "false",
        "full_len": "20",
        "fwd_start": "10",
        "fwd_end": "15",
        "rev_start": "26",
        "rev_end": "29",
        "fwd_mismatches": "1",
        "rev_mismatches": "1",
        "strand": "+",
        "original_start": 10,
        "original_end": 30,
        "product_seq": "A" * 20,
    }

    product = amplirust_rows_to_products(
        [row], [locus], "sample", enforce_locus_bounds=False
    )[0]

    assert product["amplicon_span_size_bp"] == 20
    assert product["product_size_bp"] == 19
    assert product["forward_match_length"] == 5
    assert product["reverse_match_length"] == 3


def test_legacy_selection_uses_last_matching_contig_and_forward_precedence():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
    )
    common = {
        "locus_id": "VNTR",
        "primer_error_round": 0,
        "forward_mismatches": 0,
        "reverse_mismatches": 0,
        "forward_match_length": 4,
        "reverse_match_length": 4,
    }
    products = [
        dict(common, product_id="first-small", contig="first", contig_index=0,
             orientation="forward", product_size_bp=32),
        dict(common, product_id="last-large", contig="last", contig_index=1,
             orientation="forward", product_size_bp=44),
        dict(common, product_id="last-reverse-smaller", contig="last", contig_index=1,
             orientation="reverse", product_size_bp=28),
    ]

    row = assembly_call_rows([locus], products, "S1")[0]

    assert row["evidence"] == "last-large"
    assert row["repeat_count"] == 6


def test_assembly_call_allows_n_bases_inside_legacy_product(tmp_path):
    amplirust = write_fake_amplirust(tmp_path)
    primers = tmp_path / "primers.tsv"
    primers.write_text(
        "locus_id\tforward_primer\treverse_primer\trepeat_unit_length_bp\t"
        "expected_product_size_bp\tnominal_repeat_units\n"
        "VNTR\tAAAACCCC\tGGGGTTTT\t2\t24\t4\n"
    )
    assembly = tmp_path / "assembly.fasta"
    assembly.write_text(">contig1\nAAAACCCCNNNNNNNNAAAACCCCAAAAGGGGTTTT\n")

    result = run_assembly_call(
        assembly_path=str(assembly),
        loci_path=None,
        primers_path=str(primers),
        outdir=str(tmp_path / "results"),
        sample_id="sample",
        amplirust_bin=str(amplirust),
    )

    call = read_tsv(result["calls"])[0]
    assert call["present"] == "yes"


def test_depth_distribution_selects_supported_assembly_allele():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
        expected_min_repeats=2,
        expected_max_repeats=8,
    )
    products = []
    for product_id, size in (("allele5", 40), ("allele6", 44)):
        products.append(
            {
                "locus_id": "VNTR",
                "product_id": product_id,
                "primer_error_round": 0,
                "forward_mismatches": 0,
                "reverse_mismatches": 0,
                "size_derived_repeat_count": str(5 + (size - 40) / 4),
                "product_size_bp": size,
            }
        )
    support = {
        "allele5": {"mapped_reads": 2, "mean_coverage": 2.0},
        "allele6": {"mapped_reads": 18, "mean_coverage": 18.0},
    }
    row = assembly_call_rows([locus], products, "S1", support, algorithm="novel")[0]
    assert row["repeat_count"] == 6
    assert row["evidence"] == "allele6"
    assert row["inference_method"] == "depth_weighted_product_distribution"
    assert row["allele_confidence"] > 0.85


def test_legacy_algorithm_is_default_and_uses_upstream_selection_order():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
    )
    products = [
        {
            "locus_id": "VNTR",
            "product_id": "allele6",
            "primer_error_round": 0,
            "forward_mismatches": 0,
            "reverse_mismatches": 0,
            "size_derived_repeat_count": "6",
            "product_size_bp": 44,
        },
        {
            "locus_id": "VNTR",
            "product_id": "allele5.25",
            "primer_error_round": 0,
            "forward_mismatches": 0,
            "reverse_mismatches": 0,
            "size_derived_repeat_count": "5.25",
            "product_size_bp": 41,
        },
        {
            "locus_id": "VNTR",
            "product_id": "better-size-but-later-error-round",
            "primer_error_round": 1,
            "forward_mismatches": 1,
            "reverse_mismatches": 0,
            "size_derived_repeat_count": "2",
            "product_size_bp": 28,
        },
    ]

    row = assembly_call_rows([locus], products, "S1")[0]

    assert row["evidence"] == "allele5.25"
    assert row["repeat_count_raw"] == "5.25"
    assert row["repeat_count"] == 5.5
    assert row["inference_method"] == "legacy_minimum_allele"


def test_legacy_search_bounds_cover_all_upstream_eligible_alleles():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
        expected_amplicon_min_bp=20,
        expected_amplicon_max_bp=80,
    )

    assert legacy_amplicon_bounds([locus]) == (1, 419)


def test_legacy_wide_output_uses_upstream_short_locus_names(tmp_path):
    from mlvamaps.assembly_call import write_legacy_assembly_outputs

    locus = Locus(
        locus_id="Lp03_96bp_941bp_8U",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=96,
        expected_product_size_bp=941,
        nominal_repeat_units=8,
    )
    paths = write_legacy_assembly_outputs(
        tmp_path,
        [locus],
        [],
        assembly_call_rows([locus], [], "sample"),
        "sample",
    )

    rows = list(csv.reader(paths["legacy_fingerprint"].open()))
    assert rows[0] == ["key", "Access_number", "Lp03"]


def test_fastq_read_distribution_can_call_half_allele():
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=7)
    predictions = [
        ReadPrediction(
            f"read{index}", "VNTR", 5.5, 0.95, 5, 0.05,
            "VNTR_ASV1", 0, 0, 0, 1.0,
            raw_repeat_count_estimate=5.48,
            measurement_sigma=0.12,
        )
        for index in range(20)
    ]
    row = call_loci(predictions, [locus], [], min_depth=10)[0]
    assert row["called_repeat_count"] == 5.5
    assert row["posterior_probability"] > 0.95
    assert row["call_status"] == "PASS"
    assert row["allele_distribution"].startswith("5.5:")


def test_read_predictor_retains_unrounded_measurement():
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=7)
    feature = make_repeat_feature("half", "ATG" * 5 + "A", repeat_count=5.5)
    memberships = [
        {
            "read_id": "half",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV1",
            "insertions_vs_representative": 0,
            "deletions_vs_representative": 0,
            "substitutions_vs_representative": 0,
            "aligned_repeat_sequence": "ATG",
            "aligned_representative_sequence": "ATG",
        }
    ]
    prediction = predict_read_alleles([feature], [locus], memberships)[0]
    assert prediction.predicted_repeat_count == 5.5
    assert prediction.raw_repeat_count_estimate == 5.5
    assert prediction.measurement_sigma is not None


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
    vsearch = write_fake_vsearch(tmp_path)
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
            "--vsearch-bin",
            str(vsearch),
            "--amplirust-bin",
            str(amplirust),
            "--no-locus-mapping",
            "--fastq-strategy",
            "primer",
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
    assert default_call_args.min_cluster_size == 1
    assert default_call_args.min_qscore == 17.0
    assert default_call_args.read_calling_convention == "assembly"
    assert default_call_args.sample_mode == "metagenome"
    assert default_call_args.fastq_strategy == "recruit"
    assert default_call_args.recruitment_preset is None
    assert default_call_args.recruitment_min_identity == 0.9
    assert default_call_args.recruitment_min_aligned_bp == 100
    assert default_call_args.recruitment_min_locus_margin == 10
    assert default_call_args.recruitment_database is None
    assert default_call_args.min_secondary_reads == 2
    assert default_call_args.max_confidence_depth == 25.0
    assert default_call_args.cluster_min_identity == 0.97
    assert default_call_args.min_mixture_fraction == 0.01
    assert default_call_args.vsearch_bin == "vsearch"
    assert default_call_args.amplirust_bin == "amplirust"
    assert default_call_args.minimap2_bin == "minimap2"
    assert default_call_args.no_locus_mapping is False
    assert default_call_args.min_mapping_quality == 0
    assert default_call_args.min_base_quality == 20
    assert default_call_args.min_snp_depth == 3
    assert default_call_args.min_snp_alternate_reads == 2
    assert default_call_args.min_snp_frequency == 0.2
    assert default_call_args.algorithm == "legacy"
    assert default_call_args.max_primer_mismatches == 2
    assert default_call_args.raxml_model == "DNA"

    novel_call_args = parser.parse_args(
        ["call", "primers.tsv", "assembly.fasta", "--algorithm", "novel"]
    )
    assert novel_call_args.algorithm == "novel"

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


def test_assembly_read_support_runs_minimap2(tmp_path, monkeypatch):
    reference = tmp_path / "assembly_amplicons.fasta"
    reference.write_text(">VNTR_01\nACGT\n")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read1\nACGT\n+\nIIII\n")
    sam = tmp_path / "read_support.sam"
    commands = []

    monkeypatch.setattr(
        "mlvamaps.assembly_call.check_minimap2", lambda executable: "/fake/minimap2"
    )

    def fake_run(command, stage, stdout_path=None):
        commands.append((command, stage))
        stdout_path.write_text(
            "@SQ\tSN:VNTR_01\tLN:4\n"
            "read1\t0\tVNTR_01\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        )

    monkeypatch.setattr("mlvamaps.assembly_call.run_minimap2_command", fake_run)
    depth = run_minimap2_depth(
        reference,
        reads,
        sam,
        {"VNTR_01": 4},
        threads=2,
        minimap2_preset="sr",
        executable="custom-minimap2",
    )

    assert commands == [
        (
            [
                "/fake/minimap2",
                "-a",
                "-t",
                "2",
                "-x",
                "sr",
                str(reference),
                str(reads),
            ],
            "assembly read-support mapping",
        )
    ]
    assert depth["VNTR_01"] == {"mapped_reads": 1, "mean_coverage": 1.0}


def test_minimap2_command_and_reference_relative_snp_parser(tmp_path):
    reference = tmp_path / "references.fasta"
    reads = tmp_path / "reads.fastq"
    assert build_minimap2_map_command(reference, reads, 4, preset="sr") == [
        "minimap2",
        "-a",
        "-t",
        "4",
        "-x",
        "sr",
        str(reference),
        str(reads),
    ]

    sam = tmp_path / "locus.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:VNTR_ASV1\tLN:4\n"
        "q1\t0\tVNTR_ASV1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "q2\t0\tVNTR_ASV1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "q3\t0\tVNTR_ASV1\t1\t60\t4M\t*\t0\t0\tACGA\tIIII\n"
        "q4\t0\tVNTR_ASV1\t1\t60\t4M\t*\t0\t0\tACGA\tIIII\n"
    )
    references = {
        "VNTR_ASV1": {
            "reference_name": "VNTR_ASV1",
            "locus_id": "VNTR",
            "reference_variant_id": "VNTR_ASV1",
            "reference_read_id": "q1",
            "sequence": "ACGT",
        }
    }
    queries = {
        query_id: {
            "read_id": query_id,
            "locus_id": "VNTR",
            "expected_reference": "VNTR_ASV1",
        }
        for query_id in ("q1", "q2", "q3", "q4")
    }
    summary, snps = parse_minimap2_sam(
        sam,
        references,
        queries,
        "SAMPLE",
        min_snp_depth=4,
        min_snp_alternate_reads=2,
        min_snp_frequency=0.4,
    )
    assert summary[0]["mapped_reads"] == 4
    assert summary[0]["mean_depth"] == 4.0
    assert summary[0]["coverage_percent"] == 100.0
    assert summary[0]["snp_count"] == 1
    assert snps == [
        {
            "sample_id": "SAMPLE",
            "locus_id": "VNTR",
            "reference_variant_id": "VNTR_ASV1",
            "position": 4,
            "reference_base": "T",
            "alternate_base": "A",
            "depth": 4,
            "alternate_depth": 2,
            "alternate_frequency": 0.5,
            "mean_alternate_base_quality": 40.0,
        }
    ]


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


def test_cleaned_mlvamaps_primer_tsv_is_ingestible():
    loci = read_primer_pairs("examples/seer_lab_Ba/mlvamaps_primers.example.tsv")
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
        primers_path=Path("examples/seer_lab_Ba/mlvamaps_primers.example.tsv"),
        output_path=output,
    )
    assert rows_written == 1
    assert mapped_loci == 2
    assert unmatched == []
    rows = read_tsv(output)
    assert rows[0]["profile_id"] == "UF1"
    assert rows[0]["vrrA_12bp_314bp_10U"] == "10"
    assert rows[0]["Bams13_9bp_814bp_70U"] == "72"
