import gzip
from pathlib import Path
import threading

import pytest

from mlvamaps.io import read_fastq, read_fastq_pairs
from mlvamaps.native_io import short_read_io_plan


def test_decoder_budget_and_optional_dependency_fallback(tmp_path, monkeypatch):
    first, second = tmp_path/'r1.gz', tmp_path/'r2.gz'
    for path in (first, second):
        with path.open('wb') as handle:
            handle.truncate(9 * 1024 * 1024)
    monkeypatch.setattr('mlvamaps.native_io.find_spec', lambda _: object())
    plan = short_read_io_plan((first, second), 32)
    assert plan['decoder_threads'] == [2, 2]
    assert plan['recruitment_workers'] == 25
    assert sum(plan['decoder_threads']) + plan['feeder_threads'] + plan['recruitment_workers'] + 1 == 32
    for threads in (1, 2, 4, 7, 8, 12):
        plan = short_read_io_plan((first, second), threads)
        assert plan['decoder_threads'] == [0, 0]
        assert plan['recruitment_workers'] == threads
    assert short_read_io_plan((first, None), 8)['recruitment_workers'] == 4
    first.write_bytes(gzip.compress(b'hello'))
    assert short_read_io_plan((first, None), 32)['decoder_threads'] == [0, 0]
    monkeypatch.setattr('mlvamaps.native_io.find_spec', lambda _: None)
    assert short_read_io_plan((second, None), 32)['decoder_threads'] == [0, 0]


@pytest.mark.parametrize('content', [
    b'@a/1 comment\nacgtn\n+\nIIIII\n@b\nTT\n+\n!!\n',
    b'@a\nACGT\nACGT\n+\nIIII\nIIII\n',  # HTSlib retains wrapped FASTQ support.
    b'@a\r\nACGT\r\n+\r\nIIII\r\n',
    b'',
])
def test_rapidgzip_retains_native_fastq_parsing(tmp_path, content):
    pytest.importorskip('rapidgzip')
    path = tmp_path/'reads.gz'
    path.write_bytes(gzip.compress(content))
    assert list(read_fastq(path, decompression_threads=2)) == list(read_fastq(path))


def test_rapidgzip_concatenated_members_and_early_close(tmp_path):
    pytest.importorskip('rapidgzip')
    content = b'@a\nACGT\n+\nIIII\n'
    path = tmp_path/'reads.gz'
    path.write_bytes(gzip.compress(content * 100000) + gzip.compress(content))
    assert sum(1 for _ in read_fastq(path, decompression_threads=2)) == 100001
    reader = read_fastq(path, decompression_threads=2)
    next(reader)
    reader.close()
    assert not any(t.name == 'mlvamaps-rapidgzip' for t in threading.enumerate())


def test_rapidgzip_propagates_truncation_and_pair_validation(tmp_path):
    pytest.importorskip('rapidgzip')
    content = b'@a/1\nACGT\n+\nIIII\n'
    first, second = tmp_path/'r1.gz', tmp_path/'r2.gz'
    first.write_bytes(gzip.compress(content)[:-6])
    with pytest.raises(ValueError, match='rapidgzip failed'):
        list(read_fastq(first, decompression_threads=2))
    corrupt = bytearray(gzip.compress(content))
    corrupt[-8] ^= 1  # CRC32 corruption must not silently pass through.
    first.write_bytes(corrupt)
    with pytest.raises(ValueError, match='rapidgzip failed'):
        list(read_fastq(first, decompression_threads=2))
    first.write_bytes(gzip.compress(content))
    second.write_bytes(gzip.compress(content.replace(b'a/1', b'b/2')))
    with pytest.raises(ValueError, match='IDs differ'):
        list(read_fastq_pairs(first, second, decompression_threads=(2, 2)))
    assert not any(t.name == 'mlvamaps-rapidgzip' for t in threading.enumerate())


def test_native_gzip_writer_is_standard_gzip_compatible(tmp_path, monkeypatch):
    from mlvamaps import io
    igzip = pytest.importorskip('isal.igzip')
    monkeypatch.setattr(io, '_gzip_writer', igzip)
    path = tmp_path/'output.gz'
    with io.open_text(path, 'wt') as handle:
        handle.write('sample\tsequence\nα\tACGT\n')
    with io.open_text(path, 'at') as handle:
        handle.write('next\tTT\n')
    assert gzip.decompress(path.read_bytes()).decode() == 'sample\tsequence\nα\tACGT\nnext\tTT\n'


def test_parallel_decoder_feeds_spawned_recruitment_workers(tmp_path, monkeypatch):
    pytest.importorskip('rapidgzip')
    import json
    from dataclasses import asdict
    from mlvamaps.io import write_fastq, write_tsv
    from mlvamaps.models import Locus, ReadRecord
    from mlvamaps.repeat_likelihood import panel_template
    from mlvamaps.sequence import revcomp
    from mlvamaps.short_reads import run_short_read_call
    locus = Locus('X', forward_primer='ACGTCAGTACGATCG', reverse_primer='TCGTAGCACTGATCG',
                  left_flank_sequence='GCTAGCTACGTTACAGTGCA', right_flank_sequence='CGACTGATGCATCGTACGAC',
                  repeat_motif='AATG', repeat_unit_length_bp=4, expected_max_repeats=10)
    sequence = panel_template(locus).sequence(8)
    first, second, panel = tmp_path/'r1.gz', tmp_path/'r2.gz', tmp_path/'panel.tsv'
    write_fastq((ReadRecord(f'm{i}/1', sequence, 'I'*len(sequence)) for i in range(600)), first)
    write_fastq((ReadRecord(f'm{i}/2', revcomp(sequence), 'I'*len(sequence)) for i in range(600)), second)
    write_tsv([asdict(locus)], panel, list(asdict(locus)))
    # Exercise the auto path without making a test fixture eight MiB large.
    monkeypatch.setattr('mlvamaps.native_io.short_read_io_plan',
                        lambda paths, threads: short_read_io_plan(paths, threads, minimum_bytes=0))
    outputs = []
    for threads in (1, 16):
        outputs.append(run_short_read_call(str(first), str(second), str(panel), str(tmp_path/str(threads)),
                                          's', sr_engine='repeat-likelihood', threads=threads, show_progress=False))
    for name in ('calls.tsv', 'reconstructed_loci.fasta', 'molecule_candidate_evidence.tsv',
                 'molecule_repeat_likelihoods.tsv', 'short_read_qc_summary.tsv'):
        assert (tmp_path/'1'/name).read_bytes() == (tmp_path/'16'/name).read_bytes()
    metadata = json.loads(outputs[1]['run_metadata'].read_text())
    assert metadata['performance']['io']['input_backends'] == ['rapidgzip+htslib'] * 2
