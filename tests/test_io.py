import sys
import os

TESTDIR = os.path.dirname(__file__)
DATADIR = os.path.join(TESTDIR, 'data')
LIBDIR = os.path.join(TESTDIR, '../lib')

sys.path.insert(0, LIBDIR)
from larval_gonad.io import FastqSummary


def test_FastqSummary_I1_gz():
    fname = os.path.join(DATADIR, 'I1.fastq.gz')
    fq = FastqSummary(fname)
    assert fq.num_reads == 20
    assert fq.avg_read_length == 8.0


def test_FastqSummary_I1():
    fname = os.path.join(DATADIR, 'I1.fastq')
    fq = FastqSummary(fname)
    assert fq.num_reads == 20
    assert fq.avg_read_length == 8.0


def test_FastqSummary_R1_gz():
    fname = os.path.join(DATADIR, 'R1.fastq.gz')
    fq = FastqSummary(fname, split=True)
    assert fq.num_reads == 20
    assert fq.avg_read_length == 26.0
    assert len(fq.cells) == 20


def test_FastqSummary_R1():
    fname = os.path.join(DATADIR, 'R1.fastq')
    fq = FastqSummary(fname, split=True)
    assert fq.num_reads == 20
    assert fq.avg_read_length == 26.0
    assert len(fq.cells) == 20


def test_FastqSummary_R2_gz():
    fname = os.path.join(DATADIR, 'R2.fastq.gz')
    fq = FastqSummary(fname)
    assert fq.num_reads == 20
    assert fq.avg_read_length == 98.0


def test_FastqSummary_R2():
    fname = os.path.join(DATADIR, 'R2.fastq')
    fq = FastqSummary(fname)
    assert fq.num_reads == 20
    assert fq.avg_read_length == 98.0

