#!/usr/bin/env python
"""Useful I/O related tools.

Tools for parsing input and writing output.
"""
import gzip
import sqlite3

class SingleCellFastqSummary(object):
    def __init__(self, db, I1=None, R1=None, R2=None):
        """Class to summarize information about a FASTQ formated file.

        Parameters:
        -----------

        Attributes:
        -----------
        """
        # connect to database
        self.conn = sqlite3.connect(db)
        self._create_tables()

        # Parse the files and add to database
        self.parse_fastq(I1, 'i1')
        self.parse_fastq(R1, 'r1')
        self.parse_fastq(R2, 'r2')


    def parse_fastq(self, fq, ftype='i1'):
        if fq is None:
            return

        c = self.conn.cursor()
        opener, read_type = self._opener(fq)
        with opener(fq, read_type) as fh:
            accumulator = []
            for record in self._fq_gen(fh):
                curr = []
                readID, seq, score = record
                readLen = len(seq)

                curr = [readID, seq, score, readLen]

                if ftype == 'r1':
                    bc, umi = record[1][1:14], record[1][14:]
                    curr.extend([bc, umi])

                accumulator.append(curr)

                if len(accumulator) >= 1e5:
                    if ftype == 'i1':
                        c.executemany('INSERT INTO i1 VALUES (?, ?, ?, ?)', accumulator)
                    elif ftype == 'r1':
                        c.executemany('INSERT INTO r1 VALUES (?, ?, ?, ?, ?, ?)', accumulator)
                    elif ftype == 'r2':
                        c.executemany('INSERT INTO r2 VALUES (?, ?, ?, ?)', accumulator)
                    self.conn.commit()
                    accumulator = []


    def _create_tables(self):
        c = self.conn.cursor()

        # Create tables if they don't exists
        c.execute("""CREATE TABLE IF NOT EXISTS i1
                     (readID text PRIMARY KEY, seq text, score text, readLen int)""")

        c.execute("""CREATE TABLE IF NOT EXISTS r1
                     (readID text PRIMARY KEY, seq text, score text, readLen int, barcode text, umi text)""")

        c.execute("""CREATE TABLE IF NOT EXISTS r2
                     (readID text PRIMARY KEY, seq text, score text, readLen int)""")

        # Commit Changes
        self.conn.commit()

    def _opener(self, file):
        if file.endswith('.gz'):
            fq_open = gzip.open
            read_type = 'rt'
        else:
            fq_open = open
            read_type = 'r'

        return fq_open, read_type

    def _fq_gen(self, fh):
        accumulator = []
        for i, row in enumerate(fh):
            row = row.strip()
            if i % 4 == 2:
                continue
            elif i % 4 == 3:
                accumulator.append(row)
                yield accumulator
                accumulator = []
            else:
                accumulator.append(row)

    def __str__(self):
        _str = ''
        cnt = 0
        for k, v in self.cells.items():
            _str += '{key}:\t{value}\n'.format(key=k, value=', '.join(v))
            cnt += 1
            if cnt > 20:
                _str += '...'
                break
        return _str
