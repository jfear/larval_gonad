#!/usr/bin/env python
"""Useful I/O related tools.

Tools for parsing input and writing output.
"""
import gzip
from collections import defaultdict

class FastqSummary(object):
    def __init__(self, file, split=False):
        """Class to summarize information about a FASTQ formated file.

        Parameters:
        -----------

        Attributes:
        -----------
        """
        self.cells = defaultdict(list)
        self.num_reads = 0
        self._sum_len = 0
        self._fh = self._open_file(file)
        for record in self._fq_gen():
            self.num_reads += 1
            self._sum_len += len(record[1])
            if split:
                bc, umi = record[1][1:14], record[1][14:]
                self.cells[bc].append(record[0])

        self.avg_read_length = self._sum_len / self.num_reads
        self._fh.close()
        del self._fh

    def _open_file(self, file):
        if file.endswith('.gz'):
            fq_open = gzip.open
            read_type = 'rt'
        else:
            fq_open = open
            read_type = 'r'

        return fq_open(file, read_type)

    def _fq_gen(self):
        accumulator = []
        for i, row in enumerate(self._fh):
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
