"""Parse a bam file to get cell level counts by chromosomes"""
from typing import Dict, List
from collections import defaultdict, namedtuple

import numpy as np
import pandas as pd
import pysam

CellCounts = namedtuple('CellCounts', 'cell_id chromosome number_reads')
HEADER = 'cell_id\tchromosome\tnumber_reads\n'
BAM_FILE = snakemake.input[0]
COUNT_FILE = snakemake.output[0]


def main():
    sam_file_handler = pysam.AlignmentFile(BAM_FILE, 'rb')

    chrom_counts_by_cell = []
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, '2L'))
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, '2R'))
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, '3L'))
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, '3R'))
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, '4'))
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, 'X'))
    chrom_counts_by_cell.extend(fetch_reads_from_chromosome(sam_file_handler, 'Y'))

    all_counts_sorted = sorted(chrom_counts_by_cell, key=lambda x: (x[0], x[1]))
    del chrom_counts_by_cell

    write_output_to_file(all_counts_sorted, COUNT_FILE)


def fetch_reads_from_chromosome(fh, chromosome: str) -> List[CellCounts]:
    data = defaultdict(set)
    for row in fh.fetch(chromosome):
        try:
            cell_id = row.get_tag('CB').rstrip('-1')
            umi_id = row.get_tag('UB')
            data[cell_id].add(umi_id)
        except KeyError:
            # Skip row is 'CB' or 'UB' are not there.
            pass

    return count_reads_on_chromosome(data, chromosome)


def count_reads_on_chromosome(data: dict, chromosome: str) -> List[CellCounts]:
    return [
        CellCounts(cell_id, chromosome, str(len(umi_ids)))
        for cell_id, umi_ids in data.items()
    ]


def write_output_to_file(cell_counts: List[CellCounts], file_name: str):
    with open(file_name, 'w') as fh:
        fh.write(HEADER)
        for row in cell_counts:
            fh.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    main()
