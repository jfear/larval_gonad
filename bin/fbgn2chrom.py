#!/usr/bin/env python3
"""Quick script to make a table with FBgn to chromosome."""
import os
import pandas as pd
from yaml import load

CHROMS = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chrY', 'chrM']
OUT = '../output/fbgn2chrom.tsv'

with open('../config/common.yml') as fh:
    CONFIG = load(fh)


def main():
    assembly = CONFIG['assembly']
    tag = CONFIG['tag']
    gtf = os.path.join(os.environ['REFERENCES_DIR'], assembly, tag, 'gtf',
                       f'{assembly}_{tag}.gtf')
    genes = []
    with open(gtf) as fh:
        for row in fh:
            rows = row.strip().split()

            if len(rows) == 0:
                continue

            if rows[2] == 'gene':
                genes.append((rows[0],
                              rows[9].replace('"', '').replace(';', '')))

    fbgn2chrom = pd.DataFrame(genes, columns=['chrom', 'FBgn'])
    fbgn2chrom.set_index('FBgn', inplace=True)
    fbgn2chrom = fbgn2chrom[fbgn2chrom['chrom'].isin(CHROMS)]

    fbgn2chrom.to_csv(OUT, sep='\t')


if __name__ == '__main__':
    main()
