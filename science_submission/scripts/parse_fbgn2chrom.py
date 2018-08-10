#!/usr/bin/env python3
"""Quick script to make a table with FBgn to chromosome."""

import pandas as pd


gtf = snakemake.input[0]
oname = snakemake.output[0]

genes = []
with open(gtf) as fh:
    for row in fh:
        rows = row.strip().split()

        if len(rows) == 0:
            continue

        if rows[2] == 'gene':
            genes.append((rows[0], rows[9].replace('"', '').replace(';', '')))

fbgn2chrom = pd.DataFrame(genes, columns=['chrom', 'FBgn'])
fbgn2chrom.set_index('FBgn', inplace=True)

CHROMS = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY', 'chrM']
fbgn2chrom = fbgn2chrom[fbgn2chrom['chrom'].isin(CHROMS)]

fbgn2chrom.to_parquet(oname)
