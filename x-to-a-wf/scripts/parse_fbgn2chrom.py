#!/usr/bin/env python3
"""Quick script to make a table with FBgn to chromosome."""
import pandas as pd

gtf = snakemake.input[0]
oname = snakemake.output[0]
CHROMS = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrY', 'chrM']

genes = []
with open(gtf) as fh:
    for row in fh:
        rows = row.strip().split()

        if len(rows) == 0:
            continue

        if rows[2] == 'gene':
            genes.append((rows[0], rows[9].replace('"', '').replace(';', '')))


(
    pd.DataFrame(genes, columns=['chrom', 'FBgn'])
    .set_index('FBgn')
    .query(f'chrom == {CHROMS}')
    .to_parquet(oname)
)
