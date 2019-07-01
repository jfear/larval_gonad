"""Create a table of gene locations."""
import numpy as np
import pandas as pd
import gffutils

from larval_gonad.config import config

fname = snakemake.input[0]
oname = snakemake.output[0]
chrom_order = snakemake.params.chrom_order


# Parse GTF to get gene and chrom position
db = gffutils.FeatureDB(fname)
locs = []
for gene in db.features_of_type('gene'):
    fbgn = gene.id
    chrom = gene.chrom
    pos = np.min([gene.start, gene.end])
    locs.append((fbgn, chrom, pos))

# Order by chrom and position
(
    pd.DataFrame(locs, columns=['FBgn', 'chrom', 'pos'])
    .set_index('FBgn')
    .query(f'chrom == {chrom_order}')
    .assign(chrom=lambda df: df.chrom.astype('category').cat.as_ordered().cat.reorder_categories(chrom_order))
    .sort_values(['chrom', 'pos'])
    .assign(location=lambda df: range(1, df.shape[0] + 1))
    .reset_index()
    .to_feather(oname)
)
