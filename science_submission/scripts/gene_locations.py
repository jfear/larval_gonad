"""Create a table of gene locations."""
from pathlib import Path

import numpy as np
import pandas as pd
import gffutils

from larval_gonad.config import config

db = gffutils.FeatureDB(snakemake.input[0])

# organize genes by chrom position
locs = []
for gene in db.features_of_type('gene'):
    fbgn = gene.id
    chrom = gene.chrom
    pos = np.min([gene.start, gene.end])
    locs.append((fbgn, chrom, pos))
locs = pd.DataFrame(locs, columns=['FBgn', 'chrom', 'pos']).set_index('FBgn')
locs.chrom = pd.Categorical(locs.chrom, categories=config['chrom_order'], ordered=True)

# sort and make an arbitrary location point
locs.sort_values(['chrom', 'pos'], inplace=True)
locs['location'] = range(1, locs['pos'].shape[0] + 1)

locs.to_parquet(snakemake.output[0])
