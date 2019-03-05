"""Create a dictionary mapping fbgn to gene symbol."""

from pathlib import Path
from pickle import dump

import pandas as pd

oname = snakemake.output[0]
assembly = snakemake.params.assembly
tag = snakemake.params.tag
references_dir = snakemake.params.references_dir
annot_fn = Path(references_dir, assembly, tag, 'fb_annotation', f'{assembly}_{tag}.fb_annotation')

df = pd.read_csv(annot_fn, sep='\t', index_col=1).fillna('nan')
fbgn2symbol = df['gene_symbol'].to_dict()

with open(oname, 'wb') as fo:
    dump(fbgn2symbol, fo)
