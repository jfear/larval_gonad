#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io

#%%
plt.rcParams['figure.figsize'] = (8, 8)

#%%
SAMPLE = 'testis1'
MTX = 'output/cellranger-wf/testis1/outs/raw_gene_bc_matrices/dmelr6-26/matrix.mtx'
CELL_IDS = 'output/cellranger-wf/testis1/outs/raw_gene_bc_matrices/dmelr6-26/cell_ids.tsv'
CELL_CALLS = 'output/cellselection-wf/testis1_combined_cell_calls.feather'
DOUBLETS = 'output/cellselection-wf/testis1_scrublet_dublets.txt'

#%%
dat = np.asarray(scipy.io.mmread(MTX).sum(axis=0))[0]

#%%
with open(CELL_IDS) as fh:
    cell_ids = fh.read().strip().split('\n')

#%%
cell_calls = pd.read_feather(CELL_CALLS).set_index('cell_id').sum(axis=1) >= 3
cell_calls.name = 'is_cell'

#%%
with open(DOUBLETS) as fh:
    dubs = fh.read().strip().split('\n')

#%%
df = pd.Series(data=dat, index=cell_ids).sort_values(ascending=False).rename('Total UMI').to_frame().rename_axis('cell_id')
df['Cell Rank'] = range(1, df.shape[0] + 1)

#%%
data = df.join(cell_calls).fillna(False).assign(doublet = False)
data.loc[data.index.isin(dubs), 'doublet'] = True

data['color'] = 'lightgray'
data.loc[data.is_cell, 'color'] = 'k'
data.loc[data.doublet, 'color'] = 'r'

data['sizes'] = 1
data.loc[data.doublet, 'sizes'] = 20

#%%
ax = data.plot('Cell Rank', 'Total UMI', kind='scatter', color=data.color, s=data.sizes, rasterized=True)
ax.set(xscale='symlog', yscale='symlog', title=SAMPLE);

#%%
