"""Plot tSNE"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.plotting import cluster_cmap
import common

tsne = pd.read_csv(snakemake.input.tsne, sep='\t', index_col=0)
clusters = pd.read_parquet(snakemake.input.clusters)
oname = snakemake.output[0]

colors = clusters['cluster_name'].map(cluster_cmap)
colors.name = 'colors'

dat = tsne.join(colors)
fig, ax = plt.subplots(1, 1, figsize=(2, 2))
dat.plot('tSNE_1', 'tSNE_2', c=dat['colors'], s=3, linewidth=.1, edgecolor='k', kind='scatter', ax=ax)
sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal')

# Move X and Y labels to decrease clutter
ax.set_xlabel('tSNE_1', x=.15)
ax.set_ylabel('tSNE_2', y=.15)

# Add basic labels to distinguish cell types
_color = cluster_cmap['Early Cyst Cells (5)']
ax.text(0.1, .99, 'Cyst Lineage', transform=ax.transAxes, fontdict={'weight': 'bold', 'color': _color})

_color = cluster_cmap['Spermatogonia (6)']
ax.text(1, 0.03, 'Germline Lineage', va='top', ha='right', transform=ax.transAxes,
        fontdict={'weight': 'bold', 'color': _color})

fig.subplots_adjust(left=0.08, top=.95, bottom=.1, right=1)
fig.savefig(oname)
