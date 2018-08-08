"""Plot heatmap of literature genes."""

import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels
from larval_gonad.normalization import zscore
import common

zscores = pd.read_parquet(snakemake.input.zscore)
ptrap_scores = pd.read_parquet(snakemake.input.ptrap_score)
ptrap_scores.index = ptrap_scores.index.droplevel(0)
ptrap_genes = ptrap_scores.index

oname = snakemake.output[0]

# Pull out lit genes
zscores.index = zscores.index.map(common.fbgn2symbol)
zscores = zscores.reindex(ptrap_genes)

# calculate linkages
link = linkage(zscores.values, 'average')
tree = dendrogram(link, no_plot=True)
leaves = tree['leaves']

# plot
fig = plt.figure(figsize=(2, 2))
gs = GridSpec(6, 2, width_ratios=[0.1, 1], wspace=0.04, hspace=0.2, left=.2, right=.6)
cax = fig.add_subplot(gs[:, 0])

genes = [
    'osa',
    'Mapmodulin',
    'SRPK',
    'bol',
    'Fas3',
    'cindr',
]
axes = []
for i, gene in zip(range(6), genes):
    gs00 = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[i, 1], hspace=0)
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[1, 0])

    sns.heatmap(zscores.loc[gene].to_frame().T,
                cmap=config['colors']['heatmap'],
                yticklabels=False,
                xticklabels=False,
                ax=ax1,
                cbar=False,
                vmin=-3,
                vmax=3,
                )

    sns.heatmap(ptrap_scores.loc[gene].to_frame().T,
                cmap='inferno',
                yticklabels=False,
                xticklabels=False,
                ax=ax2,
                cbar=True,
                cbar_ax=cax,
                vmin=0,
                vmax=5,
                )

    ax1.text(1.01, 0, gene, transform=ax1.transAxes, ha='left', va='center', fontsize=10)

cax.yaxis.set_ticks_position('left')
cax.yaxis.set_label_position('left')
cax.set_ylabel('Arbitrary Score')

fig.savefig(oname)
