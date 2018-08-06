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
gs = GridSpec(6, 1, hspace=0.2, right=.62)

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
    gs00 = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[i, 0], hspace=0)
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[1, 0])

    sns.heatmap(zscores.loc[gene].to_frame().T, cmap=config['colors']['heatmap'], yticklabels=False, xticklabels=False,
                ax=ax1, cbar=False, vmin=-3, vmax=3)

    sns.heatmap(ptrap_scores.loc[gene].to_frame().T, cmap='Greys', yticklabels=False, xticklabels=False,
                ax=ax2, cbar=False, vmin=0)
    ax1.text(1.01, 0, gene, transform=ax1.transAxes, ha='left', va='center', fontsize=10)
#
# #add_color_labels(axLabel)
#
# axZ.yaxis.set_ticks_position('right')
# axZ.yaxis.set_label_position('right')
# axZ.set_ylabel('')
# plt.setp(axZ.get_yticklabels(), rotation=0, fontsize=10)
#
# axP.yaxis.set_ticks_position('right')
# axP.yaxis.set_label_position('right')
# axP.set_ylabel('')
# plt.setp(axP.get_yticklabels(), rotation=0, fontsize=10)

fig.savefig(oname)
# fig = plt.figure(figsize=(8, 8))
# gs = GridSpec(2, 2, height_ratios=[1, .08], hspace=0)
# axMain = plt.subplot(gs[0, 0])
# axLabel = plt.subplot(gs[1, 0])
#
# axScore = plt.subplot(gs[0, 1])
# axScoreLabel = plt.subplot(gs[1, 1])
#
# sns.heatmap(zscores.iloc[leaves], cmap=config['colors']['heatmap'], yticklabels=True, xticklabels=False,
#             vmin=-3, vmax=3, ax=axMain, cbar=False)
#
# sns.heatmap(ptrap_scores.iloc[leaves], cmap='Greys', yticklabels=True, xticklabels=False,
#             ax=axScore, cbar=False, vmin=0)
#
# add_color_labels(axLabel)
# add_color_labels(axScoreLabel)
#
# axMain.yaxis.set_ticks_position('right')
# axMain.yaxis.set_label_position('right')
# axMain.set_ylabel('')
# plt.setp(axMain.get_yticklabels(), rotation=0, fontsize=10)
#
# axScore.yaxis.set_ticks_position('right')
# axScore.yaxis.set_label_position('right')
# axScore.set_ylabel('')
# plt.setp(axScore.get_yticklabels(), rotation=0, fontsize=10)
#
# fig.savefig(oname)
