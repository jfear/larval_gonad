"""Plot heatmap of protein trap genes.

Plots a heatmap of protein trap scores. Protein traps were selected based on
their presence in the biomarkers identified by Seurat. Protein traps were
imaged and scored by multiple people and were summarized by Miriam into a score
ranging from 0-10. There were no 10s so I scaled down to 0-8.
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels, flip_ticks
from common import fbgn2symbol

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_ptrap_scores(axMain, axLabel=None, label_size=5, **kwargs):
    ptrap_scores = pd.read_parquet('data/ptrap_scores.parquet')
    ptrap_scores.index = ptrap_scores.index.droplevel(0)

    # calculate linkages
    link = linkage(ptrap_scores.values, 'average')
    tree = dendrogram(link, no_plot=True)
    leaves = tree['leaves']

    # plot
    defaults = dict(yticklabels=True, xticklabels=False, vmin=0, vmax=8, rasterized=True, cmap='inferno',
                    cbar_kws=dict(label='Arbitrary Score'))

    defaults.update(kwargs)

    sns.heatmap(ptrap_scores.iloc[leaves], ax=axMain, **defaults)

    if axLabel is not None:
        add_color_labels(axLabel, s=label_size)

    axMain.set_ylabel('')


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(1, 2, width_ratios=[0.1, 1], hspace=0, wspace=0.2)
    axMain = fig.add_subplot(gs[0, 1])
    axCbar = fig.add_subplot(gs[0, 0])
    plot_heatmap_ptrap_scores(axMain, label_size=10, cbar_ax=axCbar, xticklabels=True)

    flip_ticks(axCbar)
    flip_ticks(axMain, pos='right')
    plt.setp(axMain.get_yticklabels(), rotation=0)

    fig.suptitle('Protein Trap Scores')

    fig.savefig(snakemake.output[0], bbox_inches='tight')
