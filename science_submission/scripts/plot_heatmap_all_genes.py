"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""

import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels, flip_ticks, cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_all_genes(axMain, axLabel, label_size=5, **kwargs):
    zscores = pd.read_parquet('../scrnaseq-wf/data/tpm_zscore.parquet')

    # calculate linkages
    link = linkage(zscores.values, 'average')
    tree = dendrogram(link, no_plot=True)
    leaves = tree['leaves']

    # plot
    defaults = dict(yticklabels=False, xticklabels=False, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)

    sns.heatmap(zscores.iloc[leaves], ax=axMain, **defaults)
    add_color_labels(axLabel, s=label_size)

    axMain.set_ylabel('')

    # Add call out box annotations
    box_defaults = dict(edgecolor='k', facecolor='none', linewidth=2, clip_on=False)
    txt_defaults = dict(xytext=(0, 8), textcoords='offset points', color='k', fontsize=8, fontweight='bold',
                        ha='center', bbox=dict(color='w', alpha=.6, pad=0))

    axMain.add_artist(plt.Rectangle((0, 9000), 1, 1000, **box_defaults))
    axMain.annotate('1', xy=(0.5, 9000), **txt_defaults)

    axMain.add_artist(plt.Rectangle((2, 2600), 1, 500, **box_defaults))
    axMain.annotate('2', xy=(2.5, 2600), **txt_defaults)

    axMain.add_artist(plt.Rectangle((3, 2100), 1, 500, **box_defaults))
    axMain.annotate('3', xy=(3.5, 2100), **txt_defaults)

    axMain.add_artist(plt.Rectangle((4, 2300), 1, 500, **box_defaults))
    axMain.annotate('4', xy=(4.5, 2300), **txt_defaults)

    axMain.add_artist(plt.Rectangle((10, 11300), 1, 1000, **box_defaults))
    axMain.annotate('5', xy=(10.5, 11300), **txt_defaults)

    axMain.add_artist(plt.Rectangle((11, 10800), 1, 1000, **box_defaults))
    axMain.annotate('6', xy=(11.5, 10800), **txt_defaults)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 2, height_ratios=[1, 0.06], width_ratios=[0.1, 1], hspace=0, wspace=0.2)
    axMain = fig.add_subplot(gs[0, 1])
    axLabel = fig.add_subplot(gs[1, 1])
    axCbar = fig.add_subplot(gs[0, 0])
    plot_heatmap_all_genes(axMain, axLabel, label_size=10, cbar_ax=axCbar)

    flip_ticks(axCbar)
    flip_ticks(axMain, pos='right')
    plt.setp(axMain.get_yticklabels(), rotation=0)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
