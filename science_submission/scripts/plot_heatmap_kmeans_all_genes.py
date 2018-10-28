"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""

import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels, flip_ticks, cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_kmeans_all_genes(axMain, axLabel, label_size=5, **kwargs):
    zscores = pd.read_parquet('../output/scrnaseq-wf/tpm_zscore.parquet', columns=config['sel_cluster_order'])

    # calculate linkages
    X = zscores.values
    km = KMeans(n_clusters=10, random_state=42)
    km.fit(X)
    _dat = zscores.assign(km=km.labels_ + 1).sort_values('km')

    # plot
    defaults = dict(yticklabels=False, xticklabels=False, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)

    sns.heatmap(_dat.drop('km', axis=1), ax=axMain, **defaults)
    add_color_labels(axLabel, s=label_size)

    axMain.set_ylabel('')

    # Add annotations
    loc = 0
    for g, dd in _dat.groupby('km'):
        loc += dd.shape[0]
        axMain.axhline(loc, color='w', lw=1)
        axMain.text(0.3, loc - (dd.shape[0] / 2), g, va='center', ha='left',
                    fontsize=10, color='k', bbox=dict(color='w', alpha=.6, ec='none', boxstyle='square,pad=.3'))


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 2, height_ratios=[1, 0.06], width_ratios=[0.1, 1], hspace=0, wspace=0.2)
    axMain = fig.add_subplot(gs[0, 1])
    axLabel = fig.add_subplot(gs[1, 1])
    axCbar = fig.add_subplot(gs[0, 0])
    plot_heatmap_kmeans_all_genes(axMain, axLabel, label_size=10, cbar_ax=axCbar)

    flip_ticks(axCbar)
    flip_ticks(axMain, pos='right')
    plt.setp(axMain.get_yticklabels(), rotation=0)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
