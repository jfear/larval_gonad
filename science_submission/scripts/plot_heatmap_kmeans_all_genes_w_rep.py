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
from larval_gonad.plotting import add_color_labels_w_rep, flip_ticks, cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_kmeans_all_genes(axMain, **kwargs):
    zscores = pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet', columns=config['sel_cluster_order_w_rep'])

    # calculate linkages
    X = zscores.values
    km = KMeans(n_clusters=10, random_state=42)
    km.fit(X)
    _dat = zscores.assign(km=km.labels_ + 1).sort_values('km')

    # plot
    defaults = dict(yticklabels=False, xticklabels=True, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)

    sns.heatmap(_dat.drop('km', axis=1), ax=axMain, **defaults)
    axMain.set_ylabel('')

    # Add annotations
    loc = 0
    for g, dd in _dat.groupby('km'):
        loc += dd.shape[0]
        axMain.axhline(loc, color='w', lw=1)
        axMain.text(-.3, loc - (dd.shape[0] / 2), g, va='center', ha='right', fontsize=10, color='k')


if __name__ == '__main__':
    fig = plt.figure(figsize=(3, 8))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01 )
    axMain = fig.add_subplot(gs[0, 0])
    axCbar = fig.add_subplot(gs[1, 0])
    plot_heatmap_kmeans_all_genes(axMain, cbar_ax=axCbar, cbar_kws={'orientation': 'horizontal'})

    axMain.xaxis.set_ticks_position('top')
    axMain.xaxis.set_label_position('top')
    plt.setp(axMain.get_xticklabels(), rotation=50, ha='left', fontsize=6)


    fig.savefig(snakemake.output[0], bbox_inches='tight')
