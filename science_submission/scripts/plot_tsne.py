"""Plot tSNE"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.plotting import cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_tsne(ax, **kwargs):

    # Get data
    tsne = pd.read_parquet('data/tsne.parquet')
    clusters = pd.read_parquet('data/clusters.parquet')

    colors = clusters['cluster_name'].map(cluster_cmap)
    colors.name = 'colors'

    dat = tsne.join(colors)

    # Plot
    defaults = dict(s=3, linewidth=.1, edgecolor='k')
    defaults.update(kwargs)

    dat.plot('tSNE_1', 'tSNE_2', c=dat['colors'], kind='scatter', ax=ax, **defaults)
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
    ax.text(0.4, 0.03, 'Germline Lineage', va='top', ha='left', transform=ax.transAxes,
            fontdict={'weight': 'bold', 'color': _color})


if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_tsne(ax, s=20)
    fig.savefig(snakemake.output[0], bbox_inches='tight')
