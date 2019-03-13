"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use('Agg')

import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

tpm = snakemake.input[0]
annotation = snakemake.params.annotation
cluster_colors = snakemake.params.cluster_colors
oname = snakemake.output[0]


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig = plt.figure(figsize=(1.8, 1.8))
    gs = GridSpec(2, 1, height_ratios=[1, .05])
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        df,
        xticklabels=True,
        yticklabels=True,
        vmin=0.75,
        vmax=1,
        square=True,
        rasterized=True,
        ax=ax,
        annot=True,
        annot_kws=dict(fontsize=4.5, fontweight='bold'),
        cbar_ax=cax,
        cbar_kws=dict(label='Spearman Correlation', orientation='horizontal', ticks=[.75, 1])
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')

    # Clean up Y axis
    ax.set_ylabel('')
    ax.yaxis.set_ticks_position('right')
    plt.setp(ax.get_yticklabels(), rotation=0, ha='left', va='center')

    # Add separation lines
    ax.axhline(4, color='w', ls='--', lw=.5)
    ax.axvline(4, color='w', ls='--', lw=.5)

    # Add Annotation
    ax.text(1.2, 2.8, 'Germline', rotation=-45, ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax.text(5.7, 7.3, 'Soma', rotation=-45, ha='center', va='center', fontsize=6,
            color=cluster_colors[4])

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    # scRNA-seq data
    df = (
        pd.read_parquet(tpm)
        .pivot_table(index='FBgn', columns='cluster', values='TPM')
        .rename(columns=annotation)
        .drop(columns='UNK')
    )

    _corr = df.corr(method='spearman')

    # calculate linkages
    link = linkage(_corr.values, 'average')
    tree = dendrogram(link, no_plot=True)
    leaves = tree['leaves']

    _corr = _corr.iloc[leaves, leaves]
    for i, j in zip(*np.tril_indices(_corr.shape[0], k=-1)):
        _corr.iloc[i, j] = np.nan

    return _corr


if __name__ == '__main__':
    main()
