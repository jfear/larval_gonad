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

bulk = snakemake.input.bulk
scrnaseq = snakemake.input.scrnaseq
oname = snakemake.output[0]


def main():
    df = get_data()

    plt.style.use('scripts/paper_1c.mplstyle')
    fig = plt.figure(figsize=(1.8, 1.8))
    gs = GridSpec(2, 1, height_ratios=[1, .05], hspace=0.1)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        df,
        xticklabels=True,
        yticklabels=True,
        vmin=0.75,
        vmax=1,
        square=True,
        rasterized=False,
        ax=ax,
        annot=True,
        annot_kws=dict(fontsize=5, fontweight='bold'),
        cbar_ax=cax,
        cbar_kws=dict(label='Spearman Correlation', orientation='horizontal', ticks=[.75, 1])
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_tick_params(pad=0, length=2)
    plt.setp(ax.get_xticklabels(), fontsize=7, rotation=90)

    # Clean up Y axis
    ax.set_ylabel('')
    ax.yaxis.set_ticks_position('right')
    ax.yaxis.set_tick_params(pad=0, length=2)
    plt.setp(ax.get_yticklabels(), fontsize=7, rotation=0, ha='left', va='center')

    # Add separation lines
    ax.axhline(3, color='w', ls='--', lw=.5)
    ax.axvline(3, color='w', ls='--', lw=.5)

    # Add Annotation
    ax.text(.5, 2.5, 'scRNA-Seq\nTriplicates', rotation=-45, ha='center', va='center', fontsize=6)
    ax.text(4, 6, 'Bulk RNA-Seq\nQuadruplicates', rotation=-45, ha='center', va='center', fontsize=6)

    # Clean up color bar
    plt.setp(cax.xaxis.label, fontsize=6)
    plt.setp(cax.get_xticklabels(), fontsize=5)
    cax.xaxis.set_tick_params(pad=0, length=2)

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    # Bulk data
    tdt = ['C1_TDT', 'C3_TDT', 'C2_TDT', 'C4_TDT']
    df_bulk = (
        pd.read_parquet(bulk, columns=tdt)
        .rename_axis('FBgn')
        .rename(columns={
            'C1_TDT': 'b1',
            'C2_TDT': 'b2',
            'C3_TDT': 'b3',
            'C4_TDT': 'b4',
        })
    )

    # scRNA-seq data
    df_sc = (
        pd.read_parquet(scrnaseq)
        .pivot_table(index='FBgn', columns='rep', values='TPM')
        .rename(columns={
            'rep1': 'sc1',
            'rep2': 'sc2',
            'rep3': 'sc3',
        })
    )

    # Merge
    df = df_sc.join(df_bulk, how='inner')
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
