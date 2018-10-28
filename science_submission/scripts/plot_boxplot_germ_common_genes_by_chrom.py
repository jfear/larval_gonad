"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import config
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_boxplot_germ_common_genes_by_chrom(gs):
    # Get list of germ cells
    germ_cells = config['cluster_order'][:4]

    # Combine autosomes into 'A'
    chroms = fbgn2chrom.copy()
    chroms.loc[chroms['chrom'] == 'chrX', 'XA'] = 'X'
    chroms.loc[chroms['chrom'] == 'chr4', 'XA'] = '4'
    chroms.loc[chroms['chrom'] == 'chrY', 'XA'] = 'Y'
    chroms.loc[chroms.chrom.isin(['chr2L', 'chr2R', 'chr3L', 'chr3R']), 'XA'] = 'A'

    # munge data
    tpm = np.log10(pd.read_parquet('../output/scrnaseq-wf/tpm.parquet', columns=germ_cells) + 1).join(chroms)
    melted = tpm.reset_index().melt(id_vars=['FBgn', 'chrom', 'XA'], var_name='Cluster', value_name='logTPM')
    melted.Cluster = pd.Categorical(melted.Cluster, categories=germ_cells, ordered=True)

    # For each gene is it expressed (>0) in more than 1/3 of germline clusters
    common = tpm.index[(tpm.iloc[:, :-2] > 0).sum(axis=1) > ((tpm.shape[1] - 2) / 3)].tolist()
    melted = melted.query(f'FBgn == {common}')

    # Create axes
    fig = plt.gcf()
    gs0 = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs, wspace=0.08)
    ax1 = fig.add_subplot(gs0[0, 0])
    ax2 = fig.add_subplot(gs0[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs0[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs0[0, 3], sharex=ax1, sharey=ax1)
    axes = [ax1, ax2, ax3, ax4]

    # Plot each cluster on different axes
    for ax, (g, dd) in zip(axes, melted.groupby('Cluster')):
        sns.boxplot('XA', 'logTPM', data=dd, notch=True, showfliers=False, order=['X', 'A', '4'], ax=ax)
        ax.set_title(g, fontsize=10)
        sns.despine(ax=ax)

    # Tweak secondary axes removing unecessary stuff
    for ax in axes[1:]:
        ax.set_ylabel('')
        ax.yaxis.set_visible(False)
        sns.despine(ax=ax, left=True)

    return axes


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(1, 1)
    plot_boxplot_germ_common_genes_by_chrom(gs[0, 0])
    fig.savefig(snakemake.output[0], bbox_inches='tight')
