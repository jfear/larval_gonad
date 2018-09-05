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
from larval_gonad.plotting import dechr
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_boxplot_all_genes_by_chrom(gs, box_kws=None, line_kws=None, stripplot_kws=None):
    box_kws = box_kws or {}
    line_kws = line_kws or {}
    stripplot_kws = stripplot_kws or {}

    # Get list of germ cells
    germ_cells = config['sel_cluster_order'][:4]

    # Skip Y chrom
    chrom_order = config['chrom_order'][:-1]

    # munge data
    tpm = np.log10(pd.read_parquet('../scrnaseq-wf/data/tpm.parquet', columns=germ_cells) + 1)\
        .join(fbgn2chrom).query(f'chrom == {chrom_order}')
    median_by_chrom = tpm.groupby('chrom').median()
    autosome_median = tpm.query('chrom == ["chr2L", "chr2R", "chr3L", "chr3R"]').median()

    melted = tpm.reset_index().melt(id_vars=['FBgn', 'chrom'], var_name='Cluster', value_name='logTPM')
    melted.Cluster = pd.Categorical(melted.Cluster, categories=germ_cells, ordered=True)

    # Create axes
    fig = plt.gcf()
    gs0 = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs, wspace=0.08)
    ax1 = fig.add_subplot(gs0[0, 0])
    ax2 = fig.add_subplot(gs0[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs0[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs0[0, 3], sharex=ax1, sharey=ax1)
    axes = [ax1, ax2, ax3, ax4]

    # Set up defaults for plotting
    box_defaults = dict(notch=True, showfliers=False, order=chrom_order, palette=config['colors']['chrom_boxplot'])
    box_defaults.update(box_kws)

    stripplot_defaults = dict(color='w', marker='d', linewidth=1, edgecolor='k', jitter=False)
    stripplot_defaults.update(stripplot_kws)

    line_defaults = dict(color='r', ls='--')
    line_defaults.update(line_kws)

    # Plot each cluster on different axes
    for ax, (g, dd) in zip(axes, melted.groupby('Cluster')):

        meds = median_by_chrom.loc[chrom_order, g].reset_index()
        auto_med = autosome_median[g]

        sns.boxplot('chrom', 'logTPM', data=dd, ax=ax, **box_defaults)

        # Add the little diamonds on the median
        sns.stripplot('chrom', g, data=meds, ax=ax, **stripplot_defaults)

        ax.axhline(auto_med, **line_defaults)
        ax.set_title(g, fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('TPM (log)')
        sns.despine(ax=ax)

    # Tweak secondary axes removing unecessary stuff
    for ax in axes[1:]:
        ax.set_ylabel('')
        ax.yaxis.set_visible(False)
        sns.despine(ax=ax, left=True)
        dechr(ax)

    return axes


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(1, 1)
    plot_boxplot_all_genes_by_chrom(gs[0, 0])
    fig.savefig(snakemake.output[0], bbox_inches='tight')
