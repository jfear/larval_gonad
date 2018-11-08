"""Plot barchart Spearman correlation with bulk RNA-Seq data."""

import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from larval_gonad.config import read_config, config
from larval_gonad.plotting import cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')

config.update(read_config('config.yaml'))


def plot_barchart_bulk_corr(gs, **kwargs):
    sc = pd.read_parquet('../output/scrnaseq-wf/tpm.parquet', columns=config['sel_cluster_order'])

    bulk = pd.read_parquet('../output/bulk-rnaseq-wf/aggregation/tpm_gene_level_counts.parquet',
                           columns=['C1_TDT', 'C2_TDT', 'C3_TDT', 'C4_TDT'])
    bulk = bulk.mean(axis=1)
    bulk.name = 'bulk'

    # I want a vector of correlation coeff with bulk.
    corr = sc.join(bulk).corr(method='spearman')
    dat = corr.bulk.drop('bulk')

    # Make colors for plotting
    colors = [cluster_cmap[x] for x in dat.index]

    # Setup borken axis
    fig = plt.gcf()
    broken = mpl.gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs, hspace=0.1)
    ax1 = fig.add_subplot(broken[0, 0])
    ax2 = fig.add_subplot(broken[1, 0], sharex=ax1)

    # plot
    dat.plot.bar(color=colors, ax=ax1)
    ax1.set_axisbelow(True)
    ax1.grid(linestyle=':', axis='y')

    dat.plot.bar(color=colors, ax=ax2)
    ax2.set_axisbelow(True)
    ax2.grid(linestyle=':', axis='y')

    # zoom axis in
    ax1.set_ylim(.7, 1)
    ax2.set_ylim(0, .3)

    # hide the spines between ax and ax2
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.set_visible(False)

    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    return ax1, ax2


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = mpl.gridspec.GridSpec(1, 1)
    ax1, ax2 = plot_barchart_bulk_corr(gs[0, 0])
    ax2.set_xlabel('Cell Cluster')
    fig.text(0.03, 0.5, 'Correlation with Bulk RNA-Seq (Spearman)', ha='center', va='center', rotation=90)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
