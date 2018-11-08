"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns


mpl.style.use('scripts/paper_1c.mplstyle')
fmt = mpl.ticker.FuncFormatter(lambda x, pos: int(np.abs(x)))


def tweak_plots(ax, x=False, title=False, legend=False):
    if x:
        ax.set_xticklabels([])
        ax.set_xlabel('')

    if title:
        ax.set_title('')

    if legend:
        ax.legend_.set_visible(False)


def _plot(up, down, locs, ax):
    dat = pd.read_csv(f'../output/scrnaseq-wf/{up}_vs_{down}.tsv', sep='\t', index_col=0)
    dat['logpval'] = -np.log10(dat.p_val_adj.replace({0: 1e-300}))
    dat.loc[dat.avg_logFC < 0, 'logpval'] *= -1
    dat = dat.join(locs, how='left')[['chrom', 'location', 'logpval']]
    dat.logpval = dat.logpval.fillna(0)

    ax.yaxis.set_major_formatter(fmt)
    for g, dd in dat.groupby('chrom'):
        ax.scatter(dd.location, dd.logpval, label=g, s=5, rasterized=True)

    ax.legend(loc='center right', bbox_to_anchor=[-.15, 0.5])
    ax.text(-0.1, 0.75, f'{up}', rotation=90, transform=ax.transAxes, ha='right', va='center')
    ax.text(-0.1, 0.25, f'{down}', rotation=90, transform=ax.transAxes, ha='right', va='center')

    ax.set_axisbelow(True)
    ax.axhline(0, color='k', alpha=.5, ls=':')
    ax.grid(axis='y', alpha=0.5, ls=':')
    sns.despine(ax=ax, left=True)


def plot_manhatten(gs):
    # Get list of germ cells
    locs = pd.read_parquet('../output/science_submission/gene_locations.parquet')
    tpm = pd.read_parquet('../output/scrnaseq-wf/tpm.parquet')
    locs = locs.reindex(tpm.index).sort_values(['chrom', 'pos'])

    # munge data
    comparisons = [
        ('gonia', 'early'),
        ('early', 'mid'),
        ('mid', 'late')
    ]

    # Create axes
    fig = plt.gcf()
    gs0 = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs, wspace=0.08)
    ax1 = fig.add_subplot(gs0[0, 0])
    ax2 = fig.add_subplot(gs0[1, 0], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs0[2, 0], sharex=ax1, sharey=ax1)
    axes = [ax1, ax2, ax3]

    for (up, down), ax in zip(comparisons, axes):
        _plot(up, down, locs, ax)

    tweak_plots(ax1, x=True, legend=True)
    tweak_plots(ax2, x=True)
    tweak_plots(ax3, legend=True)

    return axes


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(1, 1)
    plot_manhatten(gs[0, 0])
    fig.savefig(snakemake.output[0], bbox_inches='tight')
