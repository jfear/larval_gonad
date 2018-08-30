"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import config
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')


def run_test(up, down, dat):
    a = dat.loc[["chr2L", "chr2R", "chr3L", "chr3R"]].sum()
    ab = a[f'{down} bias']
    an = a.sum() - ab

    x = dat.loc["chrX"]
    xb = x[f'{down} bias']
    xn = x.sum() - xb

    t, p = fisher_exact(np.array([[xb, xn], [ab, an]]), 'less')
    return p


def _plot(up, down, ax, chroms):
    df = pd.read_csv(f'../scrnaseq-wf/data/{up}_vs_{down}.tsv', sep='\t', index_col=0).query('p_val_adj <= 0.01')
    df.index.name = 'FBgn'

    genes_up = df.query('avg_logFC > 0').index.tolist()
    genes_down = df.query('avg_logFC < 0').index.tolist()

    _chroms = chroms.copy()
    _chroms['sig'] = 'non-bias'
    _chroms.loc[chroms.index.isin(genes_up), 'sig'] = f'{up} bias'
    _chroms.loc[chroms.index.isin(genes_down), 'sig'] = f'{down} bias'

    dat = _chroms.groupby('chrom').sig.value_counts()
    dat.columns = ['count']
    dat = dat.unstack().fillna(0)

    _dat = dat.div(dat.sum(axis=1), axis=0).loc[config['chrom_order'][:5], [f'{down} bias', 'non-bias', f'{up} bias']]
    _dat.plot.bar(stacked=True, ax=ax, color=['b', 'grey', 'r'], width=.9)
    ax.set_title(f'{up} vs {down}')
    ax.legend_.set_visible(False)
    plt.setp(ax.get_xticklabels(), rotation=0)

    p = run_test(up, down, dat)
    if p <= 0.01:
        y = _dat.loc['chrX', f'{down} bias']
        ax.text(0, y, '*', color='white', va='bottom', ha='center')


def plot_stacked_bar_prop_by_chrom(gs):
    # Get list of germ cells
    germ_cells = config['cluster_order'][:4]

    # munge data
    germ_cells = config['cluster_order'][:4]
    tpm = pd.read_parquet('../scrnaseq-wf/data/tpm.parquet', columns=germ_cells)
    in_data = tpm.index[(tpm.loc[:, germ_cells].sum(axis=1) > 0)].tolist()
    chroms = fbgn2chrom.copy().query('chrom == ["chrX", "chr2L", "chr2R", "chr3L", "chr3R"]').reindex(in_data).dropna()

    # Create axes
    fig = plt.gcf()
    gs0 = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs, wspace=0.08)
    ax1 = fig.add_subplot(gs0[0, 0])
    ax2 = fig.add_subplot(gs0[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs0[0, 2], sharex=ax1, sharey=ax1)
    axes = [ax1, ax2, ax3]

    # Plot each cluster on different axes
    for ax, (up, down) in zip(axes, [('gonia', 'early'), ('early', 'mid'), ('mid', 'late')]):
        _plot(up, down, ax, chroms)

    # Tweak secondary axes
    ax1.set_ylabel('Proportion of Genes Differentially Expressed')

    sns.despine(ax=ax1)
    sns.despine(ax=ax2, left=True)
    sns.despine(ax=ax3, left=True)

    ax2.yaxis.set_visible(False)
    ax3.yaxis.set_visible(False)

    return axes


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(1, 1)
    plot_stacked_bar_prop_by_chrom(gs[0, 0])
    fig.savefig(snakemake.output[0], bbox_inches='tight')
