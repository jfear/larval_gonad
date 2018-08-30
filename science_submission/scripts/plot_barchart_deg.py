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
from larval_gonad.plotting import dechr
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')


def parse_genes(g, chroms):
    # import each upregulated dataset and make a list of genes to use
    upgenes = {
        'Spermatogonia': {
            'fname': '../scrnaseq-wf/data/gonia_vs_cytes.tsv',
            'query': 'p_val_adj <= 0.01 & avg_logFC > 0',
        },
        'Early 1º Spermatocytes': {
            'fname': '../scrnaseq-wf/data/gonia_vs_early.tsv',
            'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        },
        'Mid 1º Spermatocytes': {
            'fname': '../scrnaseq-wf/data/early_vs_mid.tsv',
            'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        },
        'Late 1º Spermatocytes': {
            'fname': '../scrnaseq-wf/data/mid_vs_late.tsv',
            'query': 'p_val_adj <= 0.01 & avg_logFC < 0',
        },
    }

    df = pd.read_csv(upgenes[g]['fname'], sep='\t')
    sig_genes = df.query(upgenes[g]['query']).set_index('primary_FBgn').index.tolist()

    # munge data for running fishers exact test
    sig_cnt = chroms.reindex(sig_genes).value_counts()
    sig_cnt.name = 'sig'

    # If there are no sig genes for a chrom make it 0
    sig_cnt = sig_cnt.reindex(chrom_order).fillna(0)

    not_sig_cnt = chroms.value_counts().reindex(chrom_order) - sig_cnt
    not_sig_cnt.name = 'ns'

    xsig = sig_cnt['chrX']
    xns = not_sig_cnt['chrX']

    foursig = sig_cnt['chr4']
    fourns = not_sig_cnt['chr4']

    asig = sig_cnt[['chr2L', 'chr2R', 'chr3L', 'chr3R']].sum()
    ans = not_sig_cnt[['chr2L', 'chr2R', 'chr3L', 'chr3R']].sum()

    # Run stats to check if depletion of DEG genes on chrX and chr4
    _, xpval = fisher_exact(np.array([[xsig, xns], [asig, ans]]), 'less')
    _, fourpval = fisher_exact(np.array([[foursig, fourns], [asig, ans]]), 'less')

    # Make proprotional for plotting
    prop = sig_cnt.div(not_sig_cnt).to_frame()
    prop.index.name = 'chrom'
    prop.columns = ['prop']

    return xpval, fourpval, prop


def plot_barchart_deg(gs):
    # Get list of germ cells
    germ_cells = config['sel_cluster_order'][:4]

    # Skip Y chrom
    global chrom_order
    chrom_order = config['chrom_order'][:-1]

    # munge data to get list of expressed genes
    tpm = np.log10(pd.read_parquet('../scrnaseq-wf/data/tpm.parquet', columns=germ_cells) + 1)\
        .join(fbgn2chrom).query(f'chrom == {chrom_order}')

    expressed = tpm.index[tpm.iloc[:, :4].sum(axis=1) > 0].tolist()
    chroms = fbgn2chrom.reindex(expressed).chrom

    # Create axes
    fig = plt.gcf()
    gs0 = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs, wspace=0.08)
    ax1 = fig.add_subplot(gs0[0, 0])
    ax2 = fig.add_subplot(gs0[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs0[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs0[0, 3], sharex=ax1, sharey=ax1)
    axes = [ax1, ax2, ax3, ax4]

    # Plot each cluster on different axes
    for ax, g, in zip(axes, germ_cells):
        xpval, fourpval, prop = parse_genes(g, chroms)

        sns.barplot('chrom', 'prop', data=prop.reset_index(), order=chrom_order,
                    palette=config['colors']['chrom_boxplot'], edgecolor='k', linewidth=1.2, ax=ax)

        ax.set_title(g, fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('Proportion Genes (DEG)')
        sns.despine(ax=ax)

        # Add astrics for chrX
        if xpval <= 0.001:
            ax.text(0, prop.prop['chrX'], '**', ha='center', va='bottom')
        elif xpval <= 0.05:
            ax.text(0, prop.prop['chrX'], '*', ha='center', va='bottom')

        # Add astrics for chr4
        if fourpval <= 0.001:
            ax.text(5, prop.prop['chr4'], '**', ha='center', va='bottom')
        elif fourpval <= 0.05:
            ax.text(5, prop.prop['chr4'], '*', ha='center', va='bottom')

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
    plot_barchart_deg(gs[0, 0])
    fig.savefig(snakemake.output[0], bbox_inches='tight')
