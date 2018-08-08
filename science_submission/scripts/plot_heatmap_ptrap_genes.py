"""Plot heatmap of literature genes."""

import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels, flip_ticks
from common import fbgn2symbol

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_ptrap_genes(gsMain, axLabel, label_size=5, expression_heatmap_kws=None, ptrap_heatmap_kws=None):
    zscores = pd.read_parquet('../scrnaseq-wf/data/tpm_zscore.parquet')
    ptrap_scores = pd.read_parquet('data/ptrap_scores.parquet')
    ptrap_scores.index = ptrap_scores.index.droplevel(0)
    ptrap_genes = ptrap_scores.index

    # Pull out ptrap genes that we have scores for
    zscores.index = zscores.index.map(fbgn2symbol)
    zscores = zscores.reindex(ptrap_genes)

    # plot
    defaults = dict(yticklabels=False, xticklabels=False, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    if isinstance(expression_heatmap_kws, dict):
        defaults.update(expression_heatmap_kws)

    # make set set of defaults for protein scores plot
    defaults2 = dict(yticklabels=False, xticklabels=False, vmin=0, vmax=8, rasterized=True,
                     cmap='inferno', cbar_kws=dict(label='Arbitrary Score'))

    if isinstance(ptrap_heatmap_kws, dict):
        defaults2.update(ptrap_heatmap_kws)

    add_color_labels(axLabel, s=label_size)

    # Iterate over genes and build
    genes = [
        'osa',
        'Mapmodulin',
        'SRPK',
        'bol',
        'Fas3',
        'cindr',
    ]

    fig = plt.gcf()
    gs0 = GridSpecFromSubplotSpec(6, 1, subplot_spec=gsMain, hspace=0.2)
    for i, gene in zip(range(6), genes):
        gs00 = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[i, 0], hspace=0)
        ax1 = fig.add_subplot(gs00[0, 0])
        ax2 = fig.add_subplot(gs00[1, 0])

        sns.heatmap(zscores.loc[gene].to_frame().T, ax=ax1, **defaults)
        sns.heatmap(ptrap_scores.loc[gene].to_frame().T, ax=ax2, **defaults2)
        ax1.text(1.01, 0, gene, transform=ax1.transAxes, ha='left', va='center', fontsize=10)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 3, height_ratios=[1, 0.06], width_ratios=[0.1, 0.1, 1], hspace=0, wspace=0.2)

    gsMain = gs[0, 2]
    axLabel = fig.add_subplot(gs[1, 2])
    axCbar1 = fig.add_subplot(gs[0, 0])
    axCbar2 = fig.add_subplot(gs[0, 1])

    plot_heatmap_ptrap_genes(gsMain, axLabel, label_size=10,
                             expression_heatmap_kws=dict(cbar_ax=axCbar1),
                             ptrap_heatmap_kws=dict(cbar_ax=axCbar2),
                             )

    flip_ticks(axCbar1)
    flip_ticks(axCbar2)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
