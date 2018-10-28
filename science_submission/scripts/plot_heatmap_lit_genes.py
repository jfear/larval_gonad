"""Plot heatmap of literature genes.

Plots the tpm normalized zscores of gene from the literature. The current gene
list can be found in `science_submission/config.yaml`.
"""

import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import read_config, config
from larval_gonad.plotting import add_color_labels, flip_ticks
from common import fbgn2symbol

mpl.style.use('scripts/paper_1c.mplstyle')

config.update(read_config('config.yaml'))


def plot_heatmap_lit_genes(axMain, axLabel, label_size=5, **kwargs):
    zscores = pd.read_parquet('../output/scrnaseq-wf/tpm_zscore.parquet', columns=config['sel_cluster_order'])
    lit_genes = config['lit_genes']

    # Pull out lit genes
    zscores.index = zscores.index.map(fbgn2symbol)
    zscores = zscores.reindex(lit_genes)

    # plot
    defaults = dict(yticklabels=True, xticklabels=False, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)

    sns.heatmap(zscores, ax=axMain, **defaults)
    add_color_labels(axLabel, s=label_size)

    axMain.set_ylabel('')


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 2, height_ratios=[1, 0.06], width_ratios=[0.1, 1], hspace=0, wspace=0.2)
    axMain = fig.add_subplot(gs[0, 1])
    axLabel = fig.add_subplot(gs[1, 1])
    axCbar = fig.add_subplot(gs[0, 0])
    plot_heatmap_lit_genes(axMain, axLabel, label_size=10, cbar_ax=axCbar)

    flip_ticks(axCbar)
    flip_ticks(axMain, pos='right')
    plt.setp(axMain.get_yticklabels(), rotation=0, fontstyle='italic')

    fig.savefig(snakemake.output[0], bbox_inches='tight')
