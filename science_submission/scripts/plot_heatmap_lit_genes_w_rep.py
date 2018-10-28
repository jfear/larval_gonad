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


def plot_heatmap_lit_genes_w_rep(axMain, **kwargs):
    zscores = pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet', columns=config['sel_cluster_order_w_rep'])
    lit_genes = config['lit_genes']

    # Pull out lit genes
    zscores.index = zscores.index.map(fbgn2symbol)
    zscores = zscores.reindex(lit_genes)

    # plot
    defaults = dict(yticklabels=True, xticklabels=True, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)
    sns.heatmap(zscores, ax=axMain, **defaults)
    axMain.set_ylabel('')


if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(4, 4))
    plot_heatmap_lit_genes_w_rep(ax, cbar=False)

    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.setp(ax.get_xticklabels(), rotation=50, ha='left', fontsize=6)

    plt.setp(ax.get_yticklabels(), rotation=0, fontstyle='italic')

    fig.savefig(snakemake.output[0], bbox_inches='tight')
