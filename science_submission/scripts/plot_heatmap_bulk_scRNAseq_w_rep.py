"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""

import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

from larval_gonad.config import config
from larval_gonad.plotting import add_color_labels, flip_ticks, cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_heatmap_bulk_scRNAseq(axMain, **kwargs):
    # Bulk data
    tdt = ['C1_TDT', 'C3_TDT', 'C2_TDT', 'C4_TDT']
    bulk = pd.read_parquet('../output/bulk-rnaseq-wf/aggregation/tpm_gene_level_counts.parquet', columns=tdt)
    bulk.columns = ['Bulk-rep1', 'Bulk-rep2', 'Bulk-rep3', 'Bulk-rep4']

    # import scRNA-seq data
    sc = pd.read_parquet('../output/scrnaseq-wf/tpm_w_rep.parquet', columns=config['sel_cluster_order_w_rep'])

    # Merge
    dat = sc.join(bulk, how='inner')
    _corr = dat.corr()

    # calculate linkages
    link = linkage(_corr.values, 'average')
    tree = dendrogram(link, no_plot=True)
    leaves = tree['leaves']

    # plot
    defaults = dict(yticklabels=False, xticklabels=False, vmin=-3, vmax=3, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)

    sns.heatmap(_corr.iloc[leaves, leaves], ax=axMain, **defaults)
    axMain.set_ylabel('')


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(1, 2, width_ratios=[0.1, 1], wspace=0.2)
    axMain = fig.add_subplot(gs[0, 1])
    axCbar = fig.add_subplot(gs[0, 0])
    plot_heatmap_bulk_scRNAseq(axMain, cbar_ax=axCbar)

    flip_ticks(axCbar)
    flip_ticks(axMain, pos='right')
    plt.setp(axMain.get_yticklabels(), rotation=0)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
