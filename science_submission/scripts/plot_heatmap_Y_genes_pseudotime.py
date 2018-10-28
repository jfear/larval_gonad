"""Plot heatmap of genes on Chromosome Y.

Plots the tpm normalized zscores of gene from the literature. The current gene
list can be found in `science_submission/config.yaml`.
"""

import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import read_config, config
from larval_gonad.plotting import add_triangle, flip_ticks
from larval_gonad.normalization import tpm, zscore
from common import fbgn2chrom, fbgn2symbol

mpl.style.use('scripts/paper_1c.mplstyle')

config.update(read_config('config.yaml'))


def plot_heatmap_Y_genes_pseudotime(axMain, axLabel, label_size=5, **kwargs):
    raw = pd.read_parquet('../output/science_submission/raw_germcells.parquet')
    gene_lengths = pd.read_csv('../output/gene_ts_lengths.tsv', sep='\t', index_col=0).gene_ts_length
    _tpm = tpm(raw, gene_lengths).dropna()
    zscores = zscore(_tpm).dropna()

    # Pull out Y genes
    zscoresY = zscores.join(fbgn2chrom).query('chrom == "chrY"').drop('chrom', axis=1)
    zscoresY.index = zscoresY.index.map(fbgn2symbol)

    # calculate linkages
    link = linkage(zscoresY.values, 'average')
    tree = dendrogram(link, no_plot=True)
    leaves = tree['leaves']

    # plot
    defaults = dict(yticklabels=True, xticklabels=False, vmin=-.52, vmax=.75, rasterized=True,
                    cmap=config['colors']['heatmap'], cbar_kws=dict(label='Normalized Read Counts (Z-score)'))

    defaults.update(kwargs)

    sns.heatmap(zscoresY.iloc[leaves], ax=axMain, **defaults)
    add_triangle(axLabel)

    axMain.set_ylabel('')


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 8))
    gs = GridSpec(2, 2, height_ratios=[1, 0.06], width_ratios=[0.1, 1], hspace=0, wspace=0.2)
    axMain = fig.add_subplot(gs[0, 1])
    axLabel = fig.add_subplot(gs[1, 1])
    axCbar = fig.add_subplot(gs[0, 0])
    plot_heatmap_Y_genes_pseudotime(axMain, axLabel, label_size=10, cbar_ax=axCbar)

    flip_ticks(axCbar)
    flip_ticks(axMain, 'right')
    plt.setp(axMain.get_yticklabels(), rotation=0)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
