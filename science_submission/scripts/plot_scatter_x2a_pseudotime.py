"""Plot heatmap of genes on Chromosome X.

Plots the tpm normalized zscores of gene from the literature. The current gene
list can be found in `science_submission/config.yaml`.
"""

import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import read_config, config
from larval_gonad.plotting import add_triangle
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')

config.update(read_config('config.yaml'))


def plot_scatter_x2a_pseudotime(axMain, axLabel, **kwargs):
    # convert autosomes to 'A'
    fbgn2chrom.chrom.replace(dict(chrX='X', chr2L='A', chr2R='A', chr3L='A', chr3R='A'), inplace=True)

    # Figure out how many reads are on X and A, scale for the number of genes.
    raw = pd.read_parquet('data/raw_germcells.parquet').join(fbgn2chrom)
    totals = raw.groupby('chrom').sum().div(fbgn2chrom.groupby('chrom').size(), axis=0).T

    # Calculate X/A ratio
    ratio = totals['X'] / totals['A']
    ratio.name = 'X2A'

    # add on pseudotime for plotting
    psTime = pd.read_parquet('data/pseudotime.parquet')
    dat = psTime.join(ratio)

    # Plot
    state_cmap = dict(zip(range(1, 6), sns.color_palette('icefire_r', n_colors=5)))
    for g, dd in dat.groupby('State'):
        dd.plot.scatter('Pseudotime', 'X2A', ax=axMain, color=state_cmap[g], **kwargs)

    axMain.set_xlabel('')
    axMain.set_ylabel('X/A (Scaled Total Read Counts)')
    axMain.set_xticklabels([])
    axMain.set_xticks([])

    add_triangle(axLabel)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(2, 1, height_ratios=[1, 0.06], hspace=0)
    axMain = fig.add_subplot(gs[0, 0])
    axLabel = fig.add_subplot(gs[1, 0])
    plot_scatter_x2a_pseudotime(axMain, axLabel)

    fig.savefig(snakemake.output[0], bbox_inches='tight')
