"""Plot barchart showing X2A ratio.

"""

import pandas as pd

import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

from larval_gonad.config import read_config, config
from larval_gonad.plotting import add_color_labels, cluster_cmap
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')

config.update(read_config('config.yaml'))


def plot_barchart_germ_x2a_ratio(ax, **kwargs):
    tpm = pd.read_parquet('../output/scrnaseq-wf/tpm.parquet').join(fbgn2chrom)
    tpm.chrom = tpm.chrom.map(dict(chr2L='A', chr2R='A', chr3L='A', chr3R='A', chrX='X', chr4='4', chrY='Y'))

    # Calculate the median expression
    meds = tpm.groupby('chrom').median()

    # Calculate x2a ratio
    ratio = meds.T['X'] / meds.T['A']
    ratio = ratio.loc[config['cluster_order']].copy()

    # Make colors for plotting
    colors = [cluster_cmap[x] for x in ratio.index]

    # plot
    ratio.plot.bar(color=colors, ax=ax)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':', axis='y')


if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(8, 4))

    plot_barchart_germ_x2a_ratio(ax)
    ax.set_xticklabels([])
    ax.set_xlabel('Cell Cluster')
    ax.set_ylabel('X to Autosome Ratio (of Medians)')

    fig.savefig(snakemake.output[0], bbox_inches='tight')
