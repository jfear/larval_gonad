"""Plot barchart Spearman correlation with bulk RNA-Seq data."""

import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

from larval_gonad.config import read_config, config
from larval_gonad.plotting import cluster_cmap

mpl.style.use('scripts/paper_1c.mplstyle')

config.update(read_config('config.yaml'))


def plot_barchart_bulk_corr(ax, **kwargs):
    sc = pd.read_parquet('../scrnaseq-wf/data/tpm.parquet')

    bulk = pd.read_parquet('../bulk-rnaseq-wf/data/aggregation/tpm_gene_level_counts.parquet',
                           columns=['C1_TDT', 'C2_TDT', 'C3_TDT', 'C4_TDT'])
    bulk = bulk.mean(axis=1)
    bulk.name = 'bulk'

    # I want a vector of correlation coeff with bulk.
    corr = sc.join(bulk).corr(method='spearman')
    dat = corr.bulk.drop('bulk')

    # Make colors for plotting
    colors = [cluster_cmap[x] for x in dat.index]

    # plot
    dat.plot.bar(color=colors, ax=ax)
    ax.set_axisbelow(True)
    ax.grid(linestyle=':', axis='y')


if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(8, 4))

    plot_barchart_bulk_corr(ax)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlabel('Cell Cluster')
    ax.set_ylabel('Correlation with Bulk RNA-Seq (Spearman)')

    fig.savefig(snakemake.output[0], bbox_inches='tight')
