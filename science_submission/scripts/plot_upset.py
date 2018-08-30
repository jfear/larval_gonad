"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import matplotlib.pyplot as plt
import upsetplot as ups

from larval_gonad.config import config
from common import fbgn2chrom

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_upset(gs):
    germ_cells = config['cluster_order'][:4]
    tpm = pd.read_parquet('../scrnaseq-wf/data/tpm.parquet', columns=germ_cells)
    xlink = fbgn2chrom.query('chrom == "chrX"').index

    # munge data
    dat = tpm > 2
    dat = dat.reindex(xlink).dropna()
    dat = dat[dat.sum(axis=1) > 0]
    dat = dat.groupby(germ_cells[::-1]).size().sort_values(ascending=False)

    return ups.plot(dat, fig, subplot_spec=gs, sort_by='cardinality', sort_sets_by=None)


if __name__ == '__main__':
    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(1, 1)
    plot_upset(gs[0, 0])
    fig.savefig(snakemake.output[0], bbox_inches='tight')
