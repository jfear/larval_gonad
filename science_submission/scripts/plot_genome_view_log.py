"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.config import config

mpl.style.use('scripts/paper_1c.mplstyle')


def plot_genome_view_log(ax):
    germ_cells = config['cluster_order'][:4]
    locs = pd.read_parquet('../output/science_submission/gene_locations.parquet')
    tpm = np.log10(pd.read_parquet('../output/scrnaseq-wf/tpm.parquet', columns=germ_cells) + 1)

    # munge data
    dat = locs.join(tpm)
    dat['color'] = dat[germ_cells].idxmax(axis=1)
    dat['value'] = dat[germ_cells].max(axis=1)
    dat.color = pd.Categorical(dat.color, categories=germ_cells, ordered=True)
    dat = dat[['location', 'value', 'color']].copy()

    cmapper = dict(zip(germ_cells, sns.color_palette('Set1')))

    for g, dd in dat.groupby('color'):
        markerline, stemlines, baseline = ax.stem(dd.location, dd.value, label=g, c=cmapper[g])
        plt.setp(markerline, markerfacecolor=cmapper[g], markeredgecolor=cmapper[g], markersize=2, alpha=.5)
        plt.setp(stemlines, color=cmapper[g], lw=1, alpha=.5)
        plt.setp(baseline, color='k', lw=1)

    ax.set_xlim(0, dat.location.max())
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))


if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize=(8, 3))
    plot_genome_view_log(ax)
    fig.savefig(snakemake.output[0], bbox_inches='tight')
