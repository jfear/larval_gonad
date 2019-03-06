"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use('Agg')

from itertools import chain
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

fname = snakemake.input[0]
oname = snakemake.output[0]
annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors
cmap = snakemake.params.cmap


def main():
    df = get_data()

    plt.style.use('scripts/paper_1c.mplstyle')
    fig = plt.figure(figsize=(4, 8))
    gs = GridSpec(2, 2, height_ratios=[.1, 1], width_ratios=[1, .02], hspace=0, wspace=0.1)
    ax = fig.add_subplot(gs[:, 0])
    cax = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        df,
        xticklabels=True,
        yticklabels=False,
        vmin=-3,
        vmax=3,
        rasterized=True,
        cmap=cmap,
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(label='Z-Score (TPM)', ticks=[-3, 0, 3])
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')
    ax.set_xticklabels(list(chain.from_iterable([('', x, '') for x in cluster_order])), ha='center', va='bottom')

    # Add additional x annotations
    yloc = 0 - (df.shape[0] * .05)
    ax.text(6, yloc, 'Germline', ha='center', fontsize=8, color=cluster_colors[0], va='bottom')
    ax.text(17, yloc, 'Somatic Cyst', ha='center', fontsize=8, color=cluster_colors[4], va='bottom')
    ax.text(24, yloc, 'Somatic Other', ha='center', fontsize=8, color=cluster_colors[8], va='bottom')
    lines = [
        plt.Line2D([0, 12], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([12, 21], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([21, 24], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([24, 27], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax.add_line(l)

    # Clean up Y axis
    ax.set_ylabel('Genes')
    # Add lines separating cell types
    for i in range(1, 9):
        ax.axvline(i * 3, color='w', ls='--')

    # Clean up color bar
    plt.setp(cax.yaxis.label, fontsize=7)
    plt.setp(cax.get_yticklabels(), fontsize=7)

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    zscores = (
        pd.read_parquet(fname)
            .reset_index()
            .assign(
            cluster=lambda df: (
                df.cluster.map(annotation)
                    .pipe(lambda x: x[x != 'UNK'])
                    .astype('category')
                    .cat.as_ordered()
                    .cat.reorder_categories(cluster_order)
            )
        )
            .assign(
            rep=lambda df: (
                df.rep.astype('category')
                    .cat.as_ordered()
                    .cat.reorder_categories(['rep1', 'rep2', 'rep3'])
            )
        )
            .sort_values(by=['cluster', 'rep'])
            .pivot_table(index='FBgn', columns=['cluster', 'rep'], values='tpm_zscore')
    )

    # calculate linkages
    link = linkage(zscores.values, 'average')
    tree = dendrogram(link, no_plot=True)
    leaves = tree['leaves']

    return zscores.iloc[leaves, :]


if __name__ == '__main__':
    main()
