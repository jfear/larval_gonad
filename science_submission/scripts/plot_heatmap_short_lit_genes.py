"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use('Agg')

from itertools import chain
from pickle import loads
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

fname = snakemake.input.zscores
fbgn2symbol = snakemake.input.fbgn2symbol

oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors
lit_genes = snakemake.params.lit_genes
cmap = snakemake.params.cmap


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig = plt.figure(figsize=(1.8, 1.5))
    gs = GridSpec(2, 1, height_ratios=[1, .05])
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        df,
        xticklabels=True,
        yticklabels=True,
        vmin=-3,
        vmax=3,
        rasterized=True,
        cmap=cmap,
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(label='Z-Score (TPM)', orientation='horizontal', ticks=[-3, 0, 3])
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')
    ax.set_xticklabels(list(chain.from_iterable([('', x, '') for x in cluster_order])), ha='center', va='bottom')

    # Add additional x annotations
    yloc = 0 - (df.shape[0] * .12)
    pad = yloc * .1
    ax.text(6, yloc + pad, 'Germline', ha='center', fontsize=6, color=cluster_colors[0], va='bottom')
    ax.text(17, yloc + pad, 'Somatic\nCyst', ha='center', fontsize=6, color=cluster_colors[4], va='bottom')
    ax.text(24, yloc + pad, 'Somatic\nOther', ha='center', fontsize=6, color=cluster_colors[8], va='bottom')
    lines = [
        plt.Line2D([0, 12], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([12, 21], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([21, 24], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([24, 27], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax.add_line(l)

    # Add lines separating cell types
    for i in range(1, 9):
        ax.axvline(i * 3, color='w', ls='--', lw=.5)

    # Clean up Y axis
    ax.set_ylabel('')
    plt.setp(ax.get_yticklabels(), fontsize=7, fontstyle='italic', rotation=0, va='center', ha='right')

    # Add lines separating lit genes
    for loc in [2, 4, 6, 8]:
        ax.axhline(loc, color='w', ls='--', lw=.5)

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

    with open(fbgn2symbol, 'rb') as fh:
        mapper = loads(fh.read())

    zscores.index = zscores.index.map(mapper)
    return zscores.reindex(lit_genes)


if __name__ == '__main__':
    main()
