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

    plt.style.use('scripts/paper_1c.mplstyle')
    fig = plt.figure(figsize=(1.8, 10))
    gs = GridSpec(2, 1, height_ratios=[1, .01], hspace=0.01)
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
        cbar_kws=dict(label='Z-Score (TPM)', ticks=[-3, 0, 3], orientation='horizontal')
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_tick_params(pad=0, length=2)
    ax.set_xticklabels(list(chain.from_iterable([('', x, '') for x in cluster_order])), ha='center', va='bottom',
                       fontsize=5.5)

    # Add additional x annotations
    yloc = 0 - (df.shape[0] * .02)
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
    ax.yaxis.set_tick_params(pad=0.1, length=2)
    plt.setp(ax.get_yticklabels(), fontsize=7, fontstyle='italic', rotation=0, va='center')

    # Add lines separating lit genes
    for loc in [12, 34, 58, 67]:
        ax.axhline(loc, color='w', ls='--', lw=.5)

    # Add additional x annotations
    ax.text(27, 6, cluster_order[0], ha='left', va='center', fontsize=6, fontweight='bold')
    ax.text(27, 23, '\n'.join(cluster_order[1:4]), ha='left', va='center', fontsize=6, fontweight='bold')
    ax.text(27, 46, '\n'.join(cluster_order[4:7]), ha='left', va='center', fontsize=6, fontweight='bold')
    ax.text(27, 62.5, cluster_order[7], ha='left', va='center', fontsize=6, fontweight='bold')
    ax.text(27, 69, cluster_order[8], ha='left', va='center', fontsize=6, fontweight='bold')


    # Clean up color bar
    plt.setp(cax.xaxis.label, fontsize=6)
    plt.setp(cax.get_xticklabels(), fontsize=5)
    cax.xaxis.set_tick_params(pad=0, length=2)

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
