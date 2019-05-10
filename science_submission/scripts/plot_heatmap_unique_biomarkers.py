"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use('Agg')

from itertools import chain
from pickle import loads
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

fname = snakemake.input.zscores
fbgn2symbol = snakemake.input.fbgn2symbol
biomarkers = snakemake.input.biomarkers

oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors
cmap = snakemake.params.cmap


def main():
    df, unique_genes_per_cluster = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig = plt.figure(figsize=(3, 4))
    gs = GridSpec(2, 1, height_ratios=[1, .02])
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        df,
        xticklabels=True,
        yticklabels=False,
        vmin=-3,
        vmax=3,
        rasterized=False,
        cmap=cmap,
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(label='Z-Score (TPM)', ticks=[-3, 0, 3], orientation='horizontal')
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')
    ax.set_xticklabels(list(chain.from_iterable([('', x, '') for x in cluster_order])), ha='center', va='bottom')

    # Add additional x annotations
    yloc = 0 - (df.shape[0] * .05)
    pad = yloc * .1
    ax.text(6, yloc + pad, 'Germline', ha='center', fontsize=6, color=cluster_colors[0], va='bottom')
    ax.text(17, yloc + pad, 'Somatic\nCyst', ha='center', fontsize=6, color=cluster_colors[4], va='bottom')
    ax.text(24, yloc + pad, 'Somatic\nOther', ha='center', fontsize=6, color=cluster_colors[8], va='bottom')
    ax.text(31, yloc + pad, 'Unknown', ha='center', fontsize=6, color=cluster_colors[-1], va='bottom')
    lines = [
        plt.Line2D([0, 12], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([12, 21], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([21, 24], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([24, 27], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
        plt.Line2D([27, 35], [yloc, yloc], color=cluster_colors[-1], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax.add_line(l)

    # Add lines separating cell types
    for i in range(1, 12):
        ax.axvline(i * 3, color='w', ls='--', lw=.5)

    # Clean up Y axis
    ax.set_ylabel('')

    # Add lines separating biomarker groups
    for loc in np.cumsum(unique_genes_per_cluster)[:-1]:
        ax.axhline(loc, color='w', ls='--', lw=.5)

    # Add additional y annotations
    loc = 0
    for i, cnt in enumerate(unique_genes_per_cluster):
        if i == 10:
            x = 38
        else:
            x = 35
        ax.text(-.5, loc + (cnt / 2), cluster_order[i], ha='right', va='center', fontweight='bold', fontsize=5.5)
        ax.text(x, loc + (cnt / 2), f'n={cnt}', ha='left', va='center', fontweight='bold', fontsize=5.5)
        loc += cnt

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    zscores = (
        pd.read_parquet(fname)
        .reset_index()
        .assign(
            cluster=lambda df: (
                df.cluster.map(annotation)
                # .pipe(lambda x: x[x != 'UNK'])
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

    # Get biomarkers by cluster
    _biomarkers = (
        pd.read_csv(biomarkers, sep='\t', usecols=['primary_FBgn', 'cluster'], index_col=0)
        .cluster
        .map(annotation)
        # .pipe(lambda x: x[x != "UNK"])
        .astype('category')
        .cat.as_ordered()
        .cat.reorder_categories(cluster_order)
        .rename_axis('FBgn')
        .to_frame()
    )

    # pull out genes that were unique to a single cluster
    unique_fbgns = _biomarkers.groupby('FBgn').size().pipe(lambda x: x[x == 1]).index
    df = _biomarkers.query(f'FBgn == {unique_fbgns.tolist()}').sort_values(by='cluster')
    fbgns = df.index
    zscores = zscores.reindex(fbgns)

    # Count the number of genes in each group
    unique_genes_per_cluster = df.cluster.value_counts().sort_index().tolist()

    # Convert to gene names
    with open(fbgn2symbol, 'rb') as fh:
        mapper = loads(fh.read())
    zscores.index = zscores.index.map(mapper)

    return zscores, unique_genes_per_cluster


if __name__ == '__main__':
    main()
