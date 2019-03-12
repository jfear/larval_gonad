"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use('Agg')

from itertools import chain, combinations
from pickle import loads
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

fname = snakemake.input.zscores
biomarkers = snakemake.input.biomarkers

oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors
cmap = snakemake.params.cmap


def main():
    df, multi_genes_per_cluster = get_data()

    plt.style.use('scripts/paper_1c.mplstyle')
    fig = plt.figure(figsize=(1.8, 2.5))
    gs = GridSpec(2, 1, height_ratios=[1, .02], hspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
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
        cbar_kws=dict(label='Z-Score (TPM)', ticks=[-3, 0, 3], orientation='horizontal')
    )

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_tick_params(pad=0, length=2)
    ax.set_xticklabels(list(chain.from_iterable([('', x, '') for x in cluster_order])), ha='center', va='bottom',
                       fontsize=5.5)

    # Add additional x annotations
    yloc = 0 - (df.shape[0] * .10)
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

    # Add lines separating biomarker groups
    for loc in np.cumsum(multi_genes_per_cluster.values)[:-1]:
        ax.axhline(loc, color='w', ls='--', lw=.5)

    # Add additional y annotations
    loc = 0
    for label, cnt in multi_genes_per_cluster.iteritems():
        if label == 'G':
            _label = f'Germ Only\n(n={cnt:,})'
        elif label == 'S':
            _label = f'Soma Only\n(n={cnt:,})'
        else:
            _label = f'Germ and\nSoma\n(n={cnt:,})'
        ax.text(-.5, loc + (cnt / 2), _label, ha='right', va='center', fontsize=5.5, fontweight='bold')
        loc += cnt

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

    # Get biomarkers by cluster
    _biomarkers = (
        pd.read_csv(biomarkers, sep='\t', usecols=['primary_FBgn', 'cluster'], index_col=0)
        .cluster
        .map(annotation)
        .pipe(lambda x: x[x != "UNK"])
        .astype('category')
        .cat.as_ordered()
        .cat.reorder_categories(cluster_order)
        .rename_axis('FBgn')
        .to_frame()
    )

    # pull out genes that are in more than one cluster
    multi_fbgns = _biomarkers.groupby('FBgn').size().pipe(lambda x: x[x > 1]).index
    df = (
        _biomarkers.query(f'FBgn == {multi_fbgns.tolist()}')
            .groupby('FBgn')
            .apply(lambda df: '|'.join(df.cluster.sort_values().values))
    )

    germ_combos = ['|'.join(x) for x in chain(
        combinations(cluster_order[:4], 2),
        combinations(cluster_order[:4], 3),
        combinations(cluster_order[:4], 4),
    )]

    soma_combos = ['|'.join(x) for x in chain(
        combinations(cluster_order[4:], 2),
        combinations(cluster_order[4:], 3),
        combinations(cluster_order[4:], 4),
        combinations(cluster_order[4:], 5),
    )]

    flags = pd.Series(index=df.index).fillna('GS')
    flags[df.isin(germ_combos)] = 'G'
    flags[df.isin(soma_combos)] = 'S'
    flags = flags.astype('category').cat.as_ordered().cat.reorder_categories(['G', 'GS', 'S'])

    fbgns = flags.sort_values().index
    zscores = zscores.reindex(fbgns).dropna()

    # Count the number of genes in each group
    multi_genes_per_cluster = flags.value_counts().sort_index()

    return zscores, multi_genes_per_cluster


if __name__ == '__main__':
    main()
