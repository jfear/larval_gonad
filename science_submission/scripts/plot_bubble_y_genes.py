import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

fname = snakemake.input.raw
fbgn2chrom = snakemake.input.fbgn2chrom
fbgn2symbol = snakemake.input.fbgn2symbol
clusters = snakemake.input.clusters
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors
cmap = snakemake.params.cmap


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(1, 1, figsize=(3, 6))

    # Get X, Y coordinates using category codes.
    xvals = df.cluster.cat.codes
    xlabels = df.cluster.cat.categories

    yvals = df.gene_symbol.cat.codes
    ylabels = df.gene_symbol.cat.categories

    # Main plot
    sc = ax.scatter(xvals, yvals, cmap='Purples_r', c=df['sum'], s=df['prop_cells'] * 100, vmin=1, vmax=1000,
                    edgecolor='k', lw=.1)
    # ax.set_axisbelow(True)
    # ax.grid(axis='y')
    sns.despine(ax=ax, left=True, bottom=True)
    ax.margins(0.02)

    # Add color bar and legend
    plt.colorbar(sc, orientation='horizontal', ticks=[1, 500, 1000], label='Total Number of Reads', fraction=.05,
                 pad=.12)

    for sizes in [5, 10, 50]:
        plt.scatter([], [], c='k', alpha=0.3, s=sizes, label=f'{sizes:,.0f} %')

    legend = plt.legend(scatterpoints=1, frameon=False, labelspacing=.6, title='Percent of Cells', loc='upper left',
                        bbox_to_anchor=[1, 1], bbox_transform=ax.transAxes)
    plt.setp(legend.get_title(), fontsize=8)

    # Clean up X
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)

    # Add additional x annotations
    yloc = -6.5
    pad = yloc * .4
    ax.text(1.5, yloc - pad, 'Germline', ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax.text(5, yloc - pad, 'Somatic\nCyst', ha='center', va='center', fontsize=6, color=cluster_colors[4])
    ax.text(7.5, yloc - pad, 'Somatic\nOther', ha='center', va='center', fontsize=6, color=cluster_colors[8])
    ax.text(10, yloc + pad, 'Unknown', ha='center', fontsize=6, color=cluster_colors[-1], va='center')
    lines = [
        plt.Line2D([0, 3], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([4, 6], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([7, 7.5], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([7.5, 8], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
        plt.Line2D([9, 11], [yloc, yloc], color=cluster_colors[-1], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax.add_line(l)

    # Clean up Y axis
    ax.set_ylabel('')
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontstyle='italic')

    fig.savefig(oname, bbox_inches='tight')


def prop_cells(x):
    return (x > 0).mean()


def get_data():
    _fbgn2symbol = pd.read_pickle(fbgn2symbol)

    df = (
        pd.read_parquet(fname)
        # pull out Y-linked genes and get their gene symbol
        .join(pd.read_parquet(fbgn2chrom))
        .query('chrom == "chrY"')
        .drop('chrom', axis=1)
        .assign(gene_symbol=lambda df: df.index.map(_fbgn2symbol))
        .set_index('gene_symbol')
        # Merge on cluster and aggregate UMI: sum and prop of cells with expression
        .T
        .join(pd.read_parquet(clusters))
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(annotation), ordered=True, categories=cluster_order))
        .rename_axis('cell_id')
        .reset_index()
        .melt(id_vars=['cell_id', 'cluster'], var_name='gene_symbol', value_name='UMI')
        .groupby(['cluster', 'gene_symbol']).agg({'UMI': ['sum', prop_cells]})
    )

    # clean up
    df.columns = df.columns.droplevel(0)
    df = df.reset_index()

    # Make gene_symbol a ordered category for easier plotting
    df.gene_symbol = pd.Categorical(df.gene_symbol, ordered=True, categories=reversed(
        sorted(df.gene_symbol.unique(), key=lambda x: x.lower())
    ))

    return df.dropna()


if __name__ == '__main__':
    main()
