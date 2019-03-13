"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt

fname = snakemake.input[0]
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(2, 2))
    df.plot(kind='bar', stacked=True, ax=ax, legend=False, width=.9)
    plt.legend(loc='upper left', bbox_to_anchor=[1, 1])

    # Clean up X axis
    ax.set_xlabel('')
    plt.setp(ax.get_xticklabels(), rotation=0)

    # Add additional x annotations
    yloc = 0 - (df.shape[0] * .025)
    pad = yloc * .3
    ax.text(1.5, yloc - pad, 'Germline', ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax.text(5, yloc - pad, 'Somatic\nCyst', ha='center', va='center', fontsize=6, color=cluster_colors[4])
    ax.text(7.5, yloc - pad, 'Somatic\nOther', ha='center', va='center', fontsize=6, color=cluster_colors[8])
    lines = [
        plt.Line2D([0, 3], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([4, 6], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([7, 7.5], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([7.5, 8], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax.add_line(l)

    # Clean up Y axis
    ax.set_ylabel('Proportion of Cells')

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    return (
        pd.read_parquet(fname)
        .assign(cluster=lambda df: df.cluster.map(annotation))
        .query('cluster != "UNK"')
        .assign(rep=lambda df: df.index.str.extract('(rep\d)', expand=False))
        .groupby(['cluster', 'rep'])
        .size()
        .unstack()
        # Calculate the proportion of cells from each rep
        .pipe(lambda df: df.div(df.sum(axis=1), axis='rows'))
        .loc[cluster_order, :]
    )


if __name__ == '__main__':
    main()
