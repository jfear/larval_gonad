"""Plot tSNE

Plots the tsne generated by `scrnaseq-wf/scripts/scrnaseq_combine_force.Rmd`.

"""
import matplotlib as mpl

mpl.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

clusters = snakemake.input.clusters
tsne = snakemake.input.tsne
oname = snakemake.output[0]
colors = snakemake.params.colors
cluster_order = snakemake.params.cluster_order
cluster_annot = snakemake.params.annotation


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, axes = plt.subplots(1, 4, figsize=(2, 6))

    for i, clus in enumerate(cluster_order[:4]):
        plot_tsne(df, clus, colors[i], axes[i])

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    return (
        pd.read_parquet(tsne)
        .join(pd.read_parquet(clusters))
        .assign(color=lambda df: pd.Categorical(df.cluster, ordered=True, categories=cluster_order + ['UNK', ]))
    )


def plot_tsne(df, title, color, ax):
    defaults = dict(edgecolor='k', lw=.02, s=3, rasterized=True)

    _df = df.query(f'cluster == "{title}"')
    ax.scatter(df.tSNE_1, df.tSNE_2, c='lightgray', alpha=.1, **defaults)
    ax.scatter(_df.tSNE_1, _df.tSNE_2, c=color, **defaults)
    sns.despine(ax=ax, left=True, bottom=True)
    plt.setp(ax, xticks=[], yticks=[], xlabel='', ylabel='', aspect='equal', xmargin=0, ymargin=0)
    ax.set_title(title)
    return ax


if __name__ == '__main__':
    main()
