import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

norm_file_name = snakemake.input.norm
clusters_file_name = snakemake.input.clusters
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(2.25, 1))

    sns.boxplot('cluster', 'prop_gene_on', data=df, palette=cluster_colors, ax=ax, notch=True, linewidth=.5,
                order=cluster_order)

    ax.set_axisbelow(True)
    ax.grid(axis='y')
    # sns.despine(ax=ax, left=True, bottom=True)

    # Clean up X axis
    ax.set_xlabel('')
    plt.setp(ax.get_xticklabels(), rotation=0)

    # Add additional x annotations
    yloc = -.03
    pad = yloc * .9
    ax.text(1.5, yloc + pad, 'Germline', ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax.text(5, yloc + pad, 'Somatic\nCyst', ha='center', va='center', fontsize=6, color=cluster_colors[4])
    ax.text(7.5, yloc + pad, 'Somatic\nOther', ha='center', va='center', fontsize=6, color=cluster_colors[8])
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
    ax.set_ylabel('Proportion of Genes On\n(Seurat Norm >= 1)')

    ax.set_title('Expression Diversity')
    fig.savefig(oname, bbox_inches='tight')


def get_data():
    norm = pd.read_parquet(norm_file_name)
    prop_on_per_cell = (
        (norm >= 1).mean()
        .rename('prop_gene_on')
        .to_frame()
        .join(pd.read_parquet(clusters_file_name))
    )

    return prop_on_per_cell


if __name__ == '__main__':
    main()
