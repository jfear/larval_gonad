import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_name = snakemake.input.dat
clusters_file_name = snakemake.input.clusters
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(2.25, 1))

    df.plot(kind='bar', stacked=True, width=.9, ax=ax)

    ax.set_axisbelow(True)
    ax.grid(axis='y')
    # sns.despine(ax=ax, left=True, bottom=True)

    # Clean up X axis
    ax.set_xlabel('')
    plt.setp(ax.get_xticklabels(), rotation=0)

    # Add additional x annotations
    yloc = -15
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
    ax.set_ylabel('Proportion of Cells')

    ax.set_title('Replicate Representation')
    ax.legend(loc='upper left', bbox_to_anchor=[1, 1])
    fig.savefig(oname, bbox_inches='tight')


def get_data():
    dat = (
        pd.read_parquet(file_name)
        .reset_index()
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(annotation), ordered=True, categories=cluster_order))
        .assign(rep=lambda df: pd.Categorical(df.rep, ordered=True, categories=['rep1', 'rep2', 'rep3']))
        .set_index(['cluster', 'rep'])
        .sort_index()
    )
    total = dat.groupby('cluster').number_of_cells.sum()

    plot_data = dat.unstack().div(total, axis='rows').fillna(0) * 100
    plot_data.columns = plot_data.columns.droplevel(0)
    return plot_data


if __name__ == '__main__':
    main()
