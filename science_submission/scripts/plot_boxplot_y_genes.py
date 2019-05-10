import matplotlib

matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

fname = snakemake.input.tpm
fbgn2chrom = snakemake.input.fbgn2chrom
clusters = snakemake.input.clusters
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(1, 1, figsize=(2.25, 1))

    sns.boxplot('cluster', 'TPM', data=df, order=cluster_order, palette=cluster_colors, ax=ax, linewidth=.5)

    # Clean up X
    ax.set_xlabel('')

    # Add additional x annotations
    yloc = -150
    pad = yloc * .4
    ax.text(1.5, yloc - pad, 'Germline', ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax.text(5, yloc - pad, 'Somatic\nCyst', ha='center', va='center', fontsize=6, color=cluster_colors[4])
    ax.text(7.5, yloc - pad, 'Somatic\nOther', ha='center', va='center', fontsize=6, color=cluster_colors[8])
    ax.text(10, yloc - pad, 'Unknown', ha='center', fontsize=6, color=cluster_colors[-1], va='center')
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
    ax.set_ylabel('Normalized Count (TPM)')

    fig.savefig(oname, bbox_inches='tight')


def prop_cells(x):
    return (x > 0).mean()


def get_data():
    return (
        pd.read_parquet(fname)
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(annotation), ordered=True, categories=cluster_order))
        .dropna()
        .join(pd.read_parquet(fbgn2chrom))
        .query('chrom == "chrY"')
        .groupby(['rep', 'cluster']).TPM.sum()
        .reset_index()
    )


if __name__ == '__main__':
    main()
