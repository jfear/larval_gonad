import matplotlib

matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

fname = snakemake.input[0]
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors


def main():
    df = get_data()

    plt.style.use('scripts/paper_1c.mplstyle')
    fig, ax = plt.subplots(figsize=(2.25, 1))

    sns.violinplot('cluster', 'log_UMI', data=df, inner=None, scale='count', order=cluster_order,
                   palette=cluster_colors, linewidth=.5, zorder=0, width=1, ax=ax)

    sns.boxplot('cluster', 'log_UMI', data=df, order=cluster_order, width=.1, linewidth=.5, showfliers=False, ax=ax,
                zorder=10, boxprops=dict(zorder=10, facecolor='w'), whiskerprops=dict(zorder=10))

    ax.set_axisbelow(True)
    ax.grid(axis='y', linewidth=.5)
    sns.despine(ax=ax, left=True, bottom=True)

    # Clean up X axis
    ax.set_xlabel('')
    ax.xaxis.set_tick_params(pad=.8, length=0)
    ax.xaxis.set_ticks_position('top')
    plt.setp(ax.get_xticklabels(), rotation=0, fontsize=5.5)

    # Add additional x annotations
    yloc = 6.5
    pad = yloc * .07
    ax.text(1.5, yloc + pad, 'Germline', ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax.text(5, yloc + pad, 'Somatic\nCyst', ha='center', va='center', fontsize=6, color=cluster_colors[4])
    ax.text(7.5, yloc + pad, 'Somatic\nOther', ha='center', va='center', fontsize=6, color=cluster_colors[8])
    lines = [
        plt.Line2D([0, 3], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([4, 6], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([7, 7.5], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([7.5, 8], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax.add_line(l)

    # Clean up Y axis
    ax.set_ylabel('Total UMI Per Cell\n(Log10)', fontsize=6)
    ax.yaxis.set_tick_params(pad=0.2, length=0)
    plt.setp(ax.get_yticklabels(), fontsize=5)

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    return (
        pd.read_parquet(fname)
            .sum()
            .rename('UMI')
            .to_frame()
            .assign(rep=lambda df: df.index.str.extract('(rep\d)', expand=False))
            .join(pd.read_parquet('../output/scrnaseq-wf/clusters.parquet'))
            .assign(cluster=lambda df: df.cluster.map(annotation))
            .query('cluster != "UNK"')
            .assign(log_UMI=lambda df: np.log10(df.UMI))
    )


if __name__ == '__main__':
    main()
