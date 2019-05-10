import matplotlib

matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

fname = snakemake.input['ratios']
pvals = snakemake.input['pvals']
oname = snakemake.output[0]

annotation = snakemake.params.annotation
cluster_order = snakemake.params.cluster_order
cluster_colors = snakemake.params.cluster_colors


def main():
    df = get_data()

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(2.25, 2), sharex=True, sharey=False, gridspec_kw=dict(hspace=0.1))

    sns.boxplot('cluster', 'x_to_a_ratio', data=df, palette=cluster_colors, ax=ax1,
                notch=True, linewidth=.5)
    ax1.set_axisbelow(True)
    ax1.axhline(1, ls='--', lw=.5, color='#b0b0b0', zorder=0)

    sns.boxplot('cluster', 'fourth_to_a_ratio', data=df, palette=cluster_colors, ax=ax2,
                notch=True, linewidth=.5)
    ax2.set_axisbelow(True)
    ax2.axhline(1, ls='--', lw=.5, color='#b0b0b0', zorder=0)

    # Clean up X axis
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    # plt.setp(ax2.get_xticklabels(), rotation=0, fontsize=5.5)

    # Add additional x annotations
    yloc = -1.2
    pad = yloc * .3
    ax2.text(1.5, yloc - pad, 'Germline', ha='center', va='center', fontsize=6, color=cluster_colors[0])
    ax2.text(5, yloc - pad, 'Somatic\nCyst', ha='center', va='center', fontsize=6, color=cluster_colors[4])
    ax2.text(7.5, yloc - pad, 'Somatic\nOther', ha='center', va='center', fontsize=6, color=cluster_colors[8])
    ax2.text(10, yloc - pad, 'Unknown', ha='center', fontsize=6, color=cluster_colors[-1], va='center')
    lines = [
        plt.Line2D([0, 3], [yloc, yloc], color=cluster_colors[0], lw=1.5, clip_on=False),
        plt.Line2D([4, 6], [yloc, yloc], color=cluster_colors[4], lw=1.5, clip_on=False),
        plt.Line2D([7, 7.5], [yloc, yloc], color=cluster_colors[7], lw=1.5, clip_on=False),
        plt.Line2D([7.5, 8], [yloc, yloc], color=cluster_colors[8], lw=1.5, clip_on=False),
        plt.Line2D([9, 11], [yloc, yloc], color=cluster_colors[-1], lw=1.5, clip_on=False),
    ]

    for l in lines:
        ax2.add_line(l)

    # Clean up Y axis
    ax1.set_ylabel('X:A')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax2.set_ylabel('4:A')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


    # Add p-value indicators
    whiskers = whisker_locations(df)
    df_pvals = get_pvals()
    for i, (cluster, dd) in enumerate(df_pvals.iterrows()):
        if dd.pval_x <= 0.05:
            pval_str = pval_to_string(dd.pval_x)
            loc = whiskers.loc[cluster, 'x_to_a_ratio']
            ax1.text(i, loc + .01, pval_str, ha='center', va='bottom', fontsize=6)

        if dd.pval_4 <= 0.05:
            pval_str = pval_to_string(dd.pval_4)
            loc = whiskers.loc[cluster, 'fourth_to_a_ratio']
            ax2.text(i, loc + .01, pval_str, ha='center', va='bottom', fontsize=6)

    fig.savefig(oname, bbox_inches='tight')


def get_data():
    return (
        pd.read_parquet(fname)
            .query('cluster != [9, 10, 11]')
            .assign(
            cluster=lambda df: pd.Categorical(df.cluster.map(annotation), ordered=True, categories=cluster_order))
    )

def pval_to_string(pval):
    if pval <= 0.001:
        return '***'
    elif pval <= 0.01:
        return '**'
    elif pval <= 0.05:
        return '*'
    else:
        return ''

def whisker_locations(df):
    """Calculate the upper whisker location for plotting astrix."""
    return (
        df.groupby('cluster').quantile([.25, .75])
            .rename_axis(['cluster', 'quantile'])
            .reset_index()
            .melt(id_vars=['cluster', 'quantile'], var_name='ratio', value_name='value')
            .pivot_table(index=['cluster', 'ratio'], columns='quantile', values='value')
            .assign(iqr=lambda df: df[0.75] - df[0.25])
            .assign(upper=lambda df: df[0.75] + (1.5 * df['iqr']))
            .upper.unstack()
    )


def get_pvals():
    """P-values from permutation test."""
    return (
        pd.read_parquet(pvals)
            .reindex(range(9))
            .rename(index=annotation)
            .reindex(cluster_order)
    )


if __name__ == '__main__':
    main()
