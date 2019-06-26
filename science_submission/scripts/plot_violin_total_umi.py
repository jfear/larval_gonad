import matplotlib

matplotlib.use('Agg')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import feather_to_cluster_rep_matrix

FNAME = snakemake.input[0]

CLUSTER_ANNOT = snakemake.params.cluster_annot
CLUSTER_ORDER = snakemake.params.cluster_order
CLUSTER_COLORS = snakemake.params.cluster_colors

ONAME = snakemake.output[0]

# Debug Settings
# FNAME = 'output/science_submission/raw_by_cluster_rep.feather'
# import yaml
# config = yaml.safe_load(open('config/common.yaml'))
# CLUSTER_ANNOT = config['cluster_annot']
# CLUSTER_ORDER = config['cluster_order']
# CLUSTER_COLORS = yaml.full_load(open('config/colors.yaml'))['clusters']


def main():
    df = (
        pd.read_feather(FNAME)
        .groupby(['FBgn', 'cluster'])
        .UMI.sum()
        .to_frame()
        .reset_index()
        .assign(log_UMI=lambda df: np.log10(df.UMI + 1))
    )

    plt.style.use('scripts/figure_styles.mplstyle')
    fig, ax = plt.subplots(figsize=(4, 3))

    sns.violinplot('cluster', 'log_UMI', data=df, inner=None, scale='count', order=CLUSTER_ORDER,
                    palette=CLUSTER_COLORS, linewidth=.5, zorder=0, width=1, ax=ax)

    sns.boxplot('cluster', 'log_UMI', data=df, order=CLUSTER_ORDER, width=.1, linewidth=.5, ax=ax,
                zorder=10, showfliers=False, boxprops=dict(zorder=10, facecolor='w'), whiskerprops=dict(zorder=10))

    ax.set_axisbelow(True)
    ax.grid(axis='y')
    ax.set(xlabel='', ylabel='Total UMI Per Cell\n(Log10)')
    plt.setp(ax.get_xticklabels(), rotation=0)

    fig.savefig(ONAME, bbox_inches='tight')


if __name__ == '__main__':
    main()
