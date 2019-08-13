import os
import shelve

import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.plotting import format_pval


def main(snake):
    with shelve.open(os.path.splitext(snake["input_file"])[0]) as db:
        df = db["data"].pipe(lambda x: x[x.ratio_type == snake["ratio_type"]])
        pvals = db["pvalues"].pipe(lambda x: x[x.ratio_type == snake["ratio_type"]])

    ylabel = YLABELS[snake["ratio_type"]]

    ax = sns.boxplot(
        x="cluster",
        y="ratio",
        data=df,
        palette=snake["cluster_colors"],
        order=snake["cluster_order"],
    )
    ax.axhline(1)
    ax.set(xlabel="Cluster ", ylabel=ylabel, title=snake["title"], ylim=(0, 3))
    for _, row in pvals.iterrows():
        format_pval(ax, row.x, row.y, row.pvalue)

    plt.savefig(snake["output_file"])


if __name__ == "__main__":
    SNAKE = dict(
        input_file=snakemake.input[0],
        output_file=snakemake.output[0],
        cluster_colors=snakemake.params["cluster_colors"],
        cluster_order=snakemake.params["cluster_order"],
        ratio_type=snakemake.wildcards.ratio_type,
        title=snakemake.params.title,
    )

    plt.style.use("scripts/figure_styles.mplstyle")

    TXT = "Normalized by Number of Genes"
    YLABELS = {
        "x_to_a_ratio": f"X / A\n{TXT}",
        "fourth_to_a_ratio": f"4 / A\n{TXT}",
        "y_to_a_ratio": f"Y / A\n{TXT}",
    }

    # Debug Settings
    # try:
    #     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
    #     print(os.getcwd())
    # except:
    #     pass
    # from larval_gonad.config import read_config
    # config = read_config('../config/config.yaml')
    # common_config = read_config('../../config/common.yaml')
    # color_config = read_config('../../config/colors.yaml')
    # snake = dict(
    #     input_file='../../output/x-to-a-wf/db/expressed.dat',
    #     output_file='',
    #     cluster_colors=color_config['clusters'],
    #     cluster_order=common_config['cluster_order'],
    #     ratio_type='fourth_to_a_ratio',
    #     title=config['gene_subset_titles']['expressed']
    # )
    # plt.style.use("figure_styles.mplstyle")

    main(SNAKE)
