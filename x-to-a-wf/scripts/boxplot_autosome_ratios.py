import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting.stats import add_pvals


TXT = "Normalized by Number of Genes"
YLABELS = {
    "x_to_a_ratio": f"X / A\n{TXT}",
    "fourth_to_a_ratio": f"4 / A\n{TXT}",
    "y_to_a_ratio": f"Y / A\n{TXT}",
}


def main():
    db = shelve_load(snakemake.input[0])
    df = db["data"].pipe(lambda x: x[x.ratio_type == snakemake.wildcards.ratio_type])
    pvals = db["pvalues"].pipe(lambda x: x[x.ratio_type == snakemake.wildcards.ratio_type])

    ylabel = YLABELS[snakemake.wildcards.ratio_type]

    ax = sns.boxplot(
        x="cluster",
        y="ratio",
        data=df,
        palette=snakemake.params.cluster_color,
        order=snakemake.params.cluster_order,
        notch=True
    )
    ax.set(xlabel="Cluster ", ylabel=ylabel, title=snakemake.params.title)
    add_pvals(pvals.x, pvals.y, pvals.pvalue, ax)

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        color_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="x-to-a-wf",
            input="../output/x-to-a-wf/db/ovary.bak",
            params=dict(
                cluster_color=color_config["clusters"],
                cluster_order=config["cluster_order"],
                title="Test",
            ),
            wildcards=dict(fbgns="ovary", ratio_type="x_to_a_ratio"),
        )

    plt.style.use("../config/figure_styles.mplstyle")
    plt.rcParams.update({
        'figure.figsize': (4, 2)
    })

    main()
