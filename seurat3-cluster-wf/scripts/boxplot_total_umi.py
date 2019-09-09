import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    df = (
        pd.read_feather(snakemake.input[0])
        .groupby(["FBgn", "cluster"])
        .UMI.sum()
        .to_frame()
        .reset_index()
    )

    ax = sns.boxplot(
        "cluster",
        "UMI",
        data=df,
        palette=snakemake.params.cluster_colors,
        notch=True,
        linewidth=0.5,
        showfliers=False,
    )
    ax.set(xlabel="", ylabel="Total UMI Per Cell (Log10)", yscale="symlog")

    fig.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        colors_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input="../output/seurat3-cluster-wf/raw_by_cluster_rep.feather",
            params=dict(cluster_colors=colors_config["clusters"]),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
