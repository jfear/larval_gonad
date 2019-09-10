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

    # If FBgns provided then subset dataset
    fbgns = snakemake.params.get("fbgns", False)
    if fbgns:
        df = df.loc[df.FBgn.isin(fbgns), :]

    ax = sns.boxplot(
        "cluster",
        "UMI",
        data=df,
        palette=snakemake.params.cluster_colors,
        notch=True,
        linewidth=0.5,
        showfliers=False,
    )
    title = snakemake.params.get("title", "")
    ax.set(xlabel="", ylabel="Total UMI Per Cell (Log10)", title=title, yscale="symlog")

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        colors_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input="../output/seurat3-cluster-wf/raw_by_cluster_rep.feather",
            params=dict(
                cluster_colors=colors_config["clusters"],
                # fbgns=["FBgn0002962", "FBgn0003525"]
            ),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
