import os
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    df = (
        pd.read_feather(snakemake.input[0])
        .set_index("cell_id")
        .groupby(["cluster", "rep"])
        .size()
        .unstack()
    )

    df_prop = df.div(df.sum(axis=1), axis="rows")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=plt.figaspect(2), sharex=True)
    df.plot.bar(stacked=True, ax=ax1)
    ax1.set(xlabel="", ylabel="Number of Cells")
    df_prop.plot.bar(stacked=True, ax=ax2)
    ax2.legend_.remove()
    ax2.set(xlabel="", ylabel="Proportion of Cells", ylim=(0, 1))

    # Rotate tick labels
    plt.setp(ax2.get_xticklabels(), rotation=0)

    # Add annotations
    rep_totals = make_text(df.sum(axis=0))
    cluster_totals = make_text(df.sum(axis=1))

    defaults = dict(va="top", ha="left", transform=ax1.transAxes, fontsize=8)
    ax1.text(0, 1, rep_totals, **defaults)
    ax1.text(0, 0.8, cluster_totals, **defaults)

    fig.savefig(snakemake.output[0])


def make_text(data):
    return re.sub(" +", ": ", data.map(lambda x: f"{x:,}").rename_axis("").to_string().strip())


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
