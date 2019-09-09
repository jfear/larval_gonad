import os
import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns


def main():
    norm = (
        pd.read_feather(snakemake.input.norm)
        .set_index("FBgn")
        .drop("gene_symbol", axis=1)
        .loc[snakemake.wildcards.FBgn]
        .rename("norm")
    )
    clusters = pd.read_feather(snakemake.input.clusters).drop("rep", axis=1).set_index("cell_id")
    df = clusters.join(norm)

    ax = sns.pointplot("cluster", "norm", data=df)
    sns.despine(ax=ax)
    ax.set_title(snakemake.wildcards.symbol, fontstyle="italic")
    ax.set(xlabel="", ylabel="Normalized Expression (by cell)")

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                norm="../output/seurat3-cluster-wf/combined_n3_normalized.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
            ),
            wildcards=dict(symbol="bol", FBgn="FBgn0011206"),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
