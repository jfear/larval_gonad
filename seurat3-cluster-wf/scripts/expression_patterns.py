import os
import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import safe_gene_name


def main():
    norm = pd.read_feather(snakemake.input.norm).set_index(["FBgn", "gene_symbol"])
    clusters = pd.read_feather(snakemake.input.clusters).drop("rep", axis=1).set_index("cell_id")

    for idx, dd in norm.iterrows():
        df = pd.concat([dd.rename("norm"), clusters], axis=1)
        ax = sns.pointplot("cluster", "norm", data=df)
        sns.despine(ax=ax)
        ax.set_title(idx[1], fontstyle="italic")
        ax.set(xlabel="", ylabel="Normalized Expression (by cell)")

        plt.savefig(snakemake.params.pattern.format(FBgn=idx[0], symbol=safe_gene_name(idx[1])))
        plt.close()

    Path(snakemake.output[0]).touch()

if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                norm="../output/seurat3-cluster-wf/combined_n3_normalized.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
            ),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
