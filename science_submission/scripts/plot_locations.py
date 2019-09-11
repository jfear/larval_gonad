import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import pickle_load

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]


def main():
    deg_fbgns = pickle_load(snakemake.input.deg)
    df = (
        pd.read_feather(snakemake.input.gene_annot)
        .set_index("FBgn")
        .assign(chrom=lambda x: x.FB_chrom)
        .assign(coordinate=lambda x: x.start)
        .loc[:, ["chrom", "coordinate"]]
        .query(f"chrom == {CHROMS}")
        .sort_values(by=["chrom", "coordinate"])
    )
    subset = df.reindex(deg_fbgns).dropna().sort_values(by=["chrom", "coordinate"])

    fig, axes = plt.subplots(7, 1, sharex=True, sharey=True)
    for ax, chrom in zip(axes.flat, CHROMS):
        sns.rugplot(df.query(f"chrom == '{chrom}'").coordinate, height=1, color="lightgray", rasterized=True, ax=ax)
        sns.rugplot(subset.query(f"chrom == '{chrom}'").coordinate, height=1, color="red", rasterized=True, ax=ax)
        sns.despine(ax=ax, left=True, bottom=True)
        ax.set_ylabel(chrom, rotation=0)
        ax.set(yticks=[])

    plt.suptitle(snakemake.wildcards.fbgns)

    fig.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        TAG = config["tag"]

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                deg="../output/seurat3-cluster-wf/combined_n3_gonia_biased.pkl",
                gene_annot=f"../references/gene_annotation_dmel_{TAG}.feather",
            ),
            wildcards=dict(fbgns="gonia"),
        )

    plt.style.use("../config/figure_styles.mplstyle")
    plt.rcParams.update({"figure.figsize": (4, 2)})

    main()
