import os

import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import pickle_load

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]


def main():
    # Get gene lists for expressed and deg
    cluster = snakemake.wildcards.direction
    expressed_fbgns = get_expressed_fbgns(cluster)
    deg_fbgns = pickle_load(snakemake.input.deg)

    # Get gene starts for different subsets
    gene_starts = get_gene_starts()
    expressed_starts = (
        gene_starts.reindex(expressed_fbgns).dropna().sort_values(by=["chrom", "start"])
    )
    deg_starts = gene_starts.reindex(deg_fbgns).dropna().sort_values(by=["chrom", "start"])

    # Plot starts for different subsets
    fig, axes = plt.subplots(7, 1, sharex=True, sharey=True)
    for ax, chrom in zip(axes.flat, CHROMS):
        defaults = dict(height=1, rasterized=True, ax=ax)
        # All genes
        sns.rugplot(gene_starts.query(f"chrom == '{chrom}'").start, color="lightgray", **defaults)
        # Expressed genes
        sns.rugplot(expressed_starts.query(f"chrom == '{chrom}'").start, color="black", **defaults)
        # DEG genes
        sns.rugplot(deg_starts.query(f"chrom == '{chrom}'").start, color="red", **defaults)

        # Clean up axes
        sns.despine(ax=ax, left=True, bottom=True)
        ax.set_ylabel(chrom, rotation=0, va="center")
        ax.set(yticks=[])

    plt.suptitle(f"{cluster}-Biased")

    fig.savefig(snakemake.output[0])


def get_expressed_fbgns(cluster):
    return (
        pd.read_feather(snakemake.input.tpm)
        .query(f"cluster == '{cluster}' & TPM >= 10")
        .FBgn.unique()
        .tolist()
    )


def get_gene_starts():
    return (
        pd.read_feather(snakemake.input.gene_annot)
        .set_index("FBgn")
        .assign(chrom=lambda x: x.FB_chrom)
        .loc[:, ["chrom", "start"]]
        .query(f"chrom == {CHROMS}")
        .sort_values(by=["chrom", "start"])
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        TAG = config["tag"]

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                tpm="../output/seurat3-cluster-wf/raw_by_cluster_rep.feather",
                deg="../output/seurat3-cluster-wf/germline_deg/GvLPS_G_biased.pkl",
                gene_annot=f"../references/gene_annotation_dmel_{TAG}.feather",
            ),
            wildcards=dict(direction="G"),
        )

    plt.style.use("../config/figure_styles.mplstyle")
    plt.rcParams.update({"figure.figsize": (4, 2)})

    main()
