import os
from itertools import chain
from more_itertools import flatten

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

from larval_gonad.config import read_config


def main():
    lit_genes = list(flatten(read_config(snakemake.input.lit_genes).values()))
    biomarkers = get_unique_biomarkers()
    zscore = pd.pivot_table(
        pd.read_feather(snakemake.input.zscore),
        values="zscore",
        index="FBgn",
        columns=["cluster", "rep"],
        aggfunc="first",
    ).reindex(biomarkers.index)

    fig = plt.figure(figsize=plt.figaspect(2))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        zscore,
        xticklabels=True,
        yticklabels=False,
        vmin=-3,
        vmax=3,
        rasterized=True,
        cmap=snakemake.params.color,
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(label="Z-Score (TPM)", ticks=[-3, 0, 3], orientation="horizontal"),
    )
    cleanup_xaxis(ax)
    cleanup_yaxis(ax, biomarkers)
    add_annotations(ax, biomarkers, lit_genes, zscore.shape[1])

    plt.savefig(snakemake.output[0])


def get_unique_biomarkers():
    biomarkers = (
        pd.read_feather(snakemake.input.biomarkers)
        .sort_values("FBgn")
        .assign(
            cluster=lambda df: pd.Categorical(
                df.cluster.map(dict(snakemake.params.cluster_annot)),
                categories=snakemake.params.cluster_order,
                ordered=True,
            )
        )
    )
    return (
        biomarkers[~biomarkers.duplicated(subset="FBgn", keep=False)]
        .sort_values("cluster")
        .set_index("FBgn")
    )


def cleanup_xaxis(ax):
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        list(chain.from_iterable([("", x, "") for x in snakemake.params.cluster_order])),
        ha="center",
        va="bottom",
    )

    # Add lines separating cell types
    for i in range(1, len(snakemake.params.cluster_order)):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    return ax


def cleanup_yaxis(ax, biomarkers):
    ax.set_ylabel("")

    # Add lines separating biomarker groups
    loc = 0
    for _, dd in biomarkers.groupby("cluster"):
        loc += dd.shape[0]
        ax.axhline(loc, color="w", ls="--", lw=0.5)

    return ax


def add_annotations(ax, biomarkers, lit_genes, cols=30):
    """Adds annotation on the left and right sides.

    Annotation includes the cluster name, the number of genes, and the gene
    names of literature genes found in that set.

    """

    loc = 0
    xloc_odd, xloc_even = (-5, cols + cols * 0.01)
    for i, (clus, dd) in enumerate(biomarkers.groupby("cluster")):
        loc, mid = take_step(loc, dd.shape[0])
        txt = f"{clus} ({dd.shape[0]:,})"
        txt += check_for_lit_genes(dd, lit_genes)

        if i % 2 == 0:
            xloc = xloc_even
        else:
            xloc = xloc_odd

        ax.text(xloc, mid, txt, ha="left", va="center", fontsize=8)


def take_step(loc, step_size):
    prev = loc
    loc += step_size
    mid = loc - ((loc - prev) / 2)
    return loc, mid


def check_for_lit_genes(df, lit_genes):
    res = []
    for fbgn, dd in df.sort_values("gene_symbol").iterrows():
        if fbgn in lit_genes:
            res.append(dd.gene_symbol)
    return "\n" + "\n".join(res)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                biomarkers="../output/seurat3-cluster-wf/combined_n3_biomarkers.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
                zscore="../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
                lit_genes="../config/literature_genes.yaml",
            ),
            params=dict(
                color="viridis",
                cluster_annot=config["cluster_annot"],
                cluster_order=config["cluster_order"],
            ),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
