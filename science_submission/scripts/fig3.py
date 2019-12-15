""" Figure 3: Autosome ratios

"""
import os
from typing import List

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import shelve_load, pickle_load
from larval_gonad import plotting
from larval_gonad.plotting.expression import plot_log_expression_by_cluster
from larval_gonad.plotting.x_to_a import plot_x_to_a, plot_4_to_a


def main():
    plt.style.use(["2c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width * 2)

    _, axes = plt.subplots(4, 2, sharex=True, gridspec_kw=dict(hspace=0.05, wspace=0.05))
    all_genes(axes[:, 0])
    common_genes(axes[:, 1])

    fix_yaxis((0.1, 2.5), [0.5, 1, 1.5, 2, 2.5], axes[0, :])
    fix_yaxis((-0.1, 3), [0, 1, 2, 3], axes[1, :])
    fix_yaxis((-0.1, 8.2), [0, 4, 8], axes[2, :])
    fix_yaxis((-0.1, 10), [0, 5, 10], axes[3, :])

    axes[0, 0].set_title("All Genes", y=0.9)
    axes[0, 1].set_title("Common Genes", y=0.9)

    for (ax, label) in zip(axes.flat, ["A", "B", "C", "D", "E", "F", "G", "H"]):
        ax.text(0.01, 0.9, label, fontweight="bold", transform=ax.transAxes)

    plt.savefig(snakemake.output[0])


def all_genes(axes: List[plt.Axes]):
    plot_x_to_a(
        snakemake.input.all, snakemake.params.cluster_color, snakemake.params.cluster_order, axes[0]
    )

    plot_4_to_a(
        snakemake.input.all, snakemake.params.cluster_color, snakemake.params.cluster_order, axes[1]
    )

    plot_log_expression_by_cluster(
        snakemake.input.tpm,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        axes[2],
        Y_FBGNS,
    )
    axes[2].set(ylabel="Y Genes\nLog(TPM)")

    plot_log_expression_by_cluster(
        snakemake.input.tpm, snakemake.params.cluster_color, snakemake.params.cluster_order, axes[3]
    )


def common_genes(axes: List[plt.Axes]):
    plot_x_to_a(
        snakemake.input.common,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        axes[0],
    )

    plot_4_to_a(
        snakemake.input.common,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        axes[1],
    )

    plot_log_expression_by_cluster(
        snakemake.input.tpm,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        axes[2],
        set(Y_FBGNS).intersection(set(COMMON_FBGNS)),
    )
    axes[2].set(ylabel="Y Genes\nLog(TPM)")

    plot_log_expression_by_cluster(
        snakemake.input.tpm,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        axes[3],
        COMMON_FBGNS,
    )


def fix_yaxis(ylim, yticks, axes):
    for ax in axes:
        ax.set(ylim=ylim, yticks=yticks)
    axes[1].set(ylabel="", yticklabels=[])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        color_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                all="../output/x-to-a-wf/db/expressed.dat",
                common="../output/x-to-a-wf/db/commonly_expressed.dat",
                tpm="../output/seurat3-cluster-wf/tpm_by_cluster.feather",
                common_fbgns="../output/cellselection-wf/commonly_expressed_genes.pkl",
            ),
            params=dict(
                cluster_color=color_config["clusters"], cluster_order=config["cluster_order"]
            ),
        )

    Y_FBGNS = (
        pd.read_feather(snakemake.input.gene_annot)
        .set_index("FBgn")
        .query("FB_chrom == 'Y'")
        .index.tolist()
    )

    COMMON_FBGNS = pickle_load(snakemake.input.common_fbgns)

    main()
