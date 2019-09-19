"""Entire panel for the MSCI figure.

This is the entire panel for the MSCI figure.

* A) X/A ratio (commonly expressed)
* B) 4/A ratio (commonly expressed)
* C) Y genes Log(TPM)
* D) Commonly expressed genes Log(TPM)

"""
import os
import pandas as pd

import matplotlib.pyplot as plt

from larval_gonad.io import pickle_load
from larval_gonad.plotting.expression import plot_log_expression_by_cluster
from larval_gonad.plotting.x_to_a import plot_x_to_a, plot_4_to_a


def main():
    y_fbgns = (
        pd.read_feather(snakemake.input.gene_annot)
        .set_index("FBgn")
        .query("FB_chrom == 'Y'")
        .index.tolist()
    )
    commonly_expressed_fbgns = pickle_load(snakemake.input.commonly_expressed)

    plt.style.use(["1c", "science_base"])
    _, (ax1, ax2, ax3, ax4) = plt.subplots(
        4, 1, figsize=(2.5, 6), sharex=True, gridspec_kw=dict(hspace=0.02)
    )

    plot_x_to_a(
        snakemake.input.autosome_ratios,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        ax1,
    )

    plot_4_to_a(
        snakemake.input.autosome_ratios,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        ax2,
    )

    plot_log_expression_by_cluster(
        snakemake.input.tpm,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        ax3,
        y_fbgns,
    )
    ax3.set_ylabel("Y\nLog(TPM)")

    plot_log_expression_by_cluster(
        snakemake.input.tpm,
        snakemake.params.cluster_color,
        snakemake.params.cluster_order,
        ax4,
        commonly_expressed_fbgns,
    )
    ax4.set_ylabel("Commonly Expressed\nLog(TPM)")

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        color_config = read_config("config/colors.yaml")
        ASSEMBLY = config["assembly"]
        TAG = config["tag"]

        snakemake = snakemake_debug(
            workdir="science-submission",
            input=dict(
                autosome_ratios="output/x-to-a-wf/db/commonly_expressed.bak",
                tpm="output/seurat3-cluster-wf/tpm_by_cluster.feather",
                commonly_expressed="output/cellselection-wf/commonly_expressed_genes.pkl",
                gene_annot=f"references/gene_annotation_{ASSEMBLY}_{TAG}.feather",
            ),
            params=dict(
                cluster_color=color_config["clusters"], cluster_order=config["cluster_order"]
            ),
        )

    main()
