import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import pickle_load

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]
COLORS = ["red", "C4", "C4", "gray", "gray", "cyan", "C2"]

def main():
    # Gene lists
    biased = pickle_load(snakemake.input.biased)

    # Data
    fbgn2chrom = pd.read_feather(snakemake.input.gene_annot, columns=['FBgn', 'FB_chrom']).set_index("FBgn").squeeze().rename("chrom")
    tpm = get_data(snakemake.wildcards.direction).join(fbgn2chrom)
    tpm_biased, tpm_not_biased = split_by_bias(tpm, biased)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    sns.boxplot("chrom", "TPM", data=tpm_biased, notch=True, order=CHROMS, palette=COLORS, ax=ax1)
    ax1.set(title=f"{snakemake.wildcards.direction}-Biased", xlabel="")
    sns.boxplot("chrom", "TPM", data=tpm_not_biased, notch=True, order=CHROMS, palette=COLORS, ax=ax2)
    ax2.set(title=f"{snakemake.wildcards.direction}-Expressed Not Biased", xlabel="")

    fig.savefig(snakemake.output[0])


def get_data(cluster):
    return (
        pd.read_feather(snakemake.input.tpm)
        .query(f"cluster == '{cluster}'")
        .drop("cluster", axis=1)
        .set_index("FBgn")
    )

def split_by_bias(tpm, fbgns):
    tpm_bias = tpm[tpm.index.isin(fbgns)]
    tpm_not_bias = tpm[~tpm.index.isin(fbgns)]
    return tpm_bias, tpm_not_bias



if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        TAG = config["tag"]

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                tpm="../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather",
                biased="../output/seurat3-cluster-wf/germline_deg/GvLPS_G_biased.pkl",
                gene_annot=f"../references/gene_annotation_dmel_{TAG}.feather",
            ),
            wildcards=dict(direction="G")
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
