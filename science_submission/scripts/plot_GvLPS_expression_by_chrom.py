import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]

def main():
    expressed = np.log1p(get_fbgns_expressed_in_G_and_LPS())
    fbgn2chrom = get_fbgn2chrom()
    tidy = make_tidy(expressed, fbgn2chrom)
    ma = make_ma_data(expressed, fbgn2chrom)

    g = sns.FacetGrid(ma, col="chrom", col_order=CHROMS, hue="direction", col_wrap=4)
    g.map(plt.scatter, "M", "A", s=8)

    fig, axes = plt.subplots(len(CHROMS), 1, figsize=(8, 12))
    for ax, (chrom, dd) in zip(axes.flat, tidy.groupby("chrom")):
        g = dd.query("cluster == 'G'").TPM.sort_values()
        lps = dd.query("cluster == 'LPS'").TPM.sort_values()
        ax.step(g, np.arange(g.shape[0]), label="G")
        ax.step(lps, np.arange(lps.shape[0]), label="LPS")
        ax.legend()
        ax.set_title(chrom)


def get_fbgns_expressed_in_G_and_LPS():
    df = (
        pd.read_feather(snakemake.input.tpm)
        .groupby(["FBgn", "cluster"])
        .TPM.mean()
        .unstack()
        .loc[:, ["G", "LPS"]]
    )
    df.columns = df.columns.tolist()
    return df
    return df[(df >= 10).all(axis=1)]


def get_fbgn2chrom():
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .rename("chrom")
        .pipe(lambda x: x[x.isin(CHROMS)])
    )


def make_tidy(expressed, fbgn2chrom):
    return (
        expressed
        .reset_index()
        .melt(id_vars="FBgn", var_name="cluster", value_name="TPM")
        .set_index("FBgn")
        .sort_index()
        .join(fbgn2chrom)
    )


def make_ma_data(expressed, fbgn2chrom):
    return (
        expressed
        .assign(M=lambda x: x.mean(axis=1))
        .assign(A=lambda x: x.G - x.LPS)
        .assign(direction=lambda x: np.where(x.A > 0, "G", "LPS"))
        .loc[:, ["M", "A", "direction"]]
        .join(fbgn2chrom)
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
                tpm="../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather",
                gene_annot=f"../references/gene_annotation_dmel_{TAG}.feather",
            ),
        )

    main()
