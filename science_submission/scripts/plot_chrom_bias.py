import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq
from larval_gonad.plotting.stats import add_pvals

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]


def main():
    bg = get_background()
    fbgn2chrom = get_fbgn2chrom(bg)
    expressed_fbgns = get_expressed_fbgns()
    counts_per_chrom = calculate_counts_per_chrom(fbgn2chrom, expressed_fbgns)
    chisq = get_short_chisq_results(counts_per_chrom)
    tidy = make_chisq_tidy(chisq)

    ax = sns.barplot("chrom", "count", hue="type", data=tidy)
    add_qvalue(chisq, ax)
    ax.set(title=snakemake.wildcards.clus, ylabel="Number Expressed Genes", ylim=(0, 1200))
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])

    plt.savefig(snakemake.output[0])


def get_background():
    return (
        pd.read_feather(snakemake.input.tpm)
        .groupby("FBgn")
        .TPM.sum()
        .pipe(lambda x: x[x > 0])
        .index.unique().tolist()
    )


def get_fbgn2chrom(bg):
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .rename("chrom")
        .pipe(lambda x: x[x.isin(CHROMS)])
        .reindex(bg).dropna()
    )


def get_expressed_fbgns():
    return (
        pd.read_feather(snakemake.input.tpm)
        .query(f"cluster == '{snakemake.wildcards.clus}'")
        .groupby("FBgn")
        .TPM.mean()
        .pipe(lambda x: x[x >= 10])
        .index.tolist()
    )


def calculate_counts_per_chrom(fbgn2chrom, expressed_fbgns):
    total = fbgn2chrom.value_counts().rename("total")
    expressed = fbgn2chrom.reindex(expressed_fbgns).dropna().value_counts().rename("expressed")
    return pd.concat([expressed, total], axis=1, sort=True).reindex(CHROMS).T.fillna(0)


def get_short_chisq_results(cnts):
    df = run_chisq(cnts).loc[("expressed", ["observed", "expected", "fdr q-value"]), :].T
    df.columns = df.columns.droplevel(0)
    df.columns = df.columns.tolist()
    df.columns.name = ""
    df.index.name = "chrom"
    return df


def make_chisq_tidy(chisq):
    return (
        chisq[["observed", "expected"]]
        .reset_index()
        .melt(id_vars="chrom", var_name="type", value_name="count")
    )


def add_qvalue(chisq, ax):
    return add_pvals(
        x=range(len(CHROMS)),
        y=chisq[["observed", "expected"]].max(axis=1),
        pval=chisq["fdr q-value"],
        ax=ax
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
            wildcards=dict(clus="LPS"),
        )

    plt.style.use("../config/figure_styles.mplstyle")
    plt.rcParams.update({"figure.figsize": (4, 2)})

    main()
