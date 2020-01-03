""" Figure WFK: boxplot version """
import os
from typing import List
from collections import ChainMap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

from larval_gonad.normalization import tpm
from larval_gonad import plotting

CHROMS = ["X", "2L", "2R", "3L", "3R", "4", "Y"]
PLOT_DEFAULTS = dict(x="chrom", order=CHROMS, notch=True)


def main():
    plt.style.use(["2c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width * 0.5)

    _, axes = plt.subplots(
        1, 2, sharex=True, sharey=True, gridspec_kw=dict(hspace=0.05, wspace=0.05)
    )
    plot_read_counts(axes)

    for ax, label in zip(axes, ["A", "B"]):
        sns.despine(ax=ax)

    plt.savefig(snakemake.output[0])


def plot_read_counts(axes: List[plt.Axes]):
    kwargs = ChainMap({}, PLOT_DEFAULTS, dict(y="TPM"))
    df = pd.concat([
        read_adult_bulk(),
        read_l3_bulk()
    ], ignore_index=True)
    plot_adult_bulk(axes[0], df, **kwargs)
    plot_L3_bulk(axes[1], df, **kwargs)

    for ax in axes[1:]:
        ax.set(ylabel="")


def plot_adult_bulk(ax: plt.Axes, df: pd.DataFrame, **kwargs):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'adult'")
    sns.boxplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **kwargs,
    )
    ax.set_xlabel("")
    ax.set_title("Adult Bulk", y=0.9)
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


def plot_L3_bulk(ax: plt.Axes, df: pd.DataFrame, **kwargs):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'larval'")
    sns.boxplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **kwargs,
    )
    ax.set_xlabel("")
    ax.set_title("L3 Bulk", y=0.9)
    ax.get_legend().remove()
    return ax


def read_adult_bulk() -> pd.DataFrame:
    gene_lens = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )
    dat = pd.pivot_table(
        pd.read_feather(snakemake.input.adult), index="FBgn", columns="sample_ID", values="Count"
    )
    _tpm = tpm(dat, gene_lens).dropna()
    return (
        _tpm.join(get_fbgn2chrom())
        .reset_index()
        .melt(id_vars=["FBgn", "chrom"], var_name="sample_ID", value_name="TPM")
        .assign(stage="adult")
        .assign(tissue=lambda x: x.sample_ID.str.extract(r".*_(\w+)_r\d", expand=False))
    )


def read_l3_bulk() -> pd.DataFrame:
    mapper = dict(TCP="testis", OCP="ovary")
    df = pd.read_csv(snakemake.input.larval, sep="\t", index_col=0).rename_axis("FBgn")
    cols = [f"l3_bulk_{x}" for x in df.columns]
    df.columns = cols
    return (
        df.join(get_fbgn2chrom(), how="inner")
        .reset_index()
        .melt(id_vars=["FBgn", "chrom"], var_name="sample_ID", value_name="TPM")
        .assign(stage="larval")
        .assign(tissue=lambda x: x.sample_ID.str.extract(r".*_(TCP|OCP)", expand=False))
        .assign(tissue=lambda x: x.tissue.map(mapper))
    )


def get_fbgn2chrom():
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .rename("chrom")
    )



def adult_stats():
    df = read_adult_bulk()
    for chrom, dd in df.groupby("chrom"):
        testis = dd.query("tissue == 'testis'")
        ovary = dd.query("tissue == 'ovary'")
        _, pval = mannwhitneyu(testis.TPM, ovary.TPM, alternative="less")
        print(chrom, pval)

    for chrom, dd in df.groupby("chrom"):
        testis = dd.query("tissue == 'testis'")
        ovary = dd.query("tissue == 'ovary'")
        _, pval = mannwhitneyu(testis.TPM, ovary.TPM, alternative="greater")
        print(chrom, pval)


def larval_stats():
    df = read_l3_bulk()
    for chrom, dd in df.groupby("chrom"):
        if chrom not in CHROMS:
            continue
        testis = dd.query("tissue == 'testis'")
        ovary = dd.query("tissue == 'ovary'")
        _, pval = mannwhitneyu(testis.TPM, ovary.TPM, alternative="less")
        print(chrom, pval)

    for chrom, dd in df.groupby("chrom"):
        if chrom not in CHROMS:
            continue
        testis = dd.query("tissue == 'testis'")
        ovary = dd.query("tissue == 'ovary'")
        _, pval = mannwhitneyu(testis.TPM, ovary.TPM, alternative="greater")
        print(chrom, pval)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        COLORS = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                adult="../output/expression-atlas-wf/w1118_gene_counts.feather",
                larval="../output/bulk-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv",
            ),
            params=dict(colors=COLORS),
        )

    main()
