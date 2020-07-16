""" Figure 3: Boxplots of Total Expression, X/AA and 44/AA, and Y Expression"""
from string import ascii_uppercase
from typing import List

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting.stats import add_pvals

plt.style.use(["2c", "science_base"])


def main():
    df = (
        load_genes()
        .pipe(add_gene_set_annotation)
        .pipe(add_total_expression)
        .pipe(add_y_expression)
        .pipe(add_ratios, snakemake.input.expressed_ratios, "All Genes")
        .pipe(
            add_ratios,
            snakemake.input.widely_expressed_ratios,
            "Widely Expressed Genes",
        )
        .pipe(tidy)
    )

    g = sns.FacetGrid(
        df,
        row="row",
        col="gene_set",
        col_order=["All Genes", "Widely Expressed Genes"],
        sharey="row",
        sharex=True,
        aspect=1.2,
        gridspec_kws=dict(wspace=0.15),
    )
    g.map(
        sns.boxplot,
        "cluster",
        "value",
        order=snakemake.params.cluster_order,
        palette=snakemake.params.cluster_color,
        notch=True,
        linewidth=0.5,
        showfliers=False,
    )
    tweak_axes(g)
    add_x_stats(g.axes[1])
    add_4th_stats(g.axes[2])
    add_y_stats(g.axes[3])

    g.savefig(snakemake.output[0])


def load_genes() -> pd.DataFrame:
    return pd.read_feather(
        snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"]
    ).rename(columns={"FB_chrom": "chrom"})


def add_gene_set_annotation(df: pd.DataFrame) -> pd.DataFrame:
    expressed = joblib.load(snakemake.input.expressed_fbgns)
    widely_expressed = joblib.load(snakemake.input.widely_expressed_fbgns)
    return pd.concat(
        [
            (df.query("FBgn in @expressed").assign(gene_set="All Genes")),
            (
                df.query("FBgn in @widely_expressed").assign(
                    gene_set="Widely Expressed Genes"
                )
            ),
        ],
        ignore_index=True,
        sort=False,
    )


def _load_log_tpm() -> pd.DataFrame:
    return np.log10(
        pd.read_feather(snakemake.input.tpm).set_index(["FBgn", "cluster"]).squeeze()
        + 1
    ).reset_index()


def add_total_expression(df: pd.DataFrame) -> pd.DataFrame:
    tpm = _load_log_tpm().rename(columns={"TPM": "Total Expression\nLog(TPM)"})
    return df.merge(tpm, on="FBgn", how="outer")


def add_y_expression(df: pd.DataFrame) -> pd.DataFrame:
    y_genes = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .query("FB_chrom == 'Y'")
        .FBgn.to_list()
    )

    tpm = (
        _load_log_tpm()
        .query("FBgn in @y_genes")
        .rename(columns={"TPM": "Y Expression\nLog(TPM)"})
    )

    return df.merge(tpm, on=["FBgn", "cluster"], how="outer")


def add_ratios(df: pd.DataFrame, gene_set_db, gene_set_name) -> pd.DataFrame:
    ratios = (
        shelve_load(gene_set_db)["data"]
        .query("ratio_type != 'y_to_a_ratio'")
        .reset_index()
        .drop("rep", axis=1)
        .set_index(["cell_id", "cluster", "ratio_type"])
        .squeeze()
        .unstack()
        .reset_index()
        .rename(columns={"x_to_a_ratio": "X/AA", "fourth_to_a_ratio": "44/AA"})
        .assign(gene_set=gene_set_name)
    )
    return pd.concat([df, ratios], ignore_index=True, sort=False)


def tidy(df: pd.DataFrame) -> pd.DataFrame:
    return pd.melt(
        df,
        id_vars=["gene_set", "cluster"],
        value_vars=[
            "Total Expression\nLog(TPM)",
            "X/AA",
            "44/AA",
            "Y Expression\nLog(TPM)",
        ],
        var_name="row",
        value_name="value",
    ).dropna()


def tweak_axes(g: sns.FacetGrid):
    g.set_titles("{col_name}")
    g.set_xlabels("")

    # Make Row Name the y-axis Label
    for row_name, row in zip(g.row_names, g.axes):
        row[0].set_ylabel(row_name)

    for label, ax in zip(ascii_uppercase, g.axes.ravel()):
        ax.text(
            0.01,
            1,
            label,
            ha="left",
            va="top",
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=12,
        )


def _add_horizontal_line(axes: List[plt.Axes], loc: float):
    for ax in axes:
        ax.axhline(loc, color="gray", ls="--")


def _add_ratio_pvalues(shelve: str, ratio_type: str, ax: plt.Axes):
    pvals = shelve_load(shelve)["pvalues"].query("ratio_type == @ratio_type")
    add_pvals(pvals.x, pvals.y, pvals.pvalue, ax, fontsize=12, fontweight="bold")


def add_x_stats(axes: List[plt.Axes]):
    _add_ratio_pvalues(snakemake.input.expressed_ratios, "x_to_a_ratio", axes[0])
    _add_ratio_pvalues(snakemake.input.widely_expressed_ratios, "x_to_a_ratio", axes[1])
    _add_horizontal_line(axes, 1)


def add_4th_stats(axes: List[plt.Axes]):
    _add_ratio_pvalues(snakemake.input.expressed_ratios, "fourth_to_a_ratio", axes[0])
    _add_ratio_pvalues(
        snakemake.input.widely_expressed_ratios, "fourth_to_a_ratio", axes[1]
    )
    _add_horizontal_line(axes, 1)


def add_y_stats(axes: List[plt.Axes]):
    # From: notebook/2020-01-03-y_stats.ipynb
    add_pvals([1], [1.45], [0.0017], axes[0], fontsize=12, fontweight="bold")


if __name__ == "__main__":
    main()
