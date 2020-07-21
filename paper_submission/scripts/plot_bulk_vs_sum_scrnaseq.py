"""Compare bulk and scRNA-seq data."""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import joblib
from scipy.stats import spearmanr

from larval_gonad.normalization import tpm
import larval_gonad.plotting  # pylint: disable=unused-import


plt.style.use("minimal")


def main():
    df = pd.concat([get_avg_bulk(), get_avg_sc()], axis=1, sort=False, join="inner")
    male_biased_fbgns = joblib.load(snakemake.input.male_biased)  # pylint: disable=unused-variable
    df_male_biased = df.query("FBgn in @male_biased_fbgns")

    ax: plt.Axes
    fig, ax = plt.subplots()
    x, y = "log_avg_tpm_sum_sc", "log_avg_tpm_bulk"
    ax.hexbin(df[x], df[y], cmap="RdBu_r", mincnt=1)
    add_reg_plot(x, y, df, "All Genes", ax, color="k")
    add_reg_plot(x, y, df_male_biased, "Male-Biased Genes", ax, color="k", ls="--")
    xlim = df[x].min(), df[x].max()
    ylim = df[y].min(), df[y].max()
    tweak_plot(ax, xlim, ylim)
    ax.set_xlabel(r"$Log_{10}\left(\overline{TPM(\sum{Single\ Cell})}\right)$")
    ax.set_ylabel(r"$Log_{10}\left(\overline{TPM(Bulk)}\right)$")
    sns.despine(ax=ax)

    fig.savefig(snakemake.output[0])


def get_avg_bulk() -> pd.DataFrame:
    """Calculate the Log10 Average TPM for Bulk RNA-Seq"""
    gene_lengths = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    counts = (
        pd.read_csv(snakemake.input.bulk, sep="\t")
        .pipe(lambda df: df[df.Geneid.str.startswith("FBgn")])  # Remove ERCCs
        .rename(columns={"Geneid": "FBgn"})
        .set_index("FBgn")
        .pipe(lambda df: df.loc[:, df.columns.str.endswith("TCP")])  # Remove ovary data
    )

    return np.log10(
        tpm(counts, gene_lengths.reindex(counts.index))
        .mean(axis=1)
        .pipe(lambda x: x[x > 0])
        .rename("log_avg_tpm_bulk")
    )


def get_avg_sc() -> pd.DataFrame:
    """Calculate the Log10 Average TPM for Summed scRNA-Seq"""
    gene_lengths = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    # I want to sum by replicate. I could do this all in pandas, but that
    # will use a lot of RAM. So I am going to do a split and combine approach.

    # Split cells by replicate
    cell_ids = {"rep1": [], "rep2": [], "rep3": []}
    for cell_id in pd.read_feather(snakemake.input.sc).set_index("FBgn").columns:
        if cell_id.startswith("rep4"):  # We did not end up using rep 4
            continue
        prefix = cell_id.split("_")[0]
        cell_ids[prefix].append(cell_id)

    # Sum cells by replicate
    counts = pd.concat(
        [
            pd.read_feather(snakemake.input.sc, columns=["FBgn"] + cells)
            .set_index("FBgn")
            .sum(axis=1)
            .rename(rep)
            for rep, cells in cell_ids.items()
        ],
        axis=1,
        sort=False,
    )

    return np.log10(
        tpm(counts, gene_lengths.reindex(counts.index))
        .mean(axis=1)
        .pipe(lambda x: x[x > 0])
        .rename("log_avg_tpm_sum_sc")
    )


def add_reg_plot(
    x: str, y: str, data: pd.DataFrame, label: str, ax: plt.Axes, **kwargs
) -> plt.Axes:
    rho, _ = spearmanr(data[x], data[y])
    stat = r" ($\rho = $" + f"{rho:.2f})"

    sns.regplot(
        x,
        y,
        data=data,
        ax=ax,
        scatter=False,
        line_kws=dict(**kwargs),
        label=label + stat,
    )
    return ax


def tweak_plot(ax: plt.Axes, xlim, ylim) -> plt.Axes:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    legend = ax.legend(frameon=False)
    plt.setp(legend.get_texts(), fontweight="bold")
    return ax


if __name__ == "__main__":
    main()

