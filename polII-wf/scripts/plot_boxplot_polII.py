"""Playing with plots for galletta"""
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_ind

from larval_gonad.plotting.stats import pval_to_string


def main():
    plt.style.use(["2c", "science_base"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(4.75, 1.16), sharey=True, gridspec_kw=dict(wspace=0.01))
    ax1.set_ylim(None, 2.3)

    # S2
    s2 = get_data(snakemake.input[0], "S2")
    sns.boxplot("chrom", "norm", data=s2, ax=ax1, palette=snakemake.params.colors, notch=True)
    ax1.set(xlabel="", title="S2 Phospho CTD / Total Pol II", ylabel="Normalized Pixel Intensity")
    sns.despine(ax=ax1)
    pval = get_pval(s2)
    add_pval_and_line(pval, ax1)

    # S5
    s5 = get_data(snakemake.input[0], "S5")
    sns.boxplot("chrom", "norm", data=s5, ax=ax2, palette=snakemake.params.colors, notch=True)
    ax2.set(xlabel="", title="S5 Phospho CTD / Total Pol II", ylabel="")
    sns.despine(ax=ax2)
    pval = get_pval(s5)
    add_pval_and_line(pval, ax2)

    fig.savefig(snakemake.output[0])


def get_data(file_name, antibody):
    return (
        pd.read_excel(file_name, sheet_name=antibody, skiprows=1)
        .melt(value_name="norm", var_name="chrom")
        .dropna()
        .assign(antibody=antibody)
    )


def get_pval(df):
    x = df.query(f'chrom == "X"').norm
    a = df.query(f'chrom == "A"').norm
    _, pval = ttest_ind(x, a, equal_var=False)
    return pval_to_string(pval)


def add_pval_and_line(pval, ax):
    ax.text(0.5, 0.99, pval, transform=ax.transAxes, va="top", ha="center")
    l = plt.Line2D([0.3, 0.7], [0.94, 0.94], transform=ax.transAxes, color="k", lw=0.8, ls="-")
    ax.add_line(l)
    return ax


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="polII-wf",
            input="../data/external/galletta/phos_over_tot_data_late_SC_only_061819.xlsx",
            params=dict(colors=["red", "grey"]),
        )

    plt.style.use("../config/figure_styles.mplstyle")

    main()
