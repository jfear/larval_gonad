"""Playing with plots for galletta"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.weightstats import ttest_ind

INPUT_FILE = snakemake.input[0]
OUTPUT_FILE = snakemake.output[0]

# Debug settings
# import os
# os.chdir('science_submission/scripts')
# INPUT_FILE = "../../data/external/galletta/phos_over_tot_data_late_SC_only_061819.xlsx"


def main():
    df = get_data()

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(4, 6), sharex=True, sharey=True, gridspec_kw=dict(hspace=0.1)
    )

    _defaults = {
        "notch": True,
        "linewidth": 0.5,
        "showfliers": False,
        "order": ["X", "A"],
        "color": "C0",
    }

    sns.boxplot("chrom", "norm", data=df.query('antibody == "S2"'), ax=ax1, **_defaults)
    sns.boxplot("chrom", "norm", data=df.query('antibody == "S5"'), ax=ax2, **_defaults)

    ax1.set(xlabel="", ylabel="S2 Phospho CTD / Total Pol II\nNormalized Pixel Intensity")
    ax2.set(xlabel="", ylabel="S5 Phospho CTD / Total Pol II\nNormalized Pixel Intensity")
    ax1.set_ylim(-.1, 2.2)

    # Add p-vals
    add_pval(get_pval(df, "S2"), ax1)
    add_pval(get_pval(df, "S5"), ax2)

    fig.savefig(OUTPUT_FILE, bbox_inches="tight")


def get_data():
    s2 = (
        pd.read_excel(INPUT_FILE, sheet_name="S2", skiprows=1)
        .melt(value_name="norm", var_name="chrom")
        .dropna()
        .assign(antibody="S2")
    )

    s5 = (
        pd.read_excel(INPUT_FILE, sheet_name="S5", skiprows=1)
        .melt(value_name="norm", var_name="chrom")
        .dropna()
        .assign(antibody="S5")
    )

    return pd.concat([s2, s5])


def add_pval(pval, ax):
    l1 = lines.Line2D([0.25, 0.75], [0.95, 0.95], color="k", transform=ax.transAxes)
    ax.add_line(l1)
    ax.text(0.5, 0.95, pval, ha="center", va="center", transform=ax.transAxes)
    return ax


def pval_to_string(pval):
    if pval <= 0.001:
        return "***"
    elif pval <= 0.01:
        return "**"
    elif pval <= 0.05:
        return "*"
    else:
        return ""


def get_pval(df, antibody):
    x = df.query(f'chrom == "X" & antibody == "{antibody}"').norm
    a = df.query(f'chrom == "A" & antibody == "{antibody}"').norm
    _, pval, _ = ttest_ind(x, a, usevar="unequal")
    return pval_to_string(pval)


if __name__ == "__main__":
    main()
