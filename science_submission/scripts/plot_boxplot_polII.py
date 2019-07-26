"""Playing with plots for galletta"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
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

    # Add p-vals
    whiskers = df.groupby(["chrom", "antibody"]).norm.apply(whisker_max)
    s2 = get_pval(df, "S2")
    add_whisker(0, 1, whiskers[("X", "S2")], whiskers[("A", "S2")], s2, ax1)

    s5 = get_pval(df, "S5")
    add_whisker(0, 1, whiskers[("X", "S5")], whiskers[("A", "S5")], s5, ax2)

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


def whisker_max(x):
    lower, upper = np.quantile(x, [0.25, 0.75])
    iqr = upper - lower
    return upper + (1.5 * iqr)


def add_whisker(x1, x2, w1, w2, pval, ax):
    y_low = max(w1, w2) + 0.01 * max(w1, w2)
    y_high = y_low + 0.1
    x_mid = x1 + ((x2 - x1) / 2)
    ax.plot([x1, x1], [y_low, y_high], c="lightgray")
    ax.plot([x2, x2], [y_low, y_high], c="lightgray")
    ax.plot([x1, x2], [y_high, y_high], c="lightgray")
    ax.text(x_mid, y_high, pval, ha="center", va="bottom")
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
