import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

FNAME = snakemake.input["ratios"]
PVALS = snakemake.input["pvals"]

CLUSTER_COLORS = snakemake.params.cluster_colors

ONAME = snakemake.output[0]


# Debug Settings
# FNAME = "output/x-to-a-wf/autosome_ratios_by_cell.feather"
# PVALS = "output/x-to-a-wf/permuted_autosome_ratio_pvalues.feather"
# import yaml
# config = yaml.safe_load(open("config/common.yaml"))
# CLUSTER_COLORS = yaml.full_load(open("config/colors.yaml"))["clusters"]


def main():
    df = pd.read_feather(FNAME).set_index("cell_id")

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(4, 6), sharex=True, sharey=False, gridspec_kw=dict(hspace=0.1)
    )

    _defaults = {
        "data": df,
        "palette": CLUSTER_COLORS,
        "notch": True,
        "linewidth": 0.5,
        "showfliers": False,
    }

    sns.boxplot("cluster", "x_to_a_ratio", ax=ax1, **_defaults)
    ax1.set_axisbelow(True)
    ax1.axhline(1, ls="--", lw=0.5, color="#b0b0b0", zorder=0)
    ax1.set(xlabel="", ylabel="X:A")
    ax1.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    sns.boxplot("cluster", "fourth_to_a_ratio", ax=ax2, **_defaults)
    ax2.set_axisbelow(True)
    ax2.axhline(1, ls="--", lw=0.5, color="#b0b0b0", zorder=0)
    ax2.set(xlabel="", ylabel="4:A")
    ax2.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    # Add p-value indicators
    whiskers = whisker_locations(df)
    df_pvals = get_pvals()
    for i, (cluster, dd) in enumerate(df_pvals.iterrows()):
        if dd.pval_x <= 0.05:
            pval_str = pval_to_string(dd.pval_x)
            loc = whiskers.loc[cluster, "x_to_a_ratio"]
            ax1.text(i, loc + 0.01, pval_str, ha="center", va="bottom", fontsize=6)

        if dd.pval_4 <= 0.05:
            pval_str = pval_to_string(dd.pval_4)
            loc = whiskers.loc[cluster, "fourth_to_a_ratio"]
            ax2.text(i, loc + 0.01, pval_str, ha="center", va="bottom", fontsize=6)

    fig.savefig(ONAME, bbox_inches="tight")


def pval_to_string(pval):
    if pval <= 0.001:
        return "***"
    elif pval <= 0.01:
        return "**"
    elif pval <= 0.05:
        return "*"
    else:
        return ""


def whisker_locations(df):
    """Calculate the upper whisker location for plotting astrix."""
    return (
        df.groupby("cluster")
        .quantile([0.25, 0.75])
        .rename_axis(["cluster", "quantile"])
        .reset_index()
        .melt(id_vars=["cluster", "quantile"], var_name="ratio", value_name="value")
        .pivot_table(index=["cluster", "ratio"], columns="quantile", values="value")
        .assign(iqr=lambda df: df[0.75] - df[0.25])
        .assign(upper=lambda df: df[0.75] + (1.5 * df["iqr"]))
        .upper.unstack()
    )


def get_pvals():
    """P-values from permutation test."""
    return pd.read_feather(PVALS).set_index("cluster")


if __name__ == "__main__":
    main()
