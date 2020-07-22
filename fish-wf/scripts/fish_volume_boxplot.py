import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ttest_rel

import larval_gonad.plotting  # pylint: disable=unused-import

plt.style.use("minimal")

BOXPLOT_KWS = dict(order=["X", "2L"], palette=snakemake.params.colors, notch=True)


def main():
    volume = pd.read_csv(snakemake.input[0])

    stacked = pd.concat(
        [
            stack_data("X", "2L", volume, "um3"),
            stack_data("X_scaled_probe", "2L_scaled_probe", volume, "um3_scaled_probe"),
            stack_data(
                "X_scaled_probe_count",
                "2L_scaled_probe_count",
                volume,
                "um3_scaled_probe_count",
            ),
        ],
        ignore_index=True,
    )

    pvals = pd.concat(
        [
            run_stats("X", "2L", volume, "um3"),
            run_stats("X_scaled_probe", "2L_scaled_probe", volume, "um3_scaled_probe"),
            run_stats(
                "X_scaled_probe_count",
                "2L_scaled_probe_count",
                volume,
                "um3_scaled_probe_count",
            ),
        ],
        ignore_index=True,
    )

    g = sns.FacetGrid(stacked, col="metric", sharey=False)
    g.map(sns.boxplot, "chrom", "volume", **BOXPLOT_KWS)
    g.set_titles("{col_name}")
    g.set_xlabels("")
    add_pvals(g, pvals)

    plt.savefig(snakemake.output[0])


def stack_data(c1: str, c2: str, df: pd.DataFrame, name: str) -> pd.DataFrame:
    _df = df[[c1, c2]]
    _df.columns = ["X" if col.startswith("X") else "2L" for col in _df.columns]
    return _df.melt(var_name="chrom", value_name="volume").assign(metric=name)


def run_stats(c1: str, c2: str, df: pd.DataFrame, name: str) -> pd.DataFrame:
    base_mean = df[[c1, c2]].sum().sum() / df.values.ravel().size
    lfc = np.log2(df[c1].mean() / df[c2].mean())
    pval = ttest_rel(df[c1], df[c2])[1]
    return pd.DataFrame(
        {
            "metric": [name],
            "baseMean": [base_mean],
            "log2FoldChange": [lfc],
            "p_value": [pval],
        }
    )


def add_pvals(g: sns.FacetGrid, pvals: pd.DataFrame):
    for ax in g.axes.flat:
        metric = ax.get_title()
        pval = pvals.query("metric == @metric")["p_value"].values[0]
        if pval <= 0.05:
            # Extend axis and add p-value.
            _, _max = ax.get_ylim()
            ax.set_ylim(None, _max + (0.05 * _max))

            ax.text(
                0.5,
                0.99,
                f"p = {pval:0.4f}",
                transform=ax.transAxes,
                va="center",
                ha="center",
            )
            l = plt.Line2D(
                [0.3, 0.7],
                [0.93, 0.93],
                transform=ax.transAxes,
                color="k",
                lw=0.9,
                ls="-",
            )
            ax.add_line(l)


if __name__ == "__main__":
    main()
