"""Plot total expression (TPM) for different gene sets."""
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting.stats import add_pvals
import larval_gonad.plotting  # pylint: disable=unused-import


plt.style.use("minimal")


def main():
    df = pd.concat(
        [
            pd.read_feather(
                snakemake.input.tau, columns=["fourth_to_a_ratio", "cluster"]
            ).assign(gene_set="Tau"),
            pd.read_feather(
                snakemake.input.tsps, columns=["fourth_to_a_ratio", "cluster"]
            ).assign(gene_set="TSPS"),
            pd.read_feather(
                snakemake.input.widely_expressed,
                columns=["fourth_to_a_ratio", "cluster"],
            ).assign(gene_set="Widely Expressed"),
        ],
        sort=False,
        ignore_index=True,
    )

    g = sns.FacetGrid(df, col="gene_set", col_order=["Tau", "TSPS", "Widely Expressed"])
    g.map(
        sns.boxplot,
        "cluster",
        "fourth_to_a_ratio",
        order=snakemake.params.order,
        palette=snakemake.params.colors,
        notch=True,
        linewidth=0.5,
        showfliers=False,
    )
    add_stats(g)
    tweak_axes(g)

    g.savefig(snakemake.output[0])


def tweak_axes(g: sns.FacetGrid):
    g.set_titles("{col_name}")
    g.set_xlabels("")
    g.set_ylabels("44/AA")

    ax: plt.Axes
    for ax in g.axes.flat:
        new_labels = [
            snakemake.params.names[label.get_text()] for label in ax.get_xticklabels()
        ]
        ax.set_xticklabels(new_labels, rotation=45)
        ax.set_axisbelow(True)
        ax.axhline(1, ls="--", color="grey")

    for ax in g.axes.flat[1:]:
        sns.despine(ax=ax, left=True)


def add_stats(g: sns.FacetGrid):
    # Tau
    df = shelve_load(snakemake.input.tau_pvals)["pvalues"].query(
        "ratio_type == 'fourth_to_a_ratio'"
    )
    add_pvals(df.x, df.y, df.pvalue, g.axes[0][0])

    # TSPS
    df = shelve_load(snakemake.input.tsps_pvals)["pvalues"].query(
        "ratio_type == 'fourth_to_a_ratio'"
    )
    add_pvals(df.x, df.y, df.pvalue, g.axes[0][1])

    # Widely Expressed
    df = shelve_load(snakemake.input.widely_expressed_pvals)["pvalues"].query(
        "ratio_type == 'fourth_to_a_ratio'"
    )
    add_pvals(df.x, df.y, df.pvalue, g.axes[0][2])


if __name__ == "__main__":
    main()
