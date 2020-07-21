""""""
from typing import Tuple, List
from collections import namedtuple
from itertools import product

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.plotting.stats import add_pvals

plt.style.use("minimal")


COLORS = ["black", "lightgray", "white"]

CHROM_ORDER = ["X", "Y", "2L", "2R", "3L", "3R", "4"]

TISSUE_MAPPER = {
    "WB": "Whole Body",
    "HD": "Head",
    "TX": "Thorax",
    "AC": "Abdomen",
    "DG": "Viscera",
    "RE": "Reproductive Tract",
    "GE": "Terminalia",
    "GO": "Gonad",
    "L3_GO": "L3 Gonad",
}

ORDER = [
    f"{species}:{tissue}"
    for tissue, species in product(TISSUE_MAPPER.keys(), ["w1118", "OreR"])
]

CHI_STR = r"$\chi^2$"

PostHoc = namedtuple("PostHoc", ["x", "y", "qval"])


def main():
    stats = pd.read_csv(snakemake.input.stats, sep="\t")
    proportions = (
        pd.read_feather(snakemake.input.data)
        .pipe(_count_number_genes_with_bias)
        .pipe(_fill_missing_Y_counts)
        .pipe(_calculate_proportion_of_genes)
        .reset_index()
        .assign(group=lambda df: df.species + ":" + df.tissue)
    )

    g = make_facet_grid(proportions)
    add_stats_to_plot(g, stats)
    make_plot_pretty(g)

    g.savefig(snakemake.output[0])


def make_facet_grid(proportions: pd.DataFrame) -> sns.FacetGrid:
    # Panel Plot species x tissue
    g = sns.FacetGrid(
        proportions,
        col="group",
        col_order=ORDER,
        col_wrap=6,
        sharex=False,
        sharey=False,
        aspect=0.5,
        height=3
    )
    g.map(_plot_stacked_bar, "chrom", "male", "NS", "female")
    g.set_titles("{col_name}")
    return g


def add_stats_to_plot(g: sns.FacetGrid, stats: pd.DataFrame):
    ax: plt.Axes
    for ax in g.axes.flat:
        species, tissue = ax.title.get_text().split(":")
        if species == "OreR" and tissue == "L3_GO":
            # There is no data here so remove from plot
            ax.set_visible(False)
            continue

        chisq, post_hoc = _wrangle_stats(
            stats.query("species == @species and tissue == @tissue")
        )
  
        ax.text(0.5, 1, f"({CHI_STR} {chisq:0.4f})", transform=ax.transAxes, ha="center", va="top")

        if chisq > 0.01:
            # Ignore stats if chi-square is not significant.
            continue

        add_pvals(
            post_hoc.x,
            post_hoc.y,
            post_hoc.qval,
            ax,
            sig_format=False,
            fontsize=12,
            fontweight="bold",
            color="w",
        )

    return g


def make_plot_pretty(g: sns.FacetGrid):
    g.despine(left=True)

    # Title
    for ax in g.axes.flat:
        species, tissue = ax.title.get_text().split(":")
        ax.set_title(f"{TISSUE_MAPPER[tissue]}\n{species}")

    # Axis labels
    g.set_xlabels("")
    g.set_ylabels("Proportion Genes")

    # Tick Labels
    xticks = list(g.axes[0].get_xticklabels())
    for ax in g.axes.flat[:-1]:  # Remember the last plot is empty
        ax.set_xticklabels(xticks, rotation=0)

    # Add Legend
    g.add_legend()
    return g


def _count_number_genes_with_bias(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["species", "tissue", "chrom"])
        .bias.value_counts()
        .rename("num_genes")
        .to_frame()
    )


def _fill_missing_Y_counts(df: pd.DataFrame) -> pd.DataFrame:
    _df = df.unstack(level="chrom").fillna(0).stack().unstack(level="bias")

    # If all 0 then add a 1 to the NS column for a prettier plot.
    all_zero = _df.sum(axis=1) == 0
    _df.loc[all_zero, ("num_genes", "NS")] = 1
    return _df.stack()


def _calculate_proportion_of_genes(df: pd.DataFrame) -> pd.DataFrame:
    cnts = df.unstack(level="bias").droplevel(level=0, axis=1).rename_axis(columns="")
    return cnts.div(cnts.sum(axis=1), axis=0).fillna(0)


def _plot_stacked_bar(chrom, male, ns, female, **kwargs):
    data = (
        pd.concat([chrom, male, ns, female], axis=1, sort=False)
        .set_index("chrom")
        .reindex(CHROM_ORDER)
        .rename(columns={"male": "Male-Bias", "female": "Female-Bias"})
    )
    ax = plt.gca()
    return data.plot.bar(
        stacked=True, ax=ax, color=COLORS, width=0.9, edgecolor="k", linewidth=0.3,
    )


def _wrangle_stats(stats: pd.DataFrame) -> Tuple[float, PostHoc]:
    table_pval = stats.table_pval.unique()[0]

    qvals = (
        stats.set_index("chrom")
        .qval_male.squeeze()
        .reindex(CHROM_ORDER)
        .fillna(1)
        .values
    )

    # Locations on the plot to add the *
    x = list(range(len(qvals)))
    y = [0.2,] * len(qvals)

    return table_pval, PostHoc(x, y, qvals)


if __name__ == "__main__":
    main()
