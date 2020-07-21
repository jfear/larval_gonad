"""Expression panel of the entire literature gene table."""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import larval_gonad.plotting  # pylint: disable=unused-import


plt.style.use("minimal")
sns.set_style("dark")

FBGNS = [
    "FBgn0039044",  # p53
    "FBgn0243486",  # rdo
    "FBgn0026573",  # ADD1
    "FBgn0083963",  # Nlg3
    "FBgn0011206",  # bol
    "FBgn0264953",  # Piezo
]

SYMBOLS = ["p53", "rdo", "ADD1", "Nlg3", "bol", "Piezo"]


def main():
    df = (
        pd.read_feather(snakemake.input.zscores)
        .query("FBgn in @FBGNS")
        .pipe(add_gene_symbol)
        .pipe(add_umap)
        .sort_values("z-score")
    )

    scatter_defaults = dict(
        s=3, linewidth=0, rasterized=True, vmin=-3, vmax=3, cmap="viridis"
    )

    g = sns.FacetGrid(
        data=df,
        col="gene_symbol",
        col_wrap=2,
        col_order=SYMBOLS,
        sharex=True,
        sharey=True,
    )
    g.map(facet_scatter, "UMAP_1", "UMAP_2", "z-score", **scatter_defaults)
    add_colorbar(g, scatter_defaults)
    tweak_axes(g)

    g.savefig(snakemake.output[0])


def add_gene_symbol(df: pd.DataFrame) -> pd.DataFrame:
    fbgn2symbol = pd.read_feather(
        snakemake.input.gene_annot, columns=["FBgn", "gene_symbol"]
    )
    return df.merge(fbgn2symbol)


def add_umap(df: pd.DataFrame) -> pd.DataFrame:
    umap = pd.read_feather(snakemake.input.umap)
    return df.melt(
        id_vars=["FBgn", "gene_symbol"], var_name="cell_id", value_name="z-score"
    ).merge(umap)


def facet_scatter(x, y, c, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    plt.scatter(x, y, c=c, **kwargs)


def add_colorbar(g: sns.FacetGrid, scatter_defaults: dict):
    g.fig.subplots_adjust(bottom=0.2)
    # cax = g.fig.add_axes([0.94, 0.25, 0.02, 0.6])
    cax = g.fig.add_axes([0.20, 0.10, 0.60, 0.02])
    points = plt.scatter([], [], c=[], **scatter_defaults)
    g.fig.colorbar(points, cax=cax, label="Z-Score (TPM)", orientation="horizontal")
    return g


def tweak_axes(g: sns.FacetGrid):
    g.set_titles("{col_name}", fontstyle="italic")
    g.set_xlabels("UMAP 1")
    g.set_ylabels("UMAP 2")
    g.despine(left=True, bottom=True)
    for ax in g.axes.ravel():
        ax.set(xticks=[-10, 0, 10], yticks=[-10, 0, 10], aspect=1)
    return g


if __name__ == "__main__":
    main()
