"""Expression panel of the entire literature gene table."""
import os

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

matplotlib.use("Agg")
plt.style.use("seaborn-white")
plt.rcParams["savefig.facecolor"] = "white"


def main():
    fbgn = snakemake.wildcards.fbgn
    gene_symbol = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .squeeze()
        .loc[fbgn]
    )

    _, ax = plt.subplots( 1, 1, figsize=(4, 4))
    panel(fbgn, gene_symbol, ax)
    plt.setp(ax.get_xticklabels(), rotation=90)

    plt.savefig(snakemake.output[0], transparent=True)


def panel(fbgn, symbol, ax):
    expression_patterns(fbgn, symbol, ax)
    axins = inset_axes(ax, 2, 2, loc=7)
    umap_zscore(fbgn, symbol, axins)


def expression_patterns(fbgn, symbol, ax):
    norm = (
        pd.read_feather(snakemake.input.norm)
        .query(f"FBgn == '{fbgn}'")
        .set_index(["FBgn", "gene_symbol"])
        .T.squeeze()
        .rename("norm")
    )
    clusters = pd.read_feather(snakemake.input.clusters, columns=["cell_id", "cluster"]).set_index(
        "cell_id"
    )
    df = clusters.join(norm)
    sns.pointplot("cluster", "norm", data=df, ax=ax)
    ax.set_title(f"{symbol}", fontstyle="italic", y=0.9)
    ax.set(xlabel="", ylabel="Normalized Expression (by cell)")
    sns.despine(ax=ax)
    return ax


def umap_zscore(fbgn: str, symbol: str, ax: plt.Axes):
    umap = pd.read_feather(snakemake.input.umap).set_index("cell_id")
    zscore = (
        pd.read_feather(snakemake.input.zscores)
        .query(f"FBgn == '{fbgn}'")
        .set_index("FBgn")
        .T.squeeze()
    )
    df = umap.join(zscore.rename("zscore"))

    # Plot
    cmap = plt.get_cmap("viridis", 512)
    norm = matplotlib.colors.Normalize(-3, 3)

    sns.scatterplot(
        x="UMAP_1",
        y="UMAP_2",
        data=df.sort_values("zscore"),
        hue="zscore",
        hue_norm=norm,
        palette=cmap,
        s=3,
        linewidth=0,
        rasterized=True,
        legend=False,
        ax=ax,
    )
    sns.despine(ax=ax, left=True, bottom=True)
    ax.set(xlabel="", xticks=[], ylabel="", yticks=[], aspect="equal")
    ax.set_facecolor('none')

    return ax


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", None):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
                zscores="../output/seurat3-cluster-wf/zscore_by_cell.feather",
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                norm="../output/seurat3-cluster-wf/combined_n3_normalized.feather",
            ),
            wildcards=dict(
                fbgn="FBgn0000504",
                symbol="dsx"
            )
        )

    main()
