"""Expression panel of 6 ptrap examples in Figure 1.

* p53
* Piezo
* Nlg3
* Bol
* Rdo
* Efa6

I also want to output a colored version of the tSNE plot.

"""
import os

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

plt.style.use("seaborn-white")
plt.rcParams["savefig.facecolor"] = "white"


def main():
    target_genes = [
        ("FBgn0039044", "p53"),
        ("FBgn0264953", "Piezo"),
        ("FBgn0083963", "Nlg3"),
        ("FBgn0011206", "bol"),
        ("FBgn0243486", "rdo"),
        ("FBgn0051158", "Efa6"),
    ]

    fig, axes = plt.subplots(
        3, 2, figsize=(8, 8), sharex=True, gridspec_kw=dict(hspace=0.1, wspace=0.1)
    )
    for gene, ax in zip(target_genes, axes.flat):
        panel(*gene, ax)

    for ax in axes[:, 1].flat:
        ax.set_ylabel("")

    plt.savefig(snakemake.output[0], bbox="tight", pad_inches=0, dpi=300)


def panel(fbgn, symbol, ax):
    expression_patterns(fbgn, symbol, ax)
    axins = inset_axes(ax, 1.5, 1.5)
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
    ax.set_title(f"{symbol}", fontstyle="italic")
    ax.set(xlabel="", ylabel="Normalized Expression (by cell)")
    sns.despine(ax=ax)
    return ax


def umap_zscore(fbgn, symbol, ax):
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

    return ax


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", None):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
                zscores="../output/seurat3-cluster-wf/zscore_by_cell.feather",
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                norm="../output/seurat3-cluster-wf/combined_n3_normalized.feather",
            ),
        )

    main()
