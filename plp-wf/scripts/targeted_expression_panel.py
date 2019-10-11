"""Expression panel of 4 target genes.

Brian G. asked me to make tSNE version of the expression panels for his
paper. He wants to focus on:

* spd-2
* asl
* sas4
* plp

I also want to output a colored version of the tSNE plot.

"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

from larval_gonad.config import read_config

plt.style.use("seaborn-white")
plt.rcParams["savefig.facecolor"] = "white"

COLOR_FILE = "../../config/colors.yaml"
CLUSTER_FILE = "../../output/seurat3-cluster-wf/combined_n3_clusters.feather"
NORM_FILE = "../../output/seurat3-cluster-wf/combined_n3_normalized.feather"
ZSCORE_FILE = "../../output/seurat3-cluster-wf/zscore_by_cell.feather"
TSNE_FILE = "../../output/plp-wf/tsne.tsv"


def main():
    tsne()

    target_genes = [
        ("FBgn0027500", "spd-2"),
        ("FBgn0261004", "asl"),
        ("FBgn0011020", "Sas-4"),
        ("FBgn0086690", "Plp"),
    ]

    fig, axes = plt.subplots(
        2, 2, figsize=(8, 8), sharex=True, gridspec_kw=dict(hspace=0.1, wspace=0.1)
    )
    for gene, ax in zip(target_genes, axes.flat):
        panel(*gene, ax)

    for ax in axes[:, 1].flat:
        ax.set_ylabel("")

    plt.savefig("../../output/plp-wf/short_panel.svg", bbox="tight", pad_inches=0, dpi=300)


def tsne():
    tsne = pd.read_csv(TSNE_FILE, sep="\t").set_index("cell_id")
    clusters = pd.read_feather(CLUSTER_FILE, columns=["cell_id", "cluster"]).set_index("cell_id")
    colors = {
        k: v for k, v in zip(clusters.cluster.cat.categories, read_config(COLOR_FILE)["clusters"])
    }
    df = tsne.join(clusters, how="inner").assign(color=lambda x: x.cluster.astype(str).map(colors))

    fig = plt.figure(figsize=(4, 4))
    plt.scatter(df.tSNE_1, df.tSNE_2, s=3, c=df.color, lw=0)
    ax = plt.gca()
    sns.despine(ax=ax, left=True, bottom=True)
    ax.set(xlabel="tSNE 1", ylabel="tSNE 2", aspect="equal", rasterized=True)
    plt.savefig("../../output/plp-wf/tSNE.svg", bbox="tight", pad_inches=0, dpi=300)


def panel(fbgn, symbol, ax):
    expression_patterns(fbgn, symbol, ax)
    axins = inset_axes(ax, 1.5, 1.5)
    tsne_zscore(fbgn, symbol, axins)


def expression_patterns(fbgn, symbol, ax):
    norm = (
        pd.read_feather(NORM_FILE)
        .query(f"FBgn == '{fbgn}'")
        .set_index(["FBgn", "gene_symbol"])
        .T.squeeze()
        .rename("norm")
    )
    clusters = pd.read_feather(CLUSTER_FILE, columns=["cell_id", "cluster"]).set_index("cell_id")
    df = clusters.join(norm)
    sns.pointplot("cluster", "norm", data=df, ax=ax)
    ax.set_title(f"{symbol}", fontstyle="italic")
    ax.set(xlabel="", ylabel="Normalized Expression (by cell)")
    sns.despine(ax=ax)
    return ax


def tsne_zscore(fbgn, symbol, ax):
    tsne = pd.read_csv(TSNE_FILE, sep="\t").set_index("cell_id")
    zscore = pd.read_feather(ZSCORE_FILE).query(f"FBgn == '{fbgn}'").set_index("FBgn").T.squeeze()
    df = tsne.join(zscore.rename("zscore"))

    # Plot
    cmap = plt.get_cmap("viridis", 512)
    norm = matplotlib.colors.Normalize(-3, 3)

    sns.scatterplot(
        x="tSNE_1",
        y="tSNE_2",
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
    import os

    try:
        os.chdir(os.path.join(os.getcwd(), "plp-wf/scripts"))
        print(os.getcwd())
    except:
        pass

    main()

