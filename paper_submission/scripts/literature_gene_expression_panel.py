"""Expression panel of the entire literature gene table."""
import os
from pathlib import Path

from more_itertools import grouper
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

plt.style.use("seaborn-white")
plt.rcParams["savefig.facecolor"] = "white"

LIT_GENES = [
    "FBgn0250816",
    "FBgn0000146",
    "FBgn0283442",
    "FBgn0000158",
    "FBgn0039044",
    "FBgn0026533",
    "FBgn0087008",
    "FBgn0010258",
    "FBgn0031091",
    "FBgn0031715",
    "FBgn0032473",
    "FBgn0034282",
    "FBgn0041103",
    "FBgn0262743",
    "FBgn0286813",
    "FBgn0004372",
    "FBgn0002842",
    "FBgn0039124",
    "FBgn0026573",
    "FBgn0002673",
    "FBgn0011206",
    "FBgn0003475",
    "FBgn0011596",
    "FBgn0012037",
    "FBgn0019828",
    "FBgn0041102",
    "FBgn0051361",
    "FBgn0000405",
    "FBgn0002862",
    "FBgn0011725",
    "FBgn0030313",
    "FBgn0034435",
    "FBgn0034739",
    "FBgn0039071",
    "FBgn0000404",
    "FBgn0031623",
    "FBgn0000416",
    "FBgn0030520",
    "FBgn0250843",
    "FBgn0261885",
    "FBgn0266521",
    "FBgn0000320",
    "FBgn0083963",
    "FBgn0051158",
    "FBgn0264494",
    "FBgn0004872",
    "FBgn0011591",
    "FBgn0038197",
    "FBgn0039754",
    "FBgn0000964",
    "FBgn0243486",
    "FBgn0001090",
    "FBgn0004108",
    "FBgn0010453",
    "FBgn0050418",
    "FBgn0264953",
    "FBgn0000636",
    "FBgn0010473",
    "FBgn0000576",
    "FBgn0004885",
    "FBgn0024288",
    "FBgn0001257",
    "FBgn0003984",
    "FBgn0004647",
    "FBgn0014163",
    "FBgn0015399",
    "FBgn0020493",
    "FBgn0027598",
    "FBgn0041182",
    "FBgn0243512",
    "FBgn0264975",
    "FBgn0265487",
]


def main():
    target_genes = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .squeeze()
        .reindex(LIT_GENES)
        .items()
    )

    cnt = 0
    for genes in grouper(16, target_genes):
        fig, axes = plt.subplots(4, 4, sharex=True, figsize=(12, 12), gridspec_kw=dict(hspace=0.05, wspace=0.1))
        for gene, ax in zip(genes, axes.flat):
            if gene is None:
                ax.set_visible(False)
            else:
                panel(*gene, ax)
                for ax in axes[:, 1:].flat:
                    ax.set_ylabel("")

                for ax in axes[-1, :].flat:
                    plt.setp(ax.get_xticklabels(), rotation=90)

        plt.savefig(snakemake.params.pattern.format(cnt=cnt), transparent=True)
        cnt += 1

    Path(snakemake.output[0]).touch()


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
    ax.set_title(f"{symbol}", fontstyle="italic", y=.9)
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
            workdir="paper_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
                zscores="../output/seurat3-cluster-wf/zscore_by_cell.feather",
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                norm="../output/seurat3-cluster-wf/combined_n3_normalized.feather",
            ),
        )

    main()
