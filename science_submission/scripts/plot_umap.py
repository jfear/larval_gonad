"""Plot UMAP"""
import matplotlib as mpl

mpl.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

UMAP = snakemake.input.umap
METADATA = snakemake.input.metadata
CLUSTER_COLORS = snakemake.params.cluster_colors
CLUSTER_ANNOT = snakemake.params.cluster_annot
CLUSTER_ORDER = snakemake.params.cluster_order

OUTPUT_FILE = snakemake.output[0]

# Debug settings
# UMAP = 'output/seurat3-cluster-wf/combined_n3_umap.feather'
# METADATA = 'output/seurat3-cluster-wf/combined_n3_metadata.feather'
# import yaml
# CLUSTER_COLORS = yaml.full_load(open('config/colors.yaml'))['clusters']
# config = yaml.safe_load(open('config/common.yaml'))
# CLUSTER_ANNOT = config['cluster_annot']
# CLUSTER_ORDER = config['cluster_order']


def main():

    color_mapper = dict(zip(CLUSTER_ORDER, CLUSTER_COLORS))

    cell_annot = (
        pd.read_feather(METADATA, columns=["cell_id", "cluster"])
        .set_index("cell_id")
        .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
        .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
        .assign(colors=lambda df: df.cluster.map(color_mapper))
    )
    cell_annot.head()

    umap = (
        pd.read_feather(UMAP).set_index("cell_id").join(cell_annot, how="right")
        .assign(UMAP_1=lambda df: df.UMAP_1 * -1)
    )

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, ax = plt.subplots()
    ax.scatter(
        umap.UMAP_1,
        umap.UMAP_2,
        c=umap.colors,
        s=3,
        linewidth=0.02,
        edgecolor="k",
        rasterized=True,
    )

    # Add text for cluster names
    for clus, row in (
        umap.groupby("cluster").agg({"UMAP_1": "mean", "UMAP_2": "mean"}).iterrows()
    ):
        plt.text(
            row.UMAP_1,
            row.UMAP_2,
            clus,
            bbox=dict(facecolor=(1, 1, 1, 0.8), edgecolor="none", pad=0.2),
            ha="center",
            va="center",
            fontsize=7,
            fontweight="bold",
        )

    # clean up plot
    plt.setp(
        ax,
        xticks=[],
        yticks=[],
        xlabel="",
        ylabel="",
        aspect="equal",
        xmargin=0,
        ymargin=0,
    )
    sns.despine(ax=ax, bottom=True, left=True)

    fig.savefig(OUTPUT_FILE, bbox_inches="tight")


if __name__ == "__main__":
    main()
