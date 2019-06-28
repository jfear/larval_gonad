"""Plot heatmap of unique biomarkers"""
import matplotlib

matplotlib.use("Agg")

from itertools import chain
from more_itertools import flatten
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, dendrogram

from larval_gonad.io import feather_to_cluster_rep_matrix

FNAME = snakemake.input.zscores
GENE_METADATA = snakemake.input.gene_metadata
BIOMARKERS = snakemake.input.biomarkers

CLUSTER_ANNOT = snakemake.params.cluster_annot
CLUSTER_ORDER = snakemake.params.cluster_order
CLUSTER_TYPES = snakemake.params.cluster_types
CMAP = snakemake.params.cmap

ONAME = snakemake.output[0]

# Debug settings
# FNAME = "output/science_submission/zscore_by_cluster_rep.feather"
# GENE_METADATA = "references/gene_annotation_dmel_r6-24.feather"
# BIOMARKERS = "output/seurat3-cluster-wf/combined_n3_biomarkers.feather"
# import yaml
# config = yaml.safe_load(open('config/common.yaml'))
# CLUSTER_ANNOT = config['cluster_annot']
# CLUSTER_ORDER = config['cluster_order']
# CLUSTER_TYPES = config['cluster_types']

# CMAP = "viridis"


def main():
    fbgn2symbol = (
        pd.read_feather(GENE_METADATA, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .to_dict()["gene_symbol"]
    )

    biomarkers = (
        pd.read_feather(BIOMARKERS, columns=["FBgn", "cluster"])
        .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
        .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
        .pipe(lambda df: df[df.FBgn.duplicated(keep=False)])
        .assign(cell_type=lambda df: df.cluster.map(dict(zip(CLUSTER_ORDER, CLUSTER_TYPES))))
    )

    # Re-label based on cell type (Germ Only, Soma Only, Germ and Soma)
    biomarkers_by_celltype = (
        biomarkers.groupby("FBgn")
        .apply(lambda df: df.cell_type.sort_values().str.cat(sep="|"))
        .apply(cell_type_mapper)
        .rename("cell_type")
        .to_frame()
        .reset_index()
    )

    zscores = feather_to_cluster_rep_matrix(FNAME).reindex(biomarkers.FBgn.unique())

    # order zscores by doing a hierarchical cluster for each cell type group.
    zscores_ordered = pd.concat(
        (
            hierarchal_cluster(zscores.reindex(v.FBgn.unique()))
            for k, v in biomarkers_by_celltype.groupby("cell_type")
        )
    )

    plt.style.use("scripts/figure_styles.mplstyle")
    fig = plt.figure(figsize=(4, 8))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        zscores_ordered,
        xticklabels=True,
        yticklabels=False,
        vmin=-3,
        vmax=3,
        rasterized=True,
        cmap=CMAP,
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(label="Z-Score (TPM)", ticks=[-3, 0, 3], orientation="horizontal"),
    )

    # Clean up X axis
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        list(chain.from_iterable([("", x, "") for x in CLUSTER_ORDER])), ha="center", va="bottom"
    )

    # Add lines separating cell types
    for i in range(1, len(CLUSTER_ORDER)):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    # Clean up Y axis
    ax.set_ylabel("")

    # Add lines separating biomarker groups
    loc = 0
    cols = zscores_ordered.shape[1]
    xloc = cols + cols * 0.01
    for clus, dd in biomarkers_by_celltype.groupby("cell_type"):
        prev = loc
        loc += dd.shape[0]
        mid = loc - ((loc - prev) / 2)
        ax.axhline(loc, color="w", ls="--", lw=0.5)
        txt = f"{clus} ({dd.shape[0]:,})"
        ax.text(xloc, mid, txt, ha="left", va="center", fontweight="bold")

    # Clean up color bar
    cax.xaxis.set_tick_params(pad=0, length=2)

    fig.savefig(ONAME, bbox_inches="tight")


def hierarchal_cluster(df):
    # cluster genes bases on expression
    link = linkage(df.values, "average")
    tree = dendrogram(link, no_plot=True)
    leaves = tree["leaves"]
    return df.iloc[leaves, :]


def cell_type_mapper(cell_type):
    if ("germ" in cell_type) and "soma" in cell_type:
        return "Germ and Soma"
    elif "germ" in cell_type:
        return "Germ Only"
    elif "soma" in cell_type:
        return "Soma Only"


if __name__ == "__main__":
    main()
