"""Plot UMAP Panels"""
import matplotlib as mpl

mpl.use("Agg")

from more_itertools import flatten, grouper
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

UMAP = snakemake.input.umap
GENE_ANNOTATION = snakemake.input["gene_annot"]
NORM_FILE = snakemake.input["norm"]

LIT_GENES = snakemake.params["lit_genes"]
FILE_PATTERN = snakemake.params["file_pattern"]


# Debug settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass
# UMAP = '../../output/seurat3-cluster-wf/combined_n3_umap.feather'
# import yaml
# LIT_GENES = yaml.safe_load(open('../../config/literature_genes.yaml'))
# GENE_ANNOTATION = '../../references/gene_annotation_dmel_r6-26.feather'
# NORM_FILE = '../../output/seurat3-cluster-wf/combined_n3_normalized.feather'


def main():
    fbgns_of_interest = get_fbgns_of_interest()
    fbgn2symbol = (
        pd.read_feather(GENE_ANNOTATION)
        .set_index("FBgn")
        .reindex(fbgns_of_interest)
        .dropna()
        .gene_symbol
        .str.lower()
        .sort_values()
    )
    fbgns_of_interest_sorted = fbgn2symbol.index.tolist()

    umap = pd.read_feather(UMAP).set_index("cell_id")
    norm = get_normalized(NORM_FILE, fbgns_of_interest_sorted).join(umap)

    plt.style.use("scripts/figure_styles.mplstyle")
    vmin, vmax = 0, 5
    cmap = sns.light_palette("blue", as_cmap=True)

    for i, grp in enumerate(grouper(18, fbgns_of_interest_sorted)):
        df = norm.query(f"FBgn == {grp}")
        g = sns.FacetGrid(df, col="gene_symbol", col_wrap=3)
        g.map(
            facet_scatter,
            "UMAP_1",
            "UMAP_2",
            "norm",
            s=5,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            rasterized=True,
            edgecolor=None,
        )
        g.set_titles("{col_name}")

        # Make space for the colorbar
        g.fig.subplots_adjust(right=0.92)

        # Define a new Axes where the colorbar will go
        cax = g.fig.add_axes([0.95, 0.80, 0.02, 0.1])

        # Get a mappable object with the same colormap as the data
        points = plt.scatter([], [], c=[], vmin=vmin, vmax=vmax, cmap=cmap)

        # Draw the colorbar
        g.fig.colorbar(points, cax=cax)

        g.savefig(FILE_PATTERN.replace("PANELID", str(i + 1)), bbox_inches="tight")


def get_fbgns_of_interest():
    return list(set(flatten([value for key, value in LIT_GENES.items()])))


def get_normalized(fname, target_fbgns):
    return (
        pd.read_feather(fname)
        .set_index("FBgn")
        .reindex(target_fbgns)
        .dropna()
        .set_index("gene_symbol", append=True)
        .stack()
        .rename("norm")
        .reset_index(level=[0, 1])
        .rename_axis("cell_id")
    )


def facet_scatter(x, y, c, **kwargs):
    """Draw scatterplot with point colors from a faceted DataFrame columns."""
    kwargs.pop("color")
    idx = np.argsort(c)
    plt.scatter(x[idx], y[idx], c=c[idx], **kwargs)


if __name__ == "__main__":
    main()
