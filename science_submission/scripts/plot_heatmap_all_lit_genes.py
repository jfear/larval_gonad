"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use("Agg")

from itertools import chain
from more_itertools import flatten
from pickle import loads
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import feather_to_cluster_rep_matrix

FNAME = snakemake.input.zscores
GENE_METADATA = snakemake.input.gene_metadata

LIT_GENES = snakemake.params.lit_genes
CMAP = snakemake.params.cmap

ONAME = snakemake.output[0]

# Debug settings
# FNAME = "output/science_submission/zscore_by_cluster_rep.feather"
# GENE_METADATA = "references/gene_annotation_dmel_r6-24.feather"
# import yaml
# LIT_GENES = yaml.safe_load(open("config/literature_genes.yaml"))
# CMAP = "viridis"


def main():
    fbgn2symbol = (
        pd.read_feather(GENE_METADATA, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .to_dict()["gene_symbol"]
    )
    zscores = (
        feather_to_cluster_rep_matrix(FNAME)
        .reindex(flatten(LIT_GENES.values()))
        .rename(fbgn2symbol)
    )

    plt.style.use("scripts/paper_1c.mplstyle")
    fig = plt.figure(figsize=(4, 8))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        zscores,
        xticklabels=True,
        yticklabels=True,
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
    ax.xaxis.set_tick_params(pad=0, length=2)
    ax.set_xticklabels(
        list(chain.from_iterable([("", x, "") for x in zscores.columns.levels[0]])),
        ha="center",
        va="bottom",
    )

    # Add lines separating cell types
    for i in range(1, 12):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    # Clean up Y axis
    ax.set_ylabel("")
    ax.yaxis.set_tick_params(pad=0.1, length=2)
    plt.setp(ax.get_yticklabels(), fontstyle="italic", fontfamily="Helvetica", rotation=0, va="center")

    # Add lines separating lit genes
    cols = zscores.shape[1]
    xloc = cols + cols * 0.01
    loc = 0
    for k, v in LIT_GENES.items():
        prev = loc
        loc += len(v)
        mid = prev + (loc - prev) / 2
        ax.axhline(loc, color="w", ls="--", lw=0.5)
        ax.text(xloc, mid, k, ha="left", va="center", fontweight="bold")

    # Clean up color bar
    cax.xaxis.set_tick_params(pad=0, length=2)

    fig.savefig(ONAME, bbox_inches="tight")


if __name__ == "__main__":
    main()
