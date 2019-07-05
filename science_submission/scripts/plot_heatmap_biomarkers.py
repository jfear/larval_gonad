"""Plot heatmap of unique biomarkers"""
import matplotlib

matplotlib.use("Agg")

from itertools import chain
from more_itertools import flatten
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import KMeans

from larval_gonad.io import feather_to_cluster_rep_matrix

FNAME = snakemake.input.zscores
GENE_METADATA = snakemake.input.gene_metadata
BIOMARKERS = snakemake.input.biomarkers

CLUSTER_ANNOT = snakemake.params.cluster_annot
CLUSTER_ORDER = snakemake.params.cluster_order
LIT_GENES = snakemake.params.lit_genes
CMAP = snakemake.params.cmap

ONAME = snakemake.output[0]

# Debug settings
# import os
# os.chdir('science_submission/scripts')
# FNAME = "../../output/science_submission/zscore_by_cluster_rep.feather"
# GENE_METADATA = "../../references/gene_annotation_dmel_r6-24.feather"
# BIOMARKERS = "../../output/seurat3-cluster-wf/combined_n3_biomarkers.feather"
# import yaml
# config = yaml.safe_load(open('../../config/common.yaml'))
# CLUSTER_ANNOT = config['cluster_annot']
# CLUSTER_ORDER = config['cluster_order']
# LIT_GENES = yaml.safe_load(open('../../config/literature_genes.yaml'))
# CMAP = "viridis"


def main():
    global fbgn2symbol
    fbgn2symbol = (
        pd.read_feather(GENE_METADATA, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .to_dict()["gene_symbol"]
    )

    biomarkers = (
        pd.read_feather(BIOMARKERS, columns=["FBgn", "cluster"])
        .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
        .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
    )

    global lit_genes
    lit_genes = pd.concat((
        pd.DataFrame({'cell_type': k}, index=pd.Index(v, name='FBgn'))
        for k, v in LIT_GENES.items()
    ))

    zscores = feather_to_cluster_rep_matrix(FNAME).reindex(biomarkers.FBgn.unique())

    km = run_kmeans(zscores)

    # Reindex and sort FBgns by KMean class
    zscores.index = pd.MultiIndex.from_arrays([zscores.index.values, km], names=['FBgn', 'K'])
    zscores_ordered = zscores.sort_index(level='K')

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

    # Clean up Y axis
    ax.set_ylabel("")

    # Add lines separating cell types (columns)
    for i in range(1, len(CLUSTER_ORDER)):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    # Add lines separating KMeans groups (rows)
    loc = 0
    for name, group in zscores_ordered.groupby(level='K'):
        loc += group.shape[0]
        ax.axhline(loc, color='w', ls='--', lw=0.5)

    # Add literature genes that belong to each KMeans group
    lit_annot = litGeneAnnotator(zscores_ordered)
    lit_annot.add_annotation(ax)

    # Clean up color bar
    cax.xaxis.set_tick_params(pad=0, length=2)

    fig.savefig(ONAME, bbox_inches="tight")


def run_kmeans(df):
    km = KMeans(n_clusters=10, random_state=42)
    km.fit(df)
    return km.labels_


def get_lit_genes(df):
    return list(
        map(
            lambda x: fbgn2symbol[x],
            (
                lit_genes
                .reindex(df.index.get_level_values('FBgn'))
                .dropna()
                .index.values
            )
        )
    )


class litGeneAnnotator:
    def __init__(self, df):
        self.xlocL = 0
        self.tlocL = self.xlocL - 5
        self.xlocR = df.shape[1]
        self.tlocR =  self.xlocR + 5
        self.grouper = df.groupby(level='K')

    def add_annotation(self, ax):
        loc = 0
        for name, group in self.grouper:
            size = group.shape[0]
            loc += size
            lit_genes_in_group = '\n'.join(get_lit_genes(group))
            if len(lit_genes_in_group) <= 1:
                continue

            # Figure out coordinates
            yloc = loc - (size / 2)
            tloc, xloc, ha = self.get_xloc(name)

            ax.annotate(
                lit_genes_in_group, 
                (xloc, yloc), 
                (tloc, yloc), 
                arrowprops={
                    'arrowstyle': 'simple',
                    'color': 'black'
                },
                ha=ha, 
                va='center', 
                fontsize=10
            )

    def get_xloc(self, name):
        if name % 2 == 0:
            tloc = self.tlocR
            self.tlocR += 5
            return tloc, self.xlocR, 'left'
        else:
            tloc = self.tlocL
            self.tlocL -= 5
            return tloc, self.xlocL, 'right'



if __name__ == "__main__":
    main()
