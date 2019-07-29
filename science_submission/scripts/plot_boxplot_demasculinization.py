"""Demasculinization of the X by cell type.

Look at how and where testis biased expressed genes are expressed in our
different cell-type groups.

"""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


BULK_DEG = snakemake.input["bulk"]
RAW = snakemake.input["raw"]
NORM = snakemake.input["norm"]
CELL_ANNOTATION = snakemake.input["cell_annotation"]
GENE_ANNOTATION = snakemake.input["gene_annotation"]

CLUSTER_ANNOT = snakemake.params["cluster_annot"]
CLUSTER_ORDER = snakemake.params["cluster_order"]

PCT_FILE = snakemake.output["pct_cells"]
NORM_FILE = snakemake.output["norm_expression"]

BOXPLOT_DEFAULTS = dict(
    showfliers=False, notch=True, order=["X", "2L", "2R", "3L", "3R", "4"], linewdith=0.5
)

# Debug settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass
# BULK_DEG = ''
# RAW = '../../output/cellselection-wf/raw.feather'
# NORM = '../../output/seurat3-cluster-wf/combined_n3_normalized.feather'
# CELL_ANNOTATION = "../../output/seurat3-cluster-wf/combined_n3_metadata.feather"
# GENE_ANNOTATION = '../../references/gene_annotation_dmel_r6-24.feather'
# import yaml
# config = yaml.safe_load(open("../../config/common.yaml"))
# CLUSTER_ANNOT = config["cluster_annot"]
# CLUSTER_ORDER = config["cluster_order"]


def main():
    plt.style.use("scripts/figure_styles.mplstyle")
    cell2cluster = get_cell_to_cluster(CELL_ANNOTATION)
    fbgn2chrom = get_fbgn2chrom(GENE_ANNOTATION)
    testis_biased_fbgns = get_testis_biased(BULK_DEG)

    # Proportion of Cells with Expression
    prop_cells = get_proportion_cells(RAW, cell2cluster, fbgn2chrom, testis_biased_fbgns)
    prop_cells_subset = filter_by_cluster(prop_cells)

    g = sns.FacetGrid(data=prop_cells_subset, row="cluster")
    g.map(sns.boxplot, "chrom", "prop_cells", **BOXPLOT_DEFAULTS)
    g.savefig(PCT_FILE, bbox_inches="tight")

    # Normalized Expression
    norm_expression = get_norm_cells(NORM, cell2cluster, fbgn2chrom, testis_biased_fbgns)
    norm_expression_subset = filter_by_cluster(norm_expression)
    avg_norm_expression = (
        norm_expression_subset.groupby(["chrom", "FBgn", "cluster"])
        .norm_expression.mean()
        .reset_index()
    )

    g = sns.FacetGrid(data=avg_norm_expression, row="cluster")
    g.map(sns.boxplot, "chrom", "norm_expression", **BOXPLOT_DEFAULTS)
    g.savefig(NORM_FILE, bbox_inches="tight")


def get_testis_biased(fname, alpha=0.01):
    # TODO: Check filtering
    return (
        pd.read_feather(fname).query(f"adj_p-val <= {alpha} & avg_logFC > 0").FBgn.values.tolist()
    )


def get_cell_to_cluster(fname):
    return (
        pd.read_feather(fname)
        .set_index("cell_id")
        .assign(
            cluster=lambda df: pd.Categorical(
                df.cluster.map(CLUSTER_ANNOT), ordered=True, categories=CLUSTER_ORDER
            )
        )
        .cluster
    )


def get_proportion_cells(fname, cell2cluster, fbgn2chrom, target_fbgns):
    return (
        pd.read_feather(fname)
        .set_index("FBgn")
        .reindex(target_fbgns)
        .dropna()
        .T.pipe(lambda x: x > 0)
        .join(cell2cluster)
        .groupby("cluster")
        .mean()
        .T.rename_axis("FBgn")
        .stack()
        .rename("prop_cells")
        .reset_index(level=-1)
        .join(fbgn2chrom, how="inner")
    )


def get_norm_cells(fname, cell2cluster, fbgn2chrom, target_fbgns):
    return (
        pd.read_feather(fname)
        .set_index("FBgn")
        .drop("gene_symbol", axis=1)
        .rename_axis("cell_id", axis=1)
        .reindex(target_fbgns)
        .dropna()
        .stack()
        .rename("norm_expression")
        .reset_index(level=0)
        .join(cell2cluster)
        .set_index("FBgn")
        .join(fbgn2chrom, how="inner")
    )


def get_fbgn2chrom(fname):
    return pd.read_feather(fname).set_index("FBgn").FB_chrom.rename("chrom")


def filter_by_cluster(df):
    _df = df.query("cluster == ['SP', 'EPS', 'PS1', 'PS2', 'PS3', 'ECY', 'CY1', 'CY2']").copy()
    clusters = _df.cluster.astype(str).replace(
        dict(PS1="PSs", PS2="PSs", PS3="PSs", ECY="CYs", CY1="CYs", CY2="CYs")
    )
    _df["cluster"] = pd.Categorical(clusters, ordered=True, categories=["SP", "EPS", "PSs", "CYs"])
    return _df


if __name__ == "__main__":
    main()
