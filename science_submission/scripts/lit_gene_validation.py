"""Cell type validation code.

This module houses the cell type validation code.
"""
import os
from functools import partial
import numpy as np
import pandas as pd

from larval_gonad.config import read_config
from larval_gonad.io import feather_to_cluster_rep_matrix
from larval_gonad.validation import GeneValidator

try:
    os.chdir(os.path.join(os.getcwd(), "src/larval_gonad"))
    print(os.getcwd())
except:
    pass


def main():
    cluster_annot = read_config("../../config/common.yaml")["cluster_annot"]
    fbgn2symbol = (
        pd.read_feather(
            "../../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]
        )
        .set_index("FBgn")
        .squeeze()
    )
    lit_genes = get_lit()
    biomarkers = get_biomarkers()
    zscores = get_zscores()

    res = []
    for fbgn, lit_gene in lit_genes.iterrows():
        if lit_gene.method == "protein":
            protein = True
        else:
            protein = False
        gene_score = GeneValidator(
            fbgn,
            lit_gene.expressed_in,
            lit_gene.missing,
            biomarkers.get(fbgn, set()),
            zscores.get(fbgn, set()),
            flag_protein=protein,
        )
        res.append(
            [fbgn, gene_score.lit_gene, gene_score.biomarker, gene_score.zscore, gene_score.score]
        )
    df = (
        pd.DataFrame(res, columns=["FBgn", "literature", "Biomarker", "Zscore", "Score"])
        .fillna(0)
        .set_index("FBgn")
        .join(fbgn2symbol)
        .set_index("gene_symbol", append=True)
    )
    df.Score.value_counts()


def get_lit():
    return (
        pd.read_csv("../../data/external/miriam/lit_gene_table.csv")
        .set_index("FBgn")
        .assign(expressed_in=lambda x: x.expressed_in.str.split("|").apply(lambda y: set(y)))
        .assign(missing=lambda x: x.not_analyzed.fillna("").str.split("|").apply(lambda y: set(y)))
        .loc[:, ["expressed_in", "missing", "method"]]
    )


def get_biomarkers():
    return (
        pd.read_feather("../../output/seurat3-cluster-wf/combined_n3_biomarkers.feather")
        .query("p_val_adj <= 0.01")
        .assign(cluster=lambda x: x.cluster.map(cluster_annot))
        .groupby("FBgn")
        .apply(lambda x: set(x.cluster))
        .rename("sc")
    )


def cell_type_above_upper_quantile(x):
    x["flag_upper"] = x.zscore > np.quantile(x.zscore, 0.75)
    flag_cluster = x.groupby(["cluster"]).flag_upper.sum() == 3
    return set(flag_cluster[flag_cluster].index.tolist())


def get_zscores():
    return (
        pd.read_feather("../../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather")
        .groupby("FBgn")
        .apply(cell_type_above_upper_quantile)
        .rename("zscore")
    )


if __name__ == "__main__":
    main()
