"""Cell type validation code.

This module houses the cell type validation code.
"""
import os
from functools import partial
import numpy as np
import pandas as pd

from larval_gonad.config import read_config
from larval_gonad.io import feather_to_cluster_rep_matrix

try:
    os.chdir(os.path.join(os.getcwd(), "src/larval_gonad"))
    print(os.getcwd())
except:
    pass


class gene_validator(object):
    def __init__(self, FBgn, lit_gene, biomarkers, zscores):
        self.fbgn = FBgn
        self.lit_gene = lit_gene
        self.biomarker = biomarkers.get(self.fbgn, set())
        self.zscore = zscores.get(self.fbgn, set())
        self.score = None
        self.validate()

    def validate(self):
        if self.lvl1():
            return

        if self.lvl2():
            return

        if self.lvl3():
            return

    def lvl1(self):
        """Lit and Biomarkers exact match"""
        if self.lit_gene == self.biomarker:
            self.score = 4
            return True

    def lvl2(self):
        """Lit and upper zscore quantail exact match"""
        if self.lit_gene == self.zscore:
            self.score = 3
            return True

    def lvl3(self):
        pass

    def __str__(self):
        return f"{self.fbgn}\t{self.lit_gene}\t{self.biomarker}\t{self.zscore}\t{self.score}"


def main():
    cluster_annot = read_config("../../config/common.yaml")["cluster_annot"]
    lit_genes = get_lit()
    biomarkers = get_biomarkers()
    zscores = get_zscores()

    for fbgn, lit_gene in lit_genes.items():
        gene_score = gene_validator(fbgn, lit_gene, biomarkers, zscores)
        print(gene_score)


def get_lit():
    return (
        pd.read_csv("../../data/external/miriam/lit_gene_table.csv")
        .set_index("FBgn")
        .assign(expressed_in=lambda x: x.expressed_in.str.split("|").apply(lambda y: set(y)))
        .expressed_in.squeeze()
        .rename("lit")
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
