"""Merge to Literature gene table

We want to add some additional columns to the literature gene table which is
a excel file. I am just going to make outputs in the correct order for easy
pasting.
"""
#%%
import os

from more_itertools import flatten
import numpy as np
import pandas as pd

from larval_gonad.config import read_config

#%%
try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass


# %%
config = read_config("../config/common.yaml")

lit_genes = read_config("../config/literature_genes.yaml")
lit_fbgns = list(flatten([v for k, v in lit_genes.items()]))

# %%
biomarkers = (
    pd.read_feather(
        "../output/seurat3-cluster-wf/combined_n3_biomarkers.feather",
        columns=["FBgn", "cluster", "p_val_adj"],
    )
    .assign(cluster=lambda x: x.cluster.map(config["cluster_annot"]))
    .assign(p_val_adj=lambda x: x.p_val_adj.map(lambda y: f"({y})"))
    .assign(value=lambda x: x.cluster.str.cat(x.p_val_adj))
    .groupby("FBgn")
    .apply(lambda x: "|".join(x.value))
    .reindex(lit_fbgns)
    .rename("Biomarkers")
    .to_frame()
)
biomarkers.to_csv("../output/notebook/2019-11-19_litgene_merge_biomarkers.tsv", sep="\t")

# %%
zscores = (
    pd.read_feather("../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather")
    .set_index(["FBgn", "rep", "cluster"])
    .unstack(level=[-1, -2])
    .reindex(lit_fbgns)
)
zscores.to_csv("../output/notebook/2019-11-19_litgene_merge_zscores.tsv", sep="\t")