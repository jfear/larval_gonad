"""Merge to ptrap gene table

We want to add some additional columns to the ptrap gene table which is
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
ptrap_fbgns = read_config("../config/ptrap_genes.yaml")["ptraps"]

#%%
fbgn2symbol = pd.read_feather("../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()

# %%
biomarkers = (
    pd.read_feather(
        "../output/seurat3-cluster-wf/combined_n3_biomarkers.feather",
        columns=["FBgn", "cluster", "p_val_adj"],
    )
    .query("p_val_adj <= 0.01")
    .assign(cluster=lambda x: x.cluster.map(config["cluster_annot"]))
    .assign(p_val_adj=lambda x: x.p_val_adj.map(lambda y: f"({y})"))
    .assign(value=lambda x: x.cluster.str.cat(x.p_val_adj))
    .groupby("FBgn")
    .apply(lambda x: "|".join(x.value))
    .reindex(ptrap_fbgns)
    .rename("Biomarkers")
    .to_frame()
)

# %%
zscores = (
    pd.read_feather("../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather")
    .set_index(["FBgn", "rep", "cluster"])
    .unstack(level=[-1, -2])
    .reindex(ptrap_fbgns)
)

#%%
pscores = pd.read_csv("../data/external/miriam/ptrap_scores_2019-11-10.tsv", sep="\t", na_values="-").set_index("FBgn")


#%%
df = (
    biomarkers.join(zscores, how="left")
    .merge(pscores.reset_index(), left_on="FBgn", right_on="FBgn", how="right")
    .set_index(["FBgn", "gene_symbol", "ptrap_id"])
)
df.to_csv("../output/notebook/2019-11-20_ptrap_scores.tsv", sep="\t")
