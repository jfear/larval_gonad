# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.0
#   kernelspec:
#     display_name: Python [conda env:larval_gonad]
#     language: python
#     name: conda-env-larval_gonad-py
# ---

# %%
import os
import sys
import re
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
# spearman correlation M1 and L1
(
    pd.read_parquet('../output/scrnaseq-wf/tpm.parquet')
    .assign(cluster=lambda df: df.cluster.map(nbconfig.short_cluster_annot))
    .query('cluster == ["M1º", "L1º"]')
    .pivot_table(index='FBgn', columns='cluster', values='TPM')
    .corr(method='spearman')
)

# %%
# Genes that are biomarkers in M1 and L1
df = (
    pd.read_csv('../output/scrnaseq-wf/scrnaseq_combine_force/biomarkers_res.0.6.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .assign(cluster=lambda df: df.cluster.map(nbconfig.short_cluster_annot))
    #.query('cluster == ["M1º", "L1º"]')
    .pivot_table(index='FBgn', columns='cluster', values='avg_logFC')
    .pipe(lambda df: ~df.isnull())
)

pd.crosstab(df["M1º"], df["L1º"])

# %%

# %%
male_sterile = (
    pd.read_csv('/home/fearjm/Downloads/male_sterile_genes.txt', header=None)
    .iloc[:, 0]
    .tolist()
)

# %%
# Genes that are expressed
df = (
    pd.read_parquet('../output/scrnaseq-wf/raw_by_cluster.parquet')
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .dropna()
    .pivot_table(index='FBgn', columns='cluster', values='UMI')
    .fillna(0)
    .pipe(lambda df: df > 0)
)

# %%
df.reindex(male_sterile).dropna().sum(axis=1).value_counts()

# %%
display(pd.crosstab(df["E1º"], df["M1º"]))
display(pd.crosstab(df["E1º"], df["L1º"]))
display(pd.crosstab(df["M1º"], df["L1º"]))

# %%
dcc = [
    'FBgn0002774',
    'FBgn0002775',
    'FBgn0005616',
    'FBgn0005617',
    'FBgn0014340',
    'FBgn0019660',
    'FBgn0019661',
    'FBgn0283442'
]



# %%
(
    pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet')
    .assign(cluster=lambda df: df.cluster.map(nbconfig.short_cluster_annot))
    .reset_index()
    .query('FBgn == "FBgn0283442"')
)



# %%
(
    pd.read_parquet('../output/scrnaseq-wf/raw_by_cluster.parquet')
    .assign(cluster=lambda df: df.cluster.map(nbconfig.short_cluster_annot))
    .query('cluster != "UNK"')
    .pivot_table(index='FBgn', columns='cluster', values='UMI')
    .reindex(columns=nbconfig.short_cluster_order)
    .reindex(dcc)
    .rename(nbconfig.fbgn2symbol)
)




# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
df.sum()[nbconfig.short_cluster_order]

# %%
(
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/mid_vs_late.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .query('p_val_adj <= 0.05')
    .assign(mid_bias=lambda df: df.avg_logFC > 0)
    .mid_bias.value_counts()
)


# %%

# %%

# %%

# %%

# %%

# %%
