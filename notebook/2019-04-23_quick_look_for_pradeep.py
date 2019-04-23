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

# %% [markdown]
# # Quick look for Pradeep

# %% [markdown]
# Pradeep send the following:
#
# > I wanted to have look at scRNAseq data for the genes listed below:
# >
# > 1. DSX (CG11094,FBgn0000504)
# > 2. tra (CG16724,FBgn0003741)
# > 3. broad (CG11491,FBgn0283451)
# > 
# > Could you help me figuring out their expression on the basis of cell type.

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
pradeep = [
    'FBgn0000504',
    'FBgn0003741',
    'FBgn0283451',
]

# %%
dat = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_w_rep.parquet')
    .assign(cluster=lambda df: df.cluster.map(nbconfig.short_cluster_annot))
    .query('cluster != "UNK"')
)

# %%
df = (
    pd.pivot_table(dat, index='FBgn', columns=['cluster', 'rep'], values='TPM')
    .reindex(pradeep)
    .assign(gene_symbol=lambda df: df.index.map(nbconfig.fbgn2symbol))
    .set_index('gene_symbol', append=True)
    .reindex(nbconfig.short_cluster_order, level='cluster', axis=1)
)

# %%
df

# %%
df.to_excel('../output/notebook/2019-04-23_quick_look_for_pradeep.xlsx')

# %%
