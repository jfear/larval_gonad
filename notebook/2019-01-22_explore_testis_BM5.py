# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 0.8.6
#   kernelspec:
#     display_name: Python [conda env:larval_gonad]
#     language: python
#     name: conda-env-larval_gonad-py
# ---

# %% [markdown]
# # Explore Testis BM5

# %% [markdown]
# First look at BM5 sample. Just trying to see what cluster is what.

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
from larval_gonad.config import read_config
from larval_gonad.notebook import Nb
from larval_gonad.normalization import tpm, zscore
from larval_gonad.scRNAseq import TSNEPlot

# %%
nb = Nb.setup_notebook(seurat_dir='../output/translocations-wf/translocation_BM5_force')

# %%
metadata = nb.seurat.get_metadata()
clusters = metadata['res.0.6']
tsne = nb.seurat.get_tsne()
df = tsne.join(clusters)

# %%
sns.lmplot('tSNE_1', 'tSNE_2', data=df, hue='res.0.6', fit_reg=False, scatter_kws=dict(s=5), palette=sns.color_palette(n_colors=11))
plt.title('Testis BM5 tSNE');

# %%
lit_genes = read_config('../science_submission/config.yaml', 'lit_genes')
lit_fbgns = [nb.symbol2fbgn[g] for g in lit_genes]

# %%
raw_cnts = nb.seurat.get_raw()
raw_by_cluster = raw_cnts.T.join(clusters).groupby('res.0.6').sum().T
raw_by_cluster.index.name = 'Fbgn'
raw_by_cluster.columns.name = 'cluster'
raw_by_cluster.columns = raw_by_cluster.columns.astype(int)

# %%
ax = sns.heatmap(zscore(raw_by_cluster.reindex(lit_fbgns)), vmin=-3, vmax=3, cmap='viridis')
yticks = [
    nb.fbgn2symbol[y.get_text()]
    for y in ax.get_yticklabels()
]

ax.set_yticklabels(yticks)
plt.title('Testis BM5 Literature Genes');

# %%


