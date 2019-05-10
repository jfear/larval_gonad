# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.3
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

# %%
# Project level imports
from larval_gonad.notebook import Nb
from larval_gonad.stats import run_chisq

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
clusters = (
    pd.read_parquet('../output/scrnaseq-wf/clusters.parquet')
    .assign(rep=lambda df: df.index.str.extract('(rep\d)', expand=False))
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order + ['UNK']))
)

# %%
number_of_cells = clusters.groupby(['cluster', 'rep']).size().rename('number_of_cells').to_frame()
number_of_cells

# %%
ct = number_of_cells.unstack().T
ct.index = ct.index.droplevel(0)
res = run_chisq(ct.T).loc[(slice(None), ['adj std residual', 'flag_sig']), :]

summary = {'flag_sig_bias': {}}
for clus, dd in res.groupby('cluster'):
    dd.index = dd.index.droplevel(0)
    for rep, ddd in dd.T.iterrows():
        summary['flag_sig_bias'][(clus, rep)] = 0
        if ddd.flag_sig:
            if ddd['adj std residual'] < 0:
                summary['flag_sig_bias'][(clus, rep)] = -1
            elif ddd['adj std residual'] > 0:
                summary['flag_sig_bias'][(clus, rep)] = 1

flag_sig_bias = (
    pd.DataFrame(summary)
    .rename_axis(['cluster', 'rep'])
    .reset_index()
    .assign(cluster=lambda df: pd.Categorical(df.cluster, ordered=True, categories=nbconfig.short_cluster_order + ['UNK']))
    .set_index(['cluster', 'rep'])
)
flag_sig_bias

# %%
summary_stats = (
    pd.read_parquet('../output/scrnaseq-wf/raw_by_cluster_w_rep.parquet')
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order + ['UNK']))
    .groupby(['cluster', 'rep']).describe()
)
summary_stats.columns = summary_stats.columns.droplevel(0)
summary_stats = summary_stats.drop('count', axis=1)
summary_stats

# %%
df = number_of_cells.join(flag_sig_bias).join(summary_stats)

# %%
df


# %%
def _calc_prop(x):
    return x / x.sum()
    
plot_dat = (
    df.groupby('cluster').number_of_cells.apply(_calc).reset_index()
    .rename(columns={'number_of_cells': 'proportion_of_cells'})
    .set_index(['cluster', 'rep'])
    .unstack()
)
plot_dat.columns = plot_dat.columns.droplevel(0)

# %%
ax = plot_dat.plot(kind='bar', stacked=True, width=.9)
ax.set(xlabel='Cluster', ylabel='Proportion of Cells', title='Replicate Representation')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
ax.legend(loc='upper left', bbox_to_anchor=[1, 1]);

# %%
df

# %%
