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
# # Missingness by Chromosome Arm

# %% [markdown]
# Missingness is still a major concern. Cameron's purposed permutation appears to be working, but we want to verify that missingness is not driving these differences. 

# %%
import os
import sys
import re
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %% [markdown]
# ## Create a list of FBgns by chromosome
#
# 1. Pull out the major autosomes and X.
# 2. Relabel autosomes as A to simplify things

# %%
# Create mapping of FBgn to X or A linked genes
autosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R']

fbgn2chrom = nbconfig.fbgn2chrom[nbconfig.fbgn2chrom.chrom.isin(autosomes + ['chrX'])].copy()
fbgn2chrom = fbgn2chrom.chrom.map(dict(chr2L='A', chr2R='A', chr3L='A', chr3R='A', chrX='X'))
fbgn2chrom.value_counts().map(lambda x: f'{x:,}')

# %% [markdown]
# ## Read in raw coverage counts for replicate 2

# %%
# read rep 2 raw data.
raw = pd.read_csv('../output/scrnaseq-wf/scrnaseq_rep2_force/raw.tsv', sep='\t', index_col=0)
raw.index.name = 'FBgn'
raw.columns.name = 'cell_id'

# %%
# Get list of cells in clusters
clusters = nbconfig.seurat.get_clusters('res.0.6').map(nbconfig.short_cluster_annot)
clusters = clusters[clusters != 'UNK'].copy()
clusters = clusters.astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order)
clusters.index.name = 'cell_id'

# %%
# only keep cells that actually are assigned to a cluster
raw = raw.reindex(columns=[x for x in clusters.index if x.startswith('rep2')]).copy()
raw_w_chrom = raw.join(fbgn2chrom, how='left')

# melt and munge
raw_melted = raw_w_chrom.reset_index().melt(id_vars=['FBgn', 'chrom'], value_name='UMI', var_name='cell_id').set_index('FBgn')

# %% [markdown]
# ## Calculate proportion of missingness of X and A

# %% [markdown]
# I consider all 0's as missing and simply calculate the proportion of 0's to total number of reads.

# %%
def get_prop_missing(df):
    x = df.query('chrom == "X"')
    x_genes = x.shape[0]
    x_missing = (x.UMI == 0).sum()
    x_prop_missing = x_missing / x_genes
    
    a = df.query('chrom == "A"')
    a_genes = a.shape[0]
    a_missing = (a.UMI == 0).sum()
    a_prop_missing = a_missing / a_genes
    
    return pd.Series([x_prop_missing, a_prop_missing], index=['x_prop_missing', 'a_prop_missing'])

# %%
# Calculate the proporiton missing by X or A
prop_missing = raw_melted.groupby('cell_id').apply(get_prop_missing)

# %%
ax = prop_missing.plot(kind='kde')
ax.set_xlabel('Proportion Missing')
plt.title('Distribution of Proportion Missing');

# %% [markdown]
# ## Missingness by cell type.

# %% [markdown]
# Next I take these cell cell missing counts and look at how they are distributed y cell type.

# %%
# look at proportion missing by cell type
prop_missing_w_clusters = prop_missing.join(clusters)

g = sns.FacetGrid(prop_missing_w_clusters, col='cluster', col_wrap=4)
g.map(sns.kdeplot, 'x_prop_missing', label='X')
g.map(sns.kdeplot, 'a_prop_missing', color='r', label='A')
g.axes[3].legend(loc=[1, 0.5])
g.set_xlabels('Proportion Missing')
plt.suptitle('Distribution of the Proportion of Missing by Cluster', va='bottom');

# %% [markdown]
# ## Logistic regression to see if your flag depeleted is driven by the missingness

# %% [markdown]
# Next I run a logistic regression comparing our flag for depleted X expression with our missingness measures.

# %%
# read in flag_depleted calls and merge on
flag_depleted = pd.read_csv('../output/notebook/2019-01-17_prototype_cell_type_permutation_test.csv', index_col=0).flag_depleted.astype(int)
df = prop_missing_w_clusters.join(flag_depleted)

# %%
df2.head()

# %%
# Simple model looking at just missingness.
result = smf.logit('flag_depleted ~ x_prop_missing + a_prop_missing', data = df).fit()
result.summary()

# %%
# Model adding in cluster identity
result = smf.logit('flag_depleted ~ x_prop_missing + a_prop_missing + cluster', data = df).fit()
result.summary()

# %%
# complex model adding interaction terms
result = smf.logit('flag_depleted ~ x_prop_missing*cluster + x_prop_missing + a_prop_missing*cluster + a_prop_missing', data = df).fit()
result.summary()

# %%


