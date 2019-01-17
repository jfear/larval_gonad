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
# # Cell Type Level Permutation Test

# %% [markdown]
# After talking with Cameron he suggested a different permutation algorithm as follows.
#
# 1. For each cell calculate Mann-Whitney U on normalized counts for X vs A. Use the p-value to determine if a cell is depleted of X.
# 2. For each cell type cluster, calculate the proportion of cells with depleted X. 
# 3. Permute cluster labels and build null distributions of proportion of depleted cells for each cluster size.
# 4. Use the corresponding null distribution to calculate emperical p-values for each cluster.
#
# The one thing that needs worked out is how to "normalize" cell counts. In principle, we have the following things we can normalize by:
# * (X) number of reads per cell
#     * I am doing on calculations within a cell so I don't think we need to account for this
# * (X) gene length
#     * 10X is 3' biased, so I don't think gene length really needs to be accounted for
# * (X) chromosome length
#     * These counts are at the gene level, so chromosome length does not come into play
# * (Maybe) number genes per chromosome
# * (Maybe) number of expressed genes per chromosome
# * (Maybe) proportion of genes expressed per chromosome (number of expressed genes per chromosome / number genes per chromosome)
#     * This combines the other two (Maybe) into a single value.
#
# I will try without normalization first and see how different clusters behave. I will then explore using the (Maybe) values to affect some kind of normalization. 
#
# I also need to decide how to incorporate replicate information.

# %%
import os
import sys
import re
from pathlib import Path
from collections import namedtuple

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

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
# melt and munge
raw_melted = raw.reset_index().melt(id_vars='FBgn', value_name='UMI').set_index('FBgn')
raw_melted_expressed = raw_melted[raw_melted.UMI > 0]
raw_melted_expressed_w_chrom = raw_melted_expressed.join(fbgn2chrom)

# %%
raw_melted_expressed_w_chrom.head()

# %% [markdown]
# ## Estimate X chromosome depletion
#
# For each cell use the Mann-Whitney U test to determine if the median X-linked gene expression is less than the median Autosome linked gene expression. Create a flag `flag_depleted` indicating that X-linked genes were depleted (i.e., Mann-Whitney was significant).
#
# *Note: I required that there are at least 100 genes for X and A, but all cells met this criteria.*
#

# %%
# Get the proportion of cells with depelete X chromosome expression
results = []
for cell_id, dd in raw_melted_expressed_w_chrom.groupby('cell_id'):
    if dd.chrom.unique().shape[0] == 1:
        continue
        
    x_genes = dd.query('chrom == "X"').UMI.values
    a_genes = dd.query('chrom == "A"').UMI.values

    if x_genes.shape[0] < 100 and a_genes.shape[0] < 100:
        results.append((cell_id, np.nan))
        continue

    stat, p_value = mannwhitneyu(x_genes, a_genes, alternative='less')
    if p_value < 0.05:
        results.append((cell_id, True))
    else:
        results.append((cell_id, False))

flag_depleted = pd.DataFrame(results, columns=['cell_id', 'flag_depleted']).set_index('cell_id').flag_depleted
flag_depleted.value_counts()

# %% [markdown]
# ## Calculate the proportion of cells with X chromosome depletion for each cell type cluster
#
# Using `flag_depleted`, I merge on cell type information and calculate the proportion of cells that showed depletion.

# %%
# Merge on cell type info
clusters = nbconfig.seurat.get_clusters('res.0.6').map(nbconfig.short_cluster_annot)
clusters = clusters[clusters != 'UNK'].copy()
clusters = clusters.astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order)

flag_depleted_w_cluster = pd.concat([flag_depleted, clusters], axis=1, join='inner')
obs_prop = flag_depleted_w_cluster.groupby('cluster').flag_depleted.sum() / flag_depleted_w_cluster.groupby('cluster').flag_depleted.size()
display(HTML('<h4>Proportion of Cells Depleted by Cluster</h4>'))
obs_prop

# %% [markdown]
# ### Permute celltype labels and create a null distribution of proportion of cell depleted
#
# Now I scramble the celltype labels and recalcualte the proportion of cells depleted for each cell type. This creates a null distribution for each celltype which controls for cluster size.

# %%
results = []
for i in range(10_000):
    _df = flag_depleted_w_cluster.copy()
    _df['cluster'] = _df['cluster'].sample(frac=1, replace=False).values
    props = _df.groupby('cluster').flag_depleted.sum() / _df.groupby('cluster').flag_depleted.size()
    results.append(props)
perm = pd.concat(results, axis=1).T.reset_index(drop=True)
perm.plot(kind='kde', title='Null Distributions')

# %% [markdown]
# ### Compare observed proportions to null distribution
#
# For each celltype compare the observed proportion of depleted cells with celltype null distribution. Calculate the empirical p-value as the proportion of permuted samples with a more extreme value.

# %%
results = []
for cluster in nbconfig.short_cluster_order:
    p_value = (perm[cluster] >= obs_prop[cluster]).sum() / perm.shape[0]
    results.append((cluster, p_value))
pd.DataFrame(results, columns=['cluster', 'p_value']).set_index('cluster')

# %% [markdown]
# ## Permutation of Germline Only

# %%
flag_depleted_w_cluster_germline_only = flag_depleted_w_cluster[flag_depleted_w_cluster.cluster.isin(['SP', 'ES', 'MS', 'LS'])]

# %% [markdown]
# ### Permute celltype labels and create a null distribution of proportion of cell depleted
#
# Now I scramble the celltype labels and recalcualte the proportion of cells depleted for each cell type. This creates a null distribution for each celltype which controls for cluster size.

# %%
results = []
for i in range(10_000):
    _df = flag_depleted_w_cluster_germline_only.copy()
    _df['cluster'] = _df['cluster'].sample(frac=1, replace=False).values
    props = _df.groupby('cluster').flag_depleted.sum() / _df.groupby('cluster').flag_depleted.size()
    results.append(props.dropna())
perm = pd.concat(results, axis=1).T.reset_index(drop=True)
perm.plot(kind='kde', title='Null Distributions')

# %% [markdown]
# ### Compare observed proportions to null distribution
#
# For each celltype compare the observed proportion of depleted cells with celltype null distribution. Calculate the empirical p-value as the proportion of permuted samples with a more extreme value.

# %%
results = []
for cluster in nbconfig.short_cluster_order[:4]:
    p_value = (perm[cluster] >= obs_prop[cluster]).sum() / perm.shape[0]
    results.append((cluster, p_value))
pd.DataFrame(results, columns=['cluster', 'p_value']).set_index('cluster')

# %%
flag_depleted.to_frame().to_csv('../output/notebook/2019-01-17_prototype_cell_type_permutation_test.csv')

# %%

