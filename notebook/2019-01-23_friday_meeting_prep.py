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
# # Friday Meeting Prep

# %% [markdown]
# I want to report on chromosome expression and X:A in this weeks Friday meeting. Here is where I am developing those plots.

# %%
import os
import sys
import re
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
def read_fbgn2chrom():
    mapper = {
        'chrX': 'X',
        'chrY': 'Y',
        'chr4': '4',
        'chr2L': 'A',
        'chr2R': 'A',
        'chr3L': 'A',
        'chr3R': 'A',
    }

    fbgn2chrom = (pd.read_csv('../output/fbgn2chrom.tsv', sep='\t', index_col=0)
                      .query('chrom != "chrM"')
                      .chrom.map(mapper)
                      .astype('category')
                      .cat.as_ordered()
                 )
    
    return fbgn2chrom.cat.reorder_categories(['X', 'A', 'Y', '4'])


def read_clusters():
    clusters = nbconfig.seurat.get_clusters('res.0.6').map(nbconfig.short_cluster_annot)
    clusters = clusters[clusters != 'UNK'].copy()
    return clusters.astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order)


def read_raw(rep2):
    raw = nbconfig.seurat.get_raw()
    if rep2:
        raw = raw.loc[:, raw.columns.str.startswith('rep2')].copy()
        
    return raw
        
    
def read_gene_length(): 
    gene_lengths = pd.read_csv('../output/gene_ts_lengths.tsv', sep='\t', index_col=0).gene_ts_length
    gene_lengths.name = 'gene_length'
    return gene_lengths
    
    
def read_tpm(rep2):
    from larval_gonad.normalization import tpm
    raw = read_raw(rep2)
    gene_lengths = read_gene_length()
    return tpm(raw, gene_lengths).dropna()
    
def get_rep(wide):    
    rep = wide.columns.str.extract('(?P<rep>rep\d)').rep
    rep.index = wide.columns
    return rep
    
def read_data(rep2=False, tpm=False):
    fbgn2chrom = read_fbgn2chrom()
    clusters = read_clusters()
    
    if tpm:
        data = read_tpm(rep2)
        value_name = 'TPM'
    else:
        data = read_raw(rep2)
        value_name = 'UMI'
    
    # Munge together
    rep = get_rep(data)
    melted = data.reset_index().melt(id_vars='FBgn', value_name=value_name)
    return melted.join(clusters, on='cell_id').join(fbgn2chrom, on='FBgn').join(rep, on='cell_id')

# %% [markdown]
# ## Data Prep

# %%
df = read_data()

# %%
fbgn2chrom = read_fbgn2chrom()
fbgn2chrom = fbgn2chrom.reindex(df.FBgn.unique())
num_genes_by_chrom = fbgn2chrom.value_counts()

# %%
total_reads_per_chrom_by_cell = df.groupby(['cell_id', 'chrom']).UMI.sum()
total_reads_per_cell = df.groupby(['cell_id']).UMI.sum()

# %%
norm_cnts = (
    total_reads_per_chrom_by_cell
        .div(num_genes_by_chrom / 1e3, level='chrom')
        .div(total_reads_per_cell / 1e3, level='cell_id')
        .to_frame()
)
norm_cnts.columns = ['norm_cnt']

norm_cnts = (
    norm_cnts
        .join(read_clusters(), on='cell_id')
        .reset_index()
)
norm_cnts = norm_cnts.join(norm_cnts.cell_id.str.extract('(?P<rep>rep\d)'))

norm_cnts.chrom = (
    norm_cnts.chrom
        .astype('category')
        .cat.as_ordered()
        .cat.reorder_categories(['X', 'A', 'Y', '4'])
)

# %% [markdown]
# ## Cell level chromosome coverage

# %%
g = sns.FacetGrid(norm_cnts, col='chrom', col_wrap=2, sharey=False)
g.map(
    sns.barplot, 
    'cluster', 
    'norm_cnt', 
    order=nbconfig.short_cluster_order, 
    palette=nbconfig.colors['clusters'],
    estimator=np.mean,
    errwidth=1,
    capsize=.2
)

# %% [markdown]
# # Rep level chromosome coverage

# %%
dat = norm_cnts.groupby(['cluster', 'rep', 'chrom']).norm_cnt.median().to_frame().reset_index()

# %%
g = sns.FacetGrid(dat, col='chrom', col_wrap=2, sharey=False)
g.map(
    sns.barplot, 
    'cluster', 
    'norm_cnt', 
    order=nbconfig.short_cluster_order, 
    palette=nbconfig.colors['clusters'],
    estimator=np.mean,
    errwidth=1,
    capsize=.2
)

# %%


# %%


# %% [markdown]
# ## Missingness is still problematic

# %%
df['missing'] = (df.UMI == 0).values

# %% [markdown]
# ### Missingness by cluster

# %%
missing_per_cell = df.groupby(['cell_id', 'cluster']).missing.sum().div(num_genes_by_chrom.sum(), level='chrom')
missing_per_cell.name = 'prop_missing'

# %%
dat = missing_per_cell.reset_index()
ax = sns.boxplot('cluster', 'prop_missing', data=dat, flierprops=dict(alpha=.5), palette=nbconfig.colors['clusters'])
#plt.setp(ax.artists, edgecolor='k', facecolor='w')
#plt.setp(ax.lines, color='k');

# %% [markdown]
# ### Missingness by cluster by chromosome

# %%
missing_per_cell_per_chrom = df.groupby(['cell_id', 'cluster', 'chrom']).missing.sum().div(num_genes_by_chrom, level='chrom')
missing_per_cell_per_chrom.name = 'prop_missing'

# %%
dat = missing_per_cell_per_chrom.reset_index()

g = sns.FacetGrid(dat, col='cluster', col_wrap=4)
g.map(
    sns.boxplot,
    'chrom',
    'prop_missing',
    order=['X', 'A', 'Y', '4'],
    flierprops=dict(alpha=.5)
)

for ax in g.axes:
    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')


# %% [markdown]
# ### Missingness is correlated between X and A

# %%
dat = missing_per_cell_per_chrom.to_frame().query('chrom == "X" | chrom == "A"').unstack()
dat.columns = ['prop_X_missing', 'prop_A_missing']
dat.reset_index(inplace=True)

# %%
def add_rho(color, marker, data):
    cluster = data.cluster.values[0]
    corr = spearmanr(data.prop_X_missing, data.prop_A_missing)[0]
    ax = plt.gca()
    ax.text(0.1, .9, f'r = {np.round(corr, 4)}', fontsize=12)
    
g = sns.lmplot(
    'prop_X_missing', 
    'prop_A_missing', 
    dat, 
    col='cluster', 
    col_wrap=4, 
    size=3, 
    scatter_kws=dict(alpha=.5),
)

g.set(xlim=(0, 1), ylim=(0, 1))
g.map_dataframe(add_rho)
g.set_xlabels('Prop X Missing')
g.set_ylabels('Prop A Missing')

# %%


# %% [markdown]
# ### Permutation Test

# %%
cell_ids = []
flags = []
for cell_id, dd in df.groupby('cell_id'):
    x_data = dd[dd.chrom == "X"].UMI.values
    x_data = x_data[x_data > 0]
    a_data = dd[dd.chrom == "X"].UMI.values
    a_data = a_data[a_data > 0]
    _, p_value = mannwhitneyu(x_data, a_data, alternative='less')
    
    if p_value <= 0.05:
        flags.append(True)
    else:
        flags.append(False)
        
    cell_ids.append(cell_id)

flag_x_lt_a = pd.Series(flags, index=cell_ids, name='flag_x_lt_a')

# %%
flag_w_cluster = pd.concat([flag_x_lt_a, read_clusters()], axis=1, sort=True)

# %%
flag_w_cluster.groupby('cluster').flag_x_lt_a.sum()

# %%


# %%


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

# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%

