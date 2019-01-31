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

# %% [markdown] {"toc-hr-collapsed": false}
# # Generate norm counts barplot for comparison with translocations.

# %% [markdown]
# I want this plot generated for w1118 samples for comparison with translocations.

# %%
import os
import sys
import re
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
from scipy.cluster.hierarchy import linkage, dendrogram
import statsmodels.formula.api as smf

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
        'chr2L': '2L',
        'chr2R': '2R',
        'chr3L': '3L',
        'chr3R': '3R',
    }

    fbgn2chrom = (pd.read_csv('../output/fbgn2chrom.tsv', sep='\t', index_col=0)
                      .query('chrom != "chrM"')
                      .chrom.map(mapper)
                      .astype('category')
                      .cat.as_ordered()
                 )
    
    return fbgn2chrom.cat.reorder_categories(['X', '2L', '2R', '3L', '3R', 'Y', '4'])


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
    return melted.join(clusters, on='cell_id').join(fbgn2chrom, on='FBgn').join(rep, on='cell_id').dropna()

# %% [markdown]
# ## Data Prep

# %%
df = read_data()
df['missing'] = (df.UMI == 0).values

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
        .cat.reorder_categories(['X', '2L', '2R', '3L', '3R', 'Y', '4'])
)

# %% [markdown] {"toc-hr-collapsed": true}
# ## Chromosome Expression

# %% [markdown]
# ### Cell level chromosome coverage

# %%
norm_cnts = norm_cnts.query('chrom == ["X", "2L", "2R"]').copy()
norm_cnts.chrom = norm_cnts.chrom.cat.remove_unused_categories()

# %%
g = sns.FacetGrid(norm_cnts, row='chrom', sharey=False, aspect=1.2)
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
g.fig.suptitle('w[1118]', fontsize=10);

# %%


