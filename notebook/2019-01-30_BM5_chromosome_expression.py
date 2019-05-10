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

# %% [markdown]
# # BM5 Chromosome Expression

# %% [markdown]
# In this project we hypothesize that the X chromosome is getting silenced prior to the autosomes. Sharvani performed a translocation experiment where she has swapped part of the X with the 2.
#
# **Do we see an increase in X expression and a decrease in 2 expression in the translocation lines?**
#
# BM5 is a sterile translocation, so we expect that the break point is near the centromere causing movement of a large chunk of the chromosome. We also have the irratiated parental stock for use as a control.

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
nbconfig = Nb.setup_notebook(seurat_dir='../output/translocations-wf/translocation_BM5_force')

# Modify nbconfig with the BM5's annotations
nbconfig.short_cluster_annot = {
    0: 'LC',
    1: 'SP',
    2: 'EC',
    3: 'LS',
    4: 'TE',
    5: 'SP',
    6: 'SP',
    7: 'MS',
    8: 'PC',
    9: 'ES',
    10: 'PC',
}

nbconfig.short_cluster_order = [
    'SP',
    'ES',
    'MS',
    'LS',
    'EC',
    'LC',
    'TE',
    'PC',
]

nbconfig.colors['clusters'] = [
    # Germline
    (0.6943944636678201, 0.07003460207612457, 0.09231833910034601),
    (0.8901960784313725, 0.18562091503267975, 0.15294117647058825),
    (0.9843752402921953, 0.4181468665897732, 0.2926566705113418),
    (0.9935870818915802, 0.8323414071510957, 0.7624913494809689),
    # soma
    (0.06251441753171857, 0.35750865051903113, 0.6429065743944637),
    (0.42274509803921567, 0.684075355632449, 0.8398923490965013),
    (0.41568627450980394, 0.23921568627450981, 0.6039215686274509),
    (0.6941176470588235, 0.34901960784313724, 0.1568627450980392),
]

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

def read_data(rep2=False, tpm=False):
    fbgn2chrom = read_fbgn2chrom()
    clusters = read_clusters()
    data = nbconfig.seurat.get_raw()
    value_name = 'UMI'
    
    # Munge together
    melted = data.reset_index().melt(id_vars='FBgn', value_name=value_name)
    return melted.join(clusters, on='cell_id').join(fbgn2chrom, on='FBgn')

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
g.fig.suptitle('BM5', fontsize=10);

# %%


