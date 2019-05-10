# -*- coding: utf-8 -*-
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
# # Identifying Cell Types of the Testis

# %%
import os
import sys
import re
from pathlib import Path
from yaml import load

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# %% [markdown]
# ## Data

# %%
# load all of my config settings
config = {}
with open('../config/common.yaml') as fh:
    config.update(load(fh.read()))
    
with open('../config/colors.yaml') as fh:
    config['colors'] = load(fh.read())
    
with open('../science_submission/config.yaml') as fh:
    config.update(load(fh.read()))

# %%
# Get list of genes from the literature
symbol2fbgn = pd.read_pickle('../output/science_submission/symbol2fbgn.pkl')
lit_genes = config['lit_genes_long']
lit_genes_fbgn = [symbol2fbgn[x] for x in config['lit_genes_long']]
print(len(lit_genes))
print(', '.join(sorted(lit_genes, key=lambda x: x.lower())))

# %%
# Get mappint of cell_id to short cluster name, remove unknown clusters
clusters = (
    pd.read_parquet('../output/scrnaseq-wf/clusters.parquet')
    .assign(cluster = lambda df: pd.Categorical(df.cluster.map(config['short_cluster_annot']), ordered=True, categories=config['short_cluster_order']))
    .dropna()
)
clusters.cluster.value_counts().sort_index().map(lambda x: f'{x:,}').rename('Cells Per Cluster').to_frame()

# %%
# Get biomarkers
resolution = config['resolution']
biomarkers = (
    pd.read_csv(f'../output/scrnaseq-wf/scrnaseq_combine_force/biomarkers_{resolution}.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .query('p_val_adj <= 0.05')
    .assign(cluster = lambda df: pd.Categorical(df.cluster.map(config['short_cluster_annot']), ordered=True, categories=config['short_cluster_order']))
    .dropna()
    .loc[:, ['gene_symbol', 'cluster']]
)
display(biomarkers.groupby('cluster').size().sort_index().map(lambda x: f'{x:,}').rename('Biomarkers per cluster').to_frame())
print(f'There are {biomarkers.index.unique().shape[0]:,} biomakrer genes')

# %%
print(f'Of the {len(lit_genes_fbgn):,} literature genes, only {biomarkers.query(f"FBgn == {lit_genes_fbgn}").index.unique().shape[0]:,} where in the biomarker list.')

# %%
# Figure out how many literature genes were called biomarker for the right cluster.
PASS = 0
TOTAL = 0
def get_lit_gene_subset(cell_type):
    if cell_type == 'gonia':
        name = 'Spermatogonia'
        idx = (0, 12)
    elif cell_type == 'primary':
        name = 'Primary Spermatocytes'
        idx = (12, 34)
    elif cell_type == 'cyst':
        name = 'Somatic Cyst Cells'
        idx = (34, 58)
    elif cell_type == 'te':
        name = 'Terminal Epithelium'
        idx = (58, 67)
    elif cell_type == 'pc':
        name = 'Pigment Cells'
        idx = (67,71)
    display(HTML(f'<h3>{name}</h3>'))
    print(', '.join(sorted(lit_genes[idx[0]:idx[1]], key=lambda x: x.lower())))
    return lit_genes_fbgn[idx[0]: idx[1]]

def check_biomarkers(fbgns, cell_type_pattern):
    subset = (
        biomarkers.query(f'FBgn == {fbgns}')
        .groupby('gene_symbol')
        .apply(lambda df: '|'.join(df.cluster.sort_values().values))
        .rename('clusters')
        .to_frame()
        .assign(lower=lambda df: df.index.str.lower())
        .sort_values(by='lower')
        .drop('lower', axis=1)
    )
    num_with_correct_cell_type = subset.clusters.str.contains(cell_type_pattern).sum()
    total_num_genes = subset.shape[0]
    
    global PASS
    global TOTAL
    PASS += num_with_correct_cell_type
    TOTAL += total_num_genes
    
    print(f'There were ({num_with_correct_cell_type:,} / {total_num_genes:,} = {num_with_correct_cell_type / total_num_genes * 100:.2f}%) literature genes that were called biomarker in the correct cluster.')
    display(subset)

fbgns = get_lit_gene_subset('gonia')
check_biomarkers(fbgns, 'SP')

fbgns = get_lit_gene_subset('primary')
check_biomarkers(fbgns, 'E1º|M1º|L1º')

fbgns = get_lit_gene_subset('cyst')
check_biomarkers(fbgns, 'EC|MC|LC')

fbgns = get_lit_gene_subset('te')
check_biomarkers(fbgns, 'TE')

fbgns = get_lit_gene_subset('pc')
check_biomarkers(fbgns, 'PC')

print(f'{PASS:,} / {TOTAL:,} = {PASS / TOTAL:,.2f}%')

# %%
