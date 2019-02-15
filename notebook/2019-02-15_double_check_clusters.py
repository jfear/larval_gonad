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
# # Double Check Clustering

# %% [markdown]
# I am still troubled by LS cluster. I just want to do some sanity checks.

# %%
import os
import sys
import re
from pathlib import Path
from itertools import combinations

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
clusters = (
    nbconfig.seurat.get_clusters('res.0.6')
    .map(nbconfig.short_cluster_annot)
    .pipe(lambda x: x[x != 'UNK'])
    .astype('category')
    .cat.as_ordered()
    .cat.reorder_categories(nbconfig.short_cluster_order)
    .rename_axis('cell_id')
)

# %%
raw = nbconfig.seurat.get_raw()

# %%
umi_by_cell = (
    raw.sum()
    .rename('UMI')
)

# %%
biomarkers = (
    nbconfig.seurat.get_biomarkers('res.0.6')
    .assign(cluster = lambda df: df.cluster.map(nbconfig.short_cluster_annot))
    .query('cluster != "UNK"')
    .assign(cluster = lambda df: df.cluster.astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order))
)

biomarkers = biomarkers.join(biomarkers.groupby('FBgn').size().rename('num_clusters'))

# %%
# Genes commented out are not present int he zscores dataset
lit_genes = [
    #GSC, spermatogonia, early spermatocytes [:12] (12) (7)
    'vas',
    'bam',
    'Phf7',
    'CG11697',
    'p53',
    #'nos',
    #'bgcn',
    #'tut',
    'Rbp9',
    'peb',
    #'tej',
    #'Marf',
    # Later spermatocytes and spermatids [12:34] (22) (18)
    'aly',
    'nht',
    'soti',
    'dj',
    'ocn',
    'can',
    'fzo',
    'bol',
    #'mle',
    #'mia',
    'CG3927',
    'sunz',
    'sowi',
    'd-cup',
    'c-cup',
    'wa-cup',
    #'p-cup',
    #'r-cup',
    'oys',
    'topi',
    'sa',
    'CG8368',
    # Enriched in CySC lineage [34:58] (24) (18)
    'tj',
    #'eya',
    'zfh1',
    'vn',
    'foxo',
    #'apt',
    'ImpL2',
    'Wnt4',
    'Nrt',
    'bnb',
    #'neur',
    'robo2',
    'EcR',
    'gbb',
    'spict',
    'puc',
    #'sev',
    'hui',
    #'sano',
    'glob1',
    'Eip93F',
    'fax',
    'kek1',
    #'so',
    # Terminal epithelia [58:67] (9) (8)
    'nord',
    'retn',
    'abd-A',
    'Abd-B',
    'Wnt2',
    'Six4',
    #'CG18628',
    'MtnA',
    'N',
    # Pigment cells [67:] (4)
    'vkg',
    'Sox100B',
    'bw',
    'ems',
]

lit_fbgns = list(map(lambda x: nbconfig.symbol2fbgn[x], lit_genes))


# %% [markdown]
# ## Number of UMI by Clusters

# %% [markdown]
# Taking a quick look at the UMI by cluster we see that late spermatocytes have the lowest read counts. Also the early spermatocytes has a really large variation.

# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=plt.figaspect(1/2))
_dat = pd.concat([umi_by_cell, clusters], axis=1, sort=True, join='inner').groupby('cluster').sum().reset_index()
sns.barplot('cluster', 'UMI', data=_dat, palette=nbconfig.colors['clusters'], ax=ax1)
ax1.set_title('Number of UMI Per Cluster')

_dat = pd.concat([umi_by_cell, clusters], axis=1, sort=True, join='inner')
sns.boxplot('cluster', 'UMI', data=_dat, showfliers=False, palette=nbconfig.colors['clusters'], ax=ax2)
ax2.set_title('Distribution of UMI Per Cluster');

# %% [markdown]
# ## Early Spermatocytes

# %%
gonia_singles = biomarkers.query('cluster == "SP" & num_clusters == 1')
gonia_singles.reindex(lit_fbgns).dropna()

# %%
early_singles = biomarkers.query('cluster == "ES" & num_clusters == 1')
early_singles.reindex(lit_fbgns).dropna()

# %%
mid_singles = biomarkers.query('cluster == "MS" & num_clusters == 1')
mid_singles.reindex(lit_fbgns).dropna()

# %%
late_singles = biomarkers.query('cluster == "LS" & num_clusters == 1')
late_singles.reindex(lit_fbgns).dropna()

# %%
late_singles

# %%


# %%
gonia_multi = biomarkers.query('cluster == "SP" & num_clusters > 1')
gonia_multi.reindex(lit_fbgns).dropna()

# %%
early_multi = biomarkers.query('cluster == "ES" & num_clusters > 1')
early_multi.reindex(lit_fbgns).dropna()

# %%
mid_multi = biomarkers.query('cluster == "MS" & num_clusters > 1')
mid_multi.reindex(lit_fbgns).dropna()

# %%
late_multi = biomarkers.query('cluster == "LS" & num_clusters > 1')
late_multi.reindex(lit_fbgns).dropna()

# %%


# %%
res = []
for gene, dd in biomarkers.query('num_clusters > 1').groupby('FBgn'):
    res.append((gene, '|'.join(dd.cluster.sort_values().values)))

multi_genes = pd.DataFrame(res, columns=['FBgn', 'clusters']).set_index('FBgn').clusters

# %%
for c1, c2 in combinations(nbconfig.short_cluster_order, 2):
    print(f'{c1}|{c2}', (multi_genes == f'{c1}|{c2}').sum())

# %%
for c1, c2, c3 in combinations(nbconfig.short_cluster_order, 3):
    print(f'{c1}|{c2}|{c3}', (multi_genes == f'{c1}|{c2}|{c3}').sum())

# %%


# %%


# %%


# %%


# %%

