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
# # Figures for BSC

# %% [markdown]
# Brian requested the following for the BSC:
#
# * heatmaps and clusters of scRNA.  
# * X, Y 4th chromosome expression patterns.  

# %%
import os
import sys
import re
from pathlib import Path
from itertools import zip_longest
from yaml import load

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    .to_frame()
    .assign(colors=lambda df: df.cluster.map(dict(zip(nbconfig.short_cluster_order, nbconfig.colors['clusters']))))
    .rename_axis('cell_id')
)

# %%
tsne = (
    nbconfig.seurat.get_tsne()
    .rename_axis('cell_id')
    .merge(clusters, on='cell_id')
)

# %%
def make_list(list_like):
    return np.array(
        list(
            zip_longest(list_like[:4], list_like[4:8], [list_like[-1]])
        )
    ).flatten().tolist()

# %%
fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(tsne.tSNE_1, tsne.tSNE_2, s=20, c=tsne.colors)

# clean up axis
ax.set_aspect('equal')
sns.despine(fig=fig, left=True, bottom=True)
plt.setp(ax, yticks=[], xticks=[]);

# legend
legend_elements = [
    #Patch(facecolor=color, edgecolor='k', label=f'{lclus} ({sclus})')
    Line2D([0], [0], marker='o', color=(1, 1, 1, 0), markeredgecolor=color, markerfacecolor=color, markersize=10, label=f'{lclus} ({sclus})')
    for sclus, lclus, color in zip(make_list(nbconfig.short_cluster_order),  make_list(nbconfig.cluster_order[:9]), make_list(nbconfig.colors['clusters'][:9]))
    if sclus is not None
]

ax.legend(handles=legend_elements, loc='lower center', ncol=4, bbox_to_anchor=[0.5, 1], facecolor=None)
for clus, row in tsne.groupby('cluster').agg({'tSNE_1': np.mean, 'tSNE_2': np.mean}).iterrows():
    plt.text(row.tSNE_1, row.tSNE_2, clus, backgroundcolor=(1, 1, 1, .9), ha='center', va='center')
plt.tight_layout()
plt.savefig('../output/notebook/2019-02-11_tsne.png')

# %%


# %%
zscores = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)

# %%
with open('../science_submission/config.yaml') as fh:
    lit_genes = load(fh.read())['lit_genes']

# %%
lit_fbgn = list(map(lambda x: nbconfig.symbol2fbgn[x], lit_genes))
lit_zscores = zscores.reindex(lit_fbgn).rename(index=nbconfig.fbgn2symbol)

# %%
lit_zscores

# %%
long_to_short = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
def process_text(txt):
    match = re.match(f'(?P<type>.*?)-(?P<rep>rep\d)', txt)
    if match['rep'] == 'rep2':
        return long_to_short[match['type']]
    return ''

# %%
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(lit_zscores, cmap='viridis', yticklabels=True, xticklabels=True, vmin=-3, vmax=3, cbar_kws=dict(label='Normalized Expression\n(z-score)'), ax=ax)

# fix up x-axis
labels = [
    process_text(l.get_text())
    for l in ax.get_xticklabels()
]
ax.set_xticklabels(labels, rotation=0, fontdict=dict(size=18), ha='center', va='bottom');
ax.set_xlabel('')
ax.xaxis.tick_top()

# fix up y-axis
labels = [
    l.get_text()
    for l in ax.get_yticklabels()
]
ax.set_yticklabels(labels, rotation=0, fontdict=dict(style='italic', size=18), va='center');
ax.set_ylabel('')

# Add cluster lines
loc = 3
for i in range(8):
    ax.axvline(loc, color='w', ls=':', lw=2)
    loc += 3
    
# Add cluster lines
ax.axhline(2, color='w', ls=':', lw=2)
ax.axhline(4, color='w', ls=':', lw=2)
ax.axhline(7, color='w', ls=':', lw=2)
ax.axhline(9, color='w', ls=':', lw=2)
ax.axhline(10, color='w', ls=':', lw=2)

# increase cbar axis
cbar = ax.collections[0].colorbar
label = cbar.ax.get_ylabel()
cbar.ax.set_ylabel(label, fontdict=dict(fontsize=18))
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('../output/notebook/2019-02-11_lit_genes.png')

# %%


# %%


# %%
from sklearn.cluster import KMeans

# %%
tpm = pd.read_parquet('../output/scrnaseq-wf/tpm_w_rep.parquet').loc[:, nbconfig.sel_cluster_order_w_rep]

# %%
X = tpm.values

# %%
kmeans = KMeans(n_clusters=9)

# %%
gene_clusters = kmeans.fit_predict(X)

# %%
zscores_kmeans = zscores.iloc[np.argsort(gene_clusters), :]

# %%
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(zscores_kmeans, cmap='viridis', yticklabels=False, xticklabels=True, vmin=-3, vmax=3, cbar_kws=dict(label='Normalized Expression\n(z-score)'), ax=ax)

# %%


# %%
# fix up x-axis
labels = [
    process_text(l.get_text())
    for l in ax.get_xticklabels()
]
ax.set_xticklabels(labels, rotation=0, fontdict=dict(size=18), ha='center', va='bottom');
ax.set_xlabel('')
ax.xaxis.tick_top()

# fix up y-axis
labels = [
    l.get_text()
    for l in ax.get_yticklabels()
]
ax.set_yticklabels(labels, rotation=0, fontdict=dict(style='italic', size=18), va='center');
ax.set_ylabel('')

# Add cluster lines
loc = 3
for i in range(8):
    ax.axvline(loc, color='w', ls=':', lw=2)
    loc += 3
    
# Add cluster lines
ax.axhline(2, color='w', ls=':', lw=2)
ax.axhline(4, color='w', ls=':', lw=2)
ax.axhline(7, color='w', ls=':', lw=2)
ax.axhline(9, color='w', ls=':', lw=2)
ax.axhline(10, color='w', ls=':', lw=2)

# increase cbar axis
cbar = ax.collections[0].colorbar
label = cbar.ax.get_ylabel()
cbar.ax.set_ylabel(label, fontdict=dict(fontsize=18))
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('../output/notebook/2019-02-11_all_genes.png')

# %%


# %%


# %%

