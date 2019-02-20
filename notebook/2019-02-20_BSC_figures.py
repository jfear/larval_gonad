# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.0
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
from collections import defaultdict

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, dendrogram

# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
sns.set_context('poster')

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

# %% [markdown]
# ## TSNE

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
plt.savefig('../output/notebook/2019-02-20_tsne.svg', bbox_inches='tight')
del tsne 
del legend_elements

# %% [markdown]
# ## Lit Genes (Short)

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
long_to_short = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
def process_text(txt):
    match = re.match(f'(?P<type>.*?)-(?P<rep>rep\d)', txt)
    if match['rep'] == 'rep2':
        return long_to_short[match['type']]
    return ''

# %%
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(lit_zscores, cmap='viridis', yticklabels=True, xticklabels=True, vmin=-3, vmax=3, cbar_kws=dict(label='Normalized Expression\n(z-score)'), ax=ax, rasterized=True)

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
plt.savefig('../output/notebook/2019-02-20_lit_genes.svg', bbox_inches='tight')
del zscores

# %% [markdown]
# ## Lit Genes (full)

# %%
zscores = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)

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


# %%
lit_fbgn = list(map(lambda x: nbconfig.symbol2fbgn[x], lit_genes))
lit_zscores = zscores.reindex(lit_fbgn).dropna().rename(index=nbconfig.fbgn2symbol)

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
sns.heatmap(lit_zscores, cmap='viridis', yticklabels=True, xticklabels=True, vmin=-3, vmax=3, cbar_kws=dict(label='Normalized Expression\n(z-score)'), ax=ax, rasterized=True)

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
ax.set_yticklabels(labels, rotation=0, fontdict=dict(style='italic', size=8), va='center');
ax.set_ylabel('')

# Add cluster lines
loc = 3
for i in range(8):
    ax.axvline(loc, color='w', ls=':', lw=2)
    loc += 3
    
# Add cluster lines
previous = 0
for annot, loc in zip(['Early\nGerm Cells', 'Late\nGerm Cells', 'CySC', 'TE', 'PC'], np.cumsum([7, 18, 18, 8, 4])):
    if annot != 'PC':
        ax.axhline(loc, color='w', ls=':', lw=2)
    text_loc = loc - np.ceil((loc - previous) / 2)
    plt.text(-6, text_loc, annot, ha='center', va='center', fontsize=14, multialignment='center')
    previous = loc

# increase cbar axis
cbar = ax.collections[0].colorbar
label = cbar.ax.get_ylabel()
cbar.ax.set_ylabel(label, fontdict=dict(fontsize=18))
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('../output/notebook/2019-02-20_all_lit_genes.svg', bbox_inches='tight')
del zscores

# %% [markdown]
# ## Unique BioMarkers

# %%
long_to_short = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
def process_text(txt):
    match = re.match(f'(?P<type>.*?)-(?P<rep>rep\d)', txt)
    if match['rep'] == 'rep2':
        return long_to_short[match['type']]
    return ''

# %%
zscores = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)

# %%
biomarkers = (
    nbconfig.seurat.get_biomarkers('res.0.6')
    .cluster.map(nbconfig.short_cluster_annot)
    .pipe(lambda x: x[x != "UNK"])
    .astype('category')
    .cat.as_ordered()
    .cat.reorder_categories(nbconfig.short_cluster_order)
)

unique_fbgns = (
    biomarkers
    .groupby('FBgn').size()
    .pipe(lambda x: x[x == 1])
    .index
)


# %%
biomarkers_by_cluster = biomarkers[biomarkers.index.isin(unique_fbgns)].sort_values()

# %%
zscores = zscores.reindex(biomarkers_by_cluster.index)

# %%
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(zscores, cmap='viridis', yticklabels=False, xticklabels=True, vmin=-3, vmax=3, cbar_kws=dict(label='Normalized Expression\n(z-score)'), ax=ax, rasterized=True)

# fix up x-axis
labels = [
    process_text(l.get_text())
    for l in ax.get_xticklabels()
]
ax.set_xticklabels(labels, rotation=0, fontdict=dict(size=18), ha='center', va='bottom');
ax.set_xlabel('')
ax.xaxis.tick_top()

# fix up y-axis
ax.set_ylabel('Gene Cluster (Unique Bio Marker)', labelpad=25)

# Add cluster lines
loc = 3
for i in range(8):
    ax.axvline(loc, color='w', ls=':', lw=2)
    loc += 3
    
# Add cluster lines
locs = biomarkers_by_cluster.value_counts().sort_index().values

previous = 0
for clus, loc in zip(nbconfig.short_cluster_order, np.cumsum(locs)):
    if clus != 'PC':
        ax.axhline(loc, color='w', ls=':', lw=2)
        
    text_loc = loc - np.ceil((loc - previous) / 2)
    plt.text(-.2, text_loc, clus, ha='right', va='center', fontsize=14)
    previous = loc

# increase cbar axis
cbar = ax.collections[0].colorbar
label = cbar.ax.get_ylabel()
cbar.ax.set_ylabel(label, fontdict=dict(fontsize=18))
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('../output/notebook/2019-02-20_unique_biomakrer_genes.svg', bbox_inches='tight')
del zscores
del biomarkers

# %% [markdown]
# ## Non-unique BioMarkers

# %%
long_to_short = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
def process_text(txt):
    match = re.match(f'(?P<type>.*?)-(?P<rep>rep\d)', txt)
    if match['rep'] == 'rep2':
        return long_to_short[match['type']]
    return ''

# %%
zscores = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)

# %%
biomarkers = (
    nbconfig.seurat.get_biomarkers('res.0.6')
    .cluster.map(nbconfig.short_cluster_annot)
    .pipe(lambda x: x[x != "UNK"])
    .astype('category')
    .cat.as_ordered()
    .cat.reorder_categories(nbconfig.short_cluster_order)
)

mapper = {
    'SP': 'Germ',
    'ES': 'Germ',
    'MS': 'Germ',
    'LS': 'Germ',
    'EC': 'Soma',
    'MC': 'Soma',
    'LC': 'Soma',
    'TE': 'Soma',
    'PC': 'Soma',
}

multi_fbgns = (
    (biomarkers.map(mapper).groupby('FBgn').value_counts() > 0)
    .unstack()
    .fillna(False)
    .assign(germ_only = lambda df: df.Germ & ~df.Soma)
    .assign(soma_only = lambda df: ~df.Germ & df.Soma)
    .assign(germ_and_soma = lambda df: df.Germ & df.Soma)
    .drop(['Germ', 'Soma'], axis=1)
    .sort_values(by='soma_only', ascending=True)
    .sort_values(by='germ_and_soma', ascending=False)
    .sort_values(by='germ_only', ascending=False)
)

# %%
ordered_genes = []
for name, dd in multi_fbgns.idxmax(axis=1).rename('names').astype('category').cat.as_ordered().cat.reorder_categories(['germ_only', 'germ_and_soma', 'soma_only']).to_frame().groupby('names'):
    _curr = zscores.reindex(dd.index)
    tree = linkage(_curr, method='average')
    leaves = dendrogram(tree, no_plot=True)['leaves']
    ordered_genes.extend(_curr.index[leaves].tolist())

# %%
zscores = zscores.reindex(ordered_genes)

# %%
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(zscores, cmap='viridis', yticklabels=False, xticklabels=True, vmin=-3, vmax=3, cbar_kws=dict(label='Normalized Expression\n(z-score)'), ax=ax, rasterized=True)

# fix up x-axis
labels = [
    process_text(l.get_text())
    for l in ax.get_xticklabels()
]
ax.set_xticklabels(labels, rotation=0, fontdict=dict(size=18), ha='center', va='bottom');
ax.set_xlabel('')
ax.xaxis.tick_top()

# fix up y-axis
ax.set_ylabel('')

# Add cluster lines
loc = 3
for i in range(8):
    ax.axvline(loc, color='w', ls=':', lw=2)
    loc += 3
    
# Add cluster lines
locs = multi_fbgns.sum()[['germ_only', 'germ_and_soma', 'soma_only']].values

previous = 0
for clus, loc in zip(['Germ Line\nOnly', 'Germline and\nSoma', 'Soma\nOnly'], np.cumsum(locs)):
    if clus != 'soma_only':
        ax.axhline(loc, color='w', ls=':', lw=2)
        
    text_loc = loc - np.ceil((loc - previous) / 2)
    plt.text(-3.5, text_loc, clus, ha='center', va='center', fontsize=14, multialignment='center')
    previous = loc

# increase cbar axis
cbar = ax.collections[0].colorbar
label = cbar.ax.get_ylabel()
cbar.ax.set_ylabel(label, fontdict=dict(fontsize=18))
cbar.ax.tick_params(labelsize=14)

# save figure
plt.savefig('../output/notebook/2019-02-20_multi_biomakrer_genes.svg', bbox_inches='tight')
del zscores
del biomarkers

# %%
germ_and_soma = multi_fbgns[multi_fbgns.germ_and_soma].index.tolist()

# %% [markdown] {"toc-hr-collapsed": false}
# ## X Y 4th Expression by Rep

# %%
num_genes_per_chrom = nbconfig.fbgn2chrom.groupby('chrom').size()
num_genes_per_chrom

# %% [markdown]
# ### TPM normalized by chrom by cluster

# %%
tpm = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)
num_genes_per_chrom = nbconfig.fbgn2chrom.groupby('chrom').size()

# %%
mapper = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
tpm_tidy = (
    tpm
    .merge(nbconfig.fbgn2chrom, on='FBgn')
    .reset_index()
    .melt(id_vars=['FBgn', 'chrom'], value_name='TPM', var_name='name')
    # split name up into cluster and replicate
    .assign(rep = lambda df: df.name.str.extract('(rep\d)').values)
    .assign(cluster = lambda df: df.name.str.extract('(?P<cluster>.*?)-rep\d').cluster.map(mapper).astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order))
    .drop('name', axis=1)
    # Sum gene counts to the chromosome arm level and normalize by the number of genes on that arm
    .groupby(['chrom', 'cluster', 'rep']).TPM.sum()
    .div(num_genes_per_chrom, level=0)
    # clean up 
    .to_frame()
    .reset_index()
    .rename({0: 'TPM'}, axis=1)
    # Remove chr fo nicer plotting
    .assign(chrom = lambda df: df.chrom.str.extract('chr(.*)').values)
)

# %%
g = (
    sns.FacetGrid(tpm_tidy, col='cluster', col_order=nbconfig.short_cluster_order, col_wrap=4, hue='cluster', palette=nbconfig.colors['clusters'])
    .map(sns.barplot, 'chrom', 'TPM', order=[x.lstrip('chr') for x in nbconfig.chrom_order], errwidth=1.8, capsize=.3, linewidth=1, edgecolor='k')
    .set_titles('{col_name}', size=18, y=.9)
    .set_ylabels('')
    .set_xlabels('Chromosome')
    .despine(left=True)
)
plt.text(0.03, 0.5, 'Normalized Chromosomal Count', rotation=90, ha='right', va='center', transform=g.fig.transFigure)
plt.subplots_adjust(hspace=0.08, wspace=0.08)
g.fig.savefig('../output/notebook/2019-02-20_bar_plot_tpm_expression_by_arm.svg', bbox_inches='tight')
del tpm
del tpm_tidy

# %% [markdown]
# ### Log TPM normalized by chrom by cluster

# %%
tpm = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)
num_genes_per_chrom = nbconfig.fbgn2chrom.groupby('chrom').size()

# %%
mapper = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
tpm_tidy = (
    tpm
    .merge(nbconfig.fbgn2chrom, on='FBgn')
    .reset_index()
    .melt(id_vars=['FBgn', 'chrom'], value_name='TPM', var_name='name')
    # split name up into cluster and replicate
    .assign(rep = lambda df: df.name.str.extract('(rep\d)').values)
    .assign(cluster = lambda df: df.name.str.extract('(?P<cluster>.*?)-rep\d').cluster.map(mapper).astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order))
    .drop('name', axis=1)
    # Sum gene counts to the chromosome arm level and normalize by the number of genes on that arm
    .groupby(['chrom', 'cluster', 'rep']).TPM.sum()
    .div(num_genes_per_chrom, level=0)
    # clean up 
    .to_frame()
    .reset_index()
    .rename({0: 'TPM'}, axis=1)
    # Remove chr fo nicer plotting
    .assign(chrom = lambda df: df.chrom.str.extract('chr(.*)').values)
    .assign(log_TPM = lambda df: np.log2(df.TPM + .7))
)

# %%
g = (
    sns.FacetGrid(tpm_tidy, col='cluster', col_order=nbconfig.short_cluster_order, col_wrap=4, hue='cluster', palette=nbconfig.colors['clusters'], ylim=(0, 8))
    .map(sns.barplot, 'chrom', 'log_TPM', order=[x.lstrip('chr') for x in nbconfig.chrom_order], errwidth=1.8, capsize=.3, linewidth=1, edgecolor='k')
    .set_titles('{col_name}', size=18, y=.9)
    .set_ylabels('')
    .set_xlabels('Chromosome')
    .despine(left=True)
)
plt.text(0.03, 0.5, 'log2(Normalized Chromosomal Count / number genes)', rotation=90, ha='right', va='center', transform=g.fig.transFigure)
plt.subplots_adjust(hspace=0.08, wspace=0.08)
g.fig.savefig('../output/notebook/2019-02-20_bar_plot_log2_tpm_expression_by_arm.svg', bbox_inches='tight')
del tpm
del tpm_tidy

# %% [markdown]
# ### Proportion of genes on

# %%
raw = (
    pd.read_parquet('../output/scrnaseq-wf/raw_by_cluster_w_rep.parquet')
    .loc[:, nbconfig.sel_cluster_order_w_rep]
)
num_genes_per_chrom = nbconfig.fbgn2chrom.groupby('chrom').size()

# %%
mapper = dict(zip(nbconfig.sel_cluster_order, nbconfig.short_cluster_order))

# %%
prop_tidy = (
    (raw > 5)
    .merge(nbconfig.fbgn2chrom, on='FBgn')
    .reset_index()
    .melt(id_vars=['FBgn', 'chrom'], value_name='flag_on', var_name='name')
    # split name up into cluster and replicate
    .assign(rep = lambda df: df.name.str.extract('(rep\d)').values)
    .assign(cluster = lambda df: df.name.str.extract('(?P<cluster>.*?)-rep\d').cluster.map(mapper).astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order))
    .drop('name', axis=1)
    # Sum gene counts to the chromosome arm level and normalize by the number of genes on that arm
    .groupby(['chrom', 'cluster', 'rep']).flag_on.sum()
    .div(num_genes_per_chrom, level=0)
    .mul(100)
    # clean up 
    .to_frame()
    .reset_index()
    .rename({0: 'Percent On'}, axis=1)
    # Remove chr fo nicer plotting
    .assign(chrom = lambda df: df.chrom.str.extract('chr(.*)').values)
)

# %%
g = (
    sns.FacetGrid(prop_tidy, col='cluster', col_order=nbconfig.short_cluster_order, col_wrap=4, hue='cluster', palette=nbconfig.colors['clusters'])
    .map(sns.barplot, 'chrom', 'Percent On', order=[x.lstrip('chr') for x in nbconfig.chrom_order], errwidth=1.8, capsize=.3, linewidth=1, edgecolor='k')
    .set_titles('{col_name}', size=18, y=.9)
    .set_ylabels('')
    .set_xlabels('Chromosome')
    .despine(left=True)
)
plt.text(0.03, 0.5, 'Percent Genes Expressed', rotation=90, ha='right', va='center', transform=g.fig.transFigure)
plt.subplots_adjust(hspace=0.08, wspace=0.08)
g.fig.savefig('../output/notebook/2019-02-20_bar_plot_percent_expressed_by_arm.svg', bbox_inches='tight')
del raw
del prop_tidy

# %% [markdown]
# ## X Y 4th to Autosome Ratios

# %%
mapper = {'chrX': 'X', 'chrY': 'Y', 'chr4': '4', 'chr2L': 'A','chr2R': 'A','chr3L': 'A','chr3R': 'A',}
fbgn2chrom = (
    pd.read_csv('../output/fbgn2chrom.tsv', sep='\t', index_col=0)
    .chrom.map(mapper)
    .dropna()
)

num_genes_per_chrom = fbgn2chrom.value_counts()
num_genes_per_chrom

# %%
raw_tidy = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet', columns=clusters.index)
    .join(fbgn2chrom, on='FBgn')
    .reset_index()
    .melt(id_vars=['FBgn', 'chrom'], value_name='UMI', var_name='cell_id')
    .join(clusters, on='cell_id')
    .drop('colors', axis=1)
)

# %%
ratios_by_cell = (
    raw_tidy.groupby(['cell_id', 'cluster', 'chrom'])
    .UMI.sum()
    .div(num_genes_per_chrom / 1e3, level='chrom')
    .unstack()
    .assign(ratio_x = lambda df: df.X / df.A)
    .assign(ratio_y = lambda df: df.Y / df.A)
    .assign(ratio_4 = lambda df: df['4'] / df.A)
    .drop(['X', 'Y', 'A', '4'], axis=1)
    .reset_index('cluster')
)

# %%
median_ratios = ratios_by_cell.groupby('cluster').median().reindex(nbconfig.short_cluster_order) + 0.0001

# %%
median_ratios

# %%
pd.concat([
    (median_ratios.loc['SP'] / median_ratios.loc['ES']).rename('SP/ES'), 
    (median_ratios.loc['SP'] / median_ratios.loc['MS']).rename('SP/MS'),
    (median_ratios.loc['SP'] / median_ratios.loc['LS']).rename('SP/LS')
], axis=1, sort=True)

# %%
1 / median_ratios

# %%
cutoff = 0.05

results = []
permuted_ratios_by_cell = ratios_by_cell.copy()
for iteration in range(10_000):
    permuted_ratios_by_cell.cluster = permuted_ratios_by_cell.cluster.sample(frac=1).values
    for clus, observed_ratios in ratios_by_cell.groupby('cluster'):
        permuted_ratios = permuted_ratios_by_cell.query(f'cluster == "{clus}"')
        _, pval_x = mannwhitneyu(observed_ratios.ratio_x, permuted_ratios.ratio_x, alternative='less')
        _, pval_y = mannwhitneyu(observed_ratios.ratio_y, permuted_ratios.ratio_y, alternative='greater')
        _, pval_4 = mannwhitneyu(observed_ratios.ratio_4, permuted_ratios.ratio_4, alternative='less')
        results.append((clus, pval_x <= cutoff, pval_y <= cutoff, pval_4 <= cutoff))

# %%
pvals = 1 - (
    pd.DataFrame(results, columns=['cluster', 'sig_x', 'sig_y', 'sig_4']).groupby('cluster')
    .mean()
    .rename(columns=dict(sig_x='pval_x_lt_a', sig_y='pval_y_gt_a', sig_4='pval_4_lt_a'))
    .loc[nbconfig.short_cluster_order, :]
)
pvals

# %%
def whisker(dat):
    low, high = np.percentile(dat, [25, 75])
    iqr = high - low
    return high + (1.5 * iqr)

# %%
def plot_pval(dat, pvals, ax):
    whiskers = (
        dat.groupby('cluster')
        .apply(whisker)
        .to_dict()
    )
    
    for i, clus in enumerate(nbconfig.short_cluster_order):
        pval = pvals.loc[clus]
        loc = whiskers[clus]
        
        if pval <= 0.001:
            ax.text(i, loc, '***', ha='center', va='bottom')
        elif pval <= 0.01:
            ax.text(i, loc, '**', ha='center', va='bottom')
        elif pval <= 0.05:
            ax.text(i, loc, '*', ha='center', va='bottom')

# %%
ratios_by_cell.head()

# %%
fig, axes = plt.subplots(3, 1, figsize=plt.figaspect(3), sharex=True, gridspec_kw=dict(hspace=.2))

_defaults = dict(x='cluster', data=ratios_by_cell, order=nbconfig.short_cluster_order, showfliers=False, palette=nbconfig.colors['clusters'], notch=True)

ax = axes[0]
sns.boxplot(y='ratio_x', ax=ax, **_defaults)
plot_pval(ratios_by_cell[['cluster', 'ratio_x']], pvals.pval_x_lt_a, ax)
sns.despine(ax=ax, left=True)
ax.set_title('X Chromsome Expression')
ax.set_xlabel('')
ax.set_ylabel('X / Autosome')
ax.axhline(1, color='gray', ls='--', zorder=0, alpha=.8, label='Median A')

ax = axes[1]
sns.boxplot(y='ratio_4', ax=ax, **_defaults)
plot_pval(ratios_by_cell[['cluster', 'ratio_4']], pvals.pval_4_lt_a, ax)
sns.despine(ax=ax, left=True)
ax.set_title('4th Chromsome Expression')
ax.set_xlabel('')
ax.set_ylabel('4 / Autosome')
ax.axhline(1, color='gray', ls='--', zorder=0, alpha=.8, label='Median A')

ax = axes[2]
sns.boxplot(y='ratio_y', ax=ax, **_defaults)
plot_pval(ratios_by_cell[['cluster', 'ratio_y']], pvals.pval_y_gt_a, ax)
sns.despine(ax=ax, left=True)
ax.set_title('Y Chromsome Expression')
ax.set_xlabel('')
ax.set_ylabel('Y / Autosome')
ax.set_ylim(None, .3)

fig.savefig('../output/notebook/2019-02-20_boxplot_autosome_ratios.svg', bbox_inches='tight')

# %%
del raw_tidy


# %%
del ratios_by_cell


# %%
