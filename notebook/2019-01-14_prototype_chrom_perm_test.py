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
# # Prototype Chromosome Permutation Test

# %% [markdown]
# I need to work out what permutation tests to do. There are several ways we can compare X:A or Y:A or 4:A. 
#
# These data are challenging because at the individual gene cell level the data are very sparse. One solution to this problem is to aggregate data to the cell type level. Unfortunately, the missingness is not completely random and shows distinct patterns for different cell types. Making me worry that aggregating to the cell type level would incorporate bias differently for each cell type. To combat this problem we are taking a permutation approach, where we can randomly sample cells and calculate a statistic that should help us capture trends inspite of the missingness. 
#
# For looking at X:A expression, I need to decide what measure to use along with what statistic. 
#
# ## Possible Measures
#
# There are a number of ways to look at X:A. Here are a few possible measures that we could use.
#
# * Compare read counts per chromosome normalized by chromosome length. (It may be hard to get chromosome counts per cell)
# * Compare read counts mapped to genic regions normalized by gene count (or total gene length).
# * Compare the proportion of genes "expressed".
# * Compare read count mapping to specific genes (i.e., housekeeping genes).
#
# The X chromosome is rather comparable in size and gene content to 2L, 2R, 3L, and 3R. I think the first three measures should all behave similarly. However, ultimately we also want to look at the behavior of Y and 4th in relation to the autosomes. These chromosomes have large differences in size and gene content, so I don't know if that will affect the usefullness of these three measures. 
#
# The last measure conceptually sounds nice, housekeeping genes are throught to behave similarly across cell types. Any change is housekeeping expression would give strong evidence. However, annotation of what is a "housekeeping" gene is challenging, and more importantly these genes are not likely to be equally distributed across chromosomes. 

# %%
import os
import sys
import re
from pathlib import Path
from itertools import combinations
import re

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb
from larval_gonad.stats import permutation_test_chrom1_lt_chrom2

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
# Shortten cluster names for nicer plots
cluster_annot = {
    0: 'LS',
    1: 'MC',
    2: 'MS',
    3: 'ES',
    4: 'LC',
    5: 'EC',
    6: 'SP',
    7: 'TE',
    8: 'PC',
}

cluster_order = ['SP', 'ES', 'MS', 'LS', 'EC', 'MC', 'LC', 'TE', 'PC']

# Get cell to cluster
clusters = nbconfig.seurat.get_clusters('res.0.6')
clusters = clusters[(clusters != 9) & (clusters != 10) & (clusters != 11)].copy()    # drop Unknown clusters
clusters = clusters.map(cluster_annot)
clusters = pd.Series(pd.Categorical(clusters.values, categories=cluster_order, ordered=True), index=pd.Index(clusters.index, name='cell_id'), name='cluster').to_frame()
clusters['rep'] = clusters.index.str.extract('(rep\d)_').values.flatten()
clusters.head()

# %%
# Get fbgn to chromosome mappings 
chroms = nbconfig.fbgn2chrom.query('chrom != "chrM"').copy()    # Drop mitochondrion
autosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R']
chroms.chrom = chroms.chrom.astype('category')
chroms.chrom = chroms.chrom.cat.reorder_categories(nbconfig.chrom_order)

# %%
chrom_cbns = [
    ('chrX', 'chr2L'), 
    ('chrX', 'chr2R'), 
    ('chrX', 'chr3L'), 
    ('chrX', 'chr3R'), 
    ('chrX', 'chrA'), 
    ('chr4', 'chr2L'), 
    ('chr4', 'chr2R'), 
    ('chr4', 'chr3L'), 
    ('chr4', 'chr3R'), 
    ('chr4', 'chrA'), 
    ('chrY', 'chr2L'), 
    ('chrY', 'chr2R'), 
    ('chrY', 'chr3L'), 
    ('chrY', 'chr3R'), 
    ('chrY', 'chrA'), 
]

# %% [markdown]
# ## Chromosome level counts

# %%
def read_chrom_cnt(rep_num):
    chrom_cnt = pd.read_csv(f'../output/scrnaseq-wf/scrnaseq_samples/testis{rep_num}_force/outs/possorted_genome_bam.bam_counts',  sep='\t', index_col=0)
    chrom_cnt.columns = ['chrom', 'UMI']
    chrom_cnt.index = pd.Index([f'rep{rep_num}_{cell_id}' for cell_id in chrom_cnt.index], name='cell_id')
    chrom_cnt.chrom = [f'chr{chrom}' for chrom in chrom_cnt.chrom]
    chrom_cnt_wide = chrom_cnt.set_index('chrom', append=True).unstack().fillna(0)
    chrom_cnt_wide.columns = chrom_cnt_wide.columns.droplevel(0)
    return chrom_cnt_wide

# %%
# get rep 1 chromosome level counts by cell
cnt1 = read_chrom_cnt(1)
cnt1 = cnt1.reindex(clusters.index).dropna()    # Only keep cells that have cluster calls
grps = cnt1.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
        else:
            chrom2 = dd[c2]
        pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
        results.append((c, c1, c2, pval))

dat = -np.log10(pd.DataFrame(results, columns=['cluster', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'chrom1', 'chrom2']) + .0001)

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
dat.plot.bar(ax=ax, legend=False)
ax.set_xlabel('Cluster Chrom1-Chrom2')
ax.set_ylabel('-log10(p-value)')
ax.set_ylim(0, None)
ax.axhline(-np.log10(0.05), color='r', ls=':', label='0.05')
ax.axhline(-np.log10(0.01), color='r', ls='-.', label='0.01')
plt.legend()
ax.set_title('Reads Mapping to Chromosome Arm')

new_labels = []
for l in ax.get_xticklabels():
    txt = l.get_text()
    clus, c1, c2 = re.match(f"\((\w\w), chr([\w\d]+), chr([\w\d]+)\)", txt).groups()
    new_labels.append(f'{clus} {c1:<2}-{c2:<2}')
ax.set_xticklabels(new_labels, fontsize=8, fontdict=dict(family='Monospace'))

loc = 4.5
for i in range(26):
    ax.axvline(loc, color='k', alpha=0.3)
    loc += 5
    
ax.add_patch(plt.Rectangle([-1, -1], width=15.5, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][0]))
ax.add_patch(plt.Rectangle([14.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][1]))
ax.add_patch(plt.Rectangle([29.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][2]))
ax.add_patch(plt.Rectangle([44.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][3]))
ax.add_patch(plt.Rectangle([59.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][4]))
ax.add_patch(plt.Rectangle([74.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][5]))
ax.add_patch(plt.Rectangle([89.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][6]))
ax.add_patch(plt.Rectangle([104.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][7]))
ax.add_patch(plt.Rectangle([119.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][8]))

# %% [markdown]
# ## Normalized Chromosome level counts

# %%
chrom_sizes = pd.read_csv('/data/LCDB/lcdb-references/dmel/r6-16/fasta/dmel_r6-16.chromsizes', sep='\t', index_col=0, header=None)
chrom_sizes.index.name = 'chrom'
chrom_sizes.columns = ['chrom_size']
chrom_sizes = chrom_sizes.reindex(nbconfig.chrom_order)
chrom_sizes = chrom_sizes.chrom_size

# %%
def read_chrom_cnt(rep_num):
    chrom_cnt = pd.read_csv(f'../output/scrnaseq-wf/scrnaseq_samples/testis{rep_num}_force/outs/possorted_genome_bam.bam_counts',  sep='\t', index_col=0)
    chrom_cnt.columns = ['chrom', 'UMI']
    chrom_cnt.index = pd.Index([f'rep{rep_num}_{cell_id}' for cell_id in chrom_cnt.index], name='cell_id')
    chrom_cnt.chrom = [f'chr{chrom}' for chrom in chrom_cnt.chrom]
    chrom_cnt_wide = chrom_cnt.set_index('chrom', append=True).unstack().fillna(0)
    chrom_cnt_wide.columns = chrom_cnt_wide.columns.droplevel(0)
    return chrom_cnt_wide

# %%
# get rep 1 chromosome level counts by cell
cnt1 = read_chrom_cnt(1)
cnt1 = cnt1.reindex(clusters.index).dropna()    # Only keep cells that have cluster calls
cnt1_chrom_length_norm = cnt1.div(chrom_sizes / 1e7) 
grps = cnt1_chrom_length_norm.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
        else:
            chrom2 = dd[c2]
        pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
        results.append((c, c1, c2, pval))

dat = -np.log10(pd.DataFrame(results, columns=['cluster', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'chrom1', 'chrom2']) + .0001)

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
dat.plot.bar(ax=ax, legend=False)
ax.set_xlabel('Cluster Chrom1-Chrom2')
ax.set_ylabel('-log10(p-value)')
ax.set_ylim(0, None)
ax.axhline(-np.log10(0.05), color='r', ls=':', label='0.05')
ax.axhline(-np.log10(0.01), color='r', ls='-.', label='0.01')
plt.legend()
ax.set_title('Reads Mapping to Chromosome Arm / (chromosome length / 1e7)')

new_labels = []
for l in ax.get_xticklabels():
    txt = l.get_text()
    clus, c1, c2 = re.match(f"\((\w\w), chr([\w\d]+), chr([\w\d]+)\)", txt).groups()
    new_labels.append(f'{clus} {c1:<2}-{c2:<2}')
ax.set_xticklabels(new_labels, fontsize=8, fontdict=dict(family='Monospace'))

loc = 4.5
for i in range(26):
    ax.axvline(loc, color='k', alpha=0.3)
    loc += 5
    
ax.add_patch(plt.Rectangle([-1, -1], width=15.5, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][0]))
ax.add_patch(plt.Rectangle([14.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][1]))
ax.add_patch(plt.Rectangle([29.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][2]))
ax.add_patch(plt.Rectangle([44.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][3]))
ax.add_patch(plt.Rectangle([59.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][4]))
ax.add_patch(plt.Rectangle([74.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][5]))
ax.add_patch(plt.Rectangle([89.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][6]))
ax.add_patch(plt.Rectangle([104.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][7]))
ax.add_patch(plt.Rectangle([119.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][8]))

# %% [markdown]
# ## All Gene Counts

# %%
raw = pd.read_parquet('../output/scrnaseq-wf/raw.parquet')

# %%
# get rep 1 chromosome level counts by cell
raw_cnts_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').sum().T
num_genes_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').size().T
norm_raw_cnts_by_chrom = raw_cnts_by_chrom.div(num_genes_by_chrom / 1e3)

# %%
grps = norm_raw_cnts_by_chrom.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
        else:
            chrom2 = dd[c2]
        pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
        results.append((c, c1, c2, pval))

dat = -np.log10(pd.DataFrame(results, columns=['cluster', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'chrom1', 'chrom2']) + .0001)

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
dat.plot.bar(ax=ax, legend=False)
ax.set_xlabel('Cluster Chrom1-Chrom2')
ax.set_ylabel('-log10(p-value)')
ax.set_ylim(0, None)
ax.axhline(-np.log10(0.05), color='r', ls=':', label='0.05')
ax.axhline(-np.log10(0.01), color='r', ls='-.', label='0.01')
plt.legend()
ax.set_title('Reads Mapping to Genes / (number of genes / 1e3)')

new_labels = []
for l in ax.get_xticklabels():
    txt = l.get_text()
    clus, c1, c2 = re.match(f"\((\w\w), chr([\w\d]+), chr([\w\d]+)\)", txt).groups()
    new_labels.append(f'{clus} {c1:<2}-{c2:<2}')
ax.set_xticklabels(new_labels, fontsize=8, fontdict=dict(family='Monospace'))

loc = 4.5
for i in range(26):
    ax.axvline(loc, color='k', alpha=0.3)
    loc += 5
    
ax.add_patch(plt.Rectangle([-1, -1], width=15.5, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][0]))
ax.add_patch(plt.Rectangle([14.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][1]))
ax.add_patch(plt.Rectangle([29.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][2]))
ax.add_patch(plt.Rectangle([44.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][3]))
ax.add_patch(plt.Rectangle([59.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][4]))
ax.add_patch(plt.Rectangle([74.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][5]))
ax.add_patch(plt.Rectangle([89.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][6]))
ax.add_patch(plt.Rectangle([104.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][7]))
ax.add_patch(plt.Rectangle([119.5, -1], width=15, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][8]))

# %%



# %%



# %%



# %% [markdown]
# ## Commonly Expressed Gene Counts

# %%
from larval_gonad.x_to_a import commonly_expressed

# %%
raw = pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
expressed = commonly_expressed(data=raw)
raw = raw.reindex(expressed)
nbconfig.fbgn2chrom.reindex(expressed).chrom.value_counts()

# %%
chrom_cbns_no_Y = [
    ('chrX', 'chr2L'), 
    ('chrX', 'chr2R'), 
    ('chrX', 'chr3L'), 
    ('chrX', 'chr3R'), 
    ('chrX', 'chrA'), 
    ('chr4', 'chr2L'), 
    ('chr4', 'chr2R'), 
    ('chr4', 'chr3L'), 
    ('chr4', 'chr3R'), 
    ('chr4', 'chrA'), 
]

# %%
# get rep 1 chromosome level counts by cell
raw_cnts_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').sum().T
num_genes_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').size().T
norm_raw_cnts_by_chrom = raw_cnts_by_chrom.div(num_genes_by_chrom / 1e1)

# %%
grps = norm_raw_cnts_by_chrom.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns_no_Y:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
        else:
            chrom2 = dd[c2]
        pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
        results.append((c, c1, c2, pval))

dat = -np.log10(pd.DataFrame(results, columns=['cluster', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'chrom1', 'chrom2']) + .0001)

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
dat.plot.bar(ax=ax, legend=False)
ax.set_xlabel('Cluster Chrom1-Chrom2')
ax.set_ylabel('-log10(p-value)')
ax.set_ylim(0, None)
ax.axhline(-np.log10(0.05), color='r', ls=':', label='0.05')
ax.axhline(-np.log10(0.01), color='r', ls='-.', label='0.01')
plt.legend()
ax.set_title('Reads Mapping to Commonly Expressed Genes / (number of genes / 10)')

new_labels = []
for l in ax.get_xticklabels():
    txt = l.get_text()
    clus, c1, c2 = re.match(f"\((\w\w), chr([\w\d]+), chr([\w\d]+)\)", txt).groups()
    new_labels.append(f'{clus} {c1:<2}-{c2:<2}')
ax.set_xticklabels(new_labels, fontsize=8, fontdict=dict(family='Monospace'))

loc = 4.5
for i in range(26):
    ax.axvline(loc, color='k', alpha=0.3)
    loc += 5
    
ax.add_patch(plt.Rectangle([-1, -1], width=10.5, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][0]))
ax.add_patch(plt.Rectangle([9.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][1]))
ax.add_patch(plt.Rectangle([19.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][2]))
ax.add_patch(plt.Rectangle([29.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][3]))
ax.add_patch(plt.Rectangle([39.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][4]))
ax.add_patch(plt.Rectangle([49.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][5]))
ax.add_patch(plt.Rectangle([59.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][6]))
ax.add_patch(plt.Rectangle([69.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][7]))
ax.add_patch(plt.Rectangle([79.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][8]))

# %%



# %%



# %%



# %%



# %% [markdown]
# ## Tau Gene Counts

# %%
tau = pd.read_csv('../output/notebook/2018-02-05_tau_modENCODE_tau.tsv', sep='\t', index_col=0, header=None)
tau.index.name = 'FBgn'
tau.columns = ['tau']
tau = tau.tau.dropna()

# %%
tau.plot(kind='kde')

# %%
housekeeping = tau[tau <= 0.85].index.tolist()

# %%
raw = pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
raw = raw.reindex(housekeeping).dropna()

# %%
raw.join(nbconfig.fbgn2chrom).chrom.value_counts()

# %%
chrom_cbns_no_Y = [
    ('chrX', 'chr2L'), 
    ('chrX', 'chr2R'), 
    ('chrX', 'chr3L'), 
    ('chrX', 'chr3R'), 
    ('chrX', 'chrA'), 
    ('chr4', 'chr2L'), 
    ('chr4', 'chr2R'), 
    ('chr4', 'chr3L'), 
    ('chr4', 'chr3R'), 
    ('chr4', 'chrA'), 
]

# %%
# get rep 1 chromosome level counts by cell
raw_cnts_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').sum().T
num_genes_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').size().T
norm_raw_cnts_by_chrom = raw_cnts_by_chrom.div(num_genes_by_chrom / 1e3)

# %%
grps = norm_raw_cnts_by_chrom.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns_no_Y:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
        else:
            chrom2 = dd[c2]
        pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
        results.append((c, c1, c2, pval))

dat = -np.log10(pd.DataFrame(results, columns=['cluster', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'chrom1', 'chrom2']) + .0001)

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
dat.plot.bar(ax=ax, legend=False)
ax.set_xlabel('Cluster Chrom1-Chrom2')
ax.set_ylabel('-log10(p-value)')
ax.set_ylim(0, None)
ax.axhline(-np.log10(0.05), color='r', ls=':', label='0.05')
ax.axhline(-np.log10(0.01), color='r', ls='-.', label='0.01')
plt.legend()
ax.set_title('Reads Mapping to "tau <=-0.85" Genes / (number of genes / 1e3)')

new_labels = []
for l in ax.get_xticklabels():
    txt = l.get_text()
    clus, c1, c2 = re.match(f"\((\w\w), chr([\w\d]+), chr([\w\d]+)\)", txt).groups()
    new_labels.append(f'{clus} {c1:<2}-{c2:<2}')
ax.set_xticklabels(new_labels, fontsize=8, fontdict=dict(family='Monospace'))

loc = 4.5
for i in range(26):
    ax.axvline(loc, color='k', alpha=0.3)
    loc += 5
    
ax.add_patch(plt.Rectangle([-1, -1], width=10.5, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][0]))
ax.add_patch(plt.Rectangle([9.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][1]))
ax.add_patch(plt.Rectangle([19.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][2]))
ax.add_patch(plt.Rectangle([29.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][3]))
ax.add_patch(plt.Rectangle([39.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][4]))
ax.add_patch(plt.Rectangle([49.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][5]))
ax.add_patch(plt.Rectangle([59.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][6]))
ax.add_patch(plt.Rectangle([69.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][7]))
ax.add_patch(plt.Rectangle([79.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][8]))

# %%



# %%



# %% [markdown]
# ## TSPS Gene Counts

# %%
tsps = pd.read_csv('../output/notebook/2018-02-05_tau_modENCODE_tsps.tsv', sep='\t', index_col=0, header=None)
tsps.index.name = 'FBgn'
tsps.columns = ['tsps']
tsps = tsps.tsps.dropna()

# %%
tsps.plot(kind='kde')

# %%
housekeeping = tsps[tsps <= 1].index.tolist()

# %%
raw = pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
raw = raw.reindex(housekeeping).dropna()

# %%
raw.join(nbconfig.fbgn2chrom).chrom.value_counts()

# %%
chrom_cbns_no_Y = [
    ('chrX', 'chr2L'), 
    ('chrX', 'chr2R'), 
    ('chrX', 'chr3L'), 
    ('chrX', 'chr3R'), 
    ('chrX', 'chrA'), 
    ('chr4', 'chr2L'), 
    ('chr4', 'chr2R'), 
    ('chr4', 'chr3L'), 
    ('chr4', 'chr3R'), 
    ('chr4', 'chrA'), 
]

# %%
# get rep 1 chromosome level counts by cell
raw_cnts_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').sum().T
num_genes_by_chrom = raw.join(nbconfig.fbgn2chrom).groupby('chrom').size().T
norm_raw_cnts_by_chrom = raw_cnts_by_chrom.div(num_genes_by_chrom / 1e3)

# %%
grps = norm_raw_cnts_by_chrom.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns_no_Y:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
        else:
            chrom2 = dd[c2]
        pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
        results.append((c, c1, c2, pval))

dat = -np.log10(pd.DataFrame(results, columns=['cluster', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'chrom1', 'chrom2']) + .0001)

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 8))
dat.plot.bar(ax=ax, legend=False)
ax.set_xlabel('Cluster Chrom1-Chrom2')
ax.set_ylabel('-log10(p-value)')
ax.set_ylim(0, None)
ax.axhline(-np.log10(0.05), color='r', ls=':', label='0.05')
ax.axhline(-np.log10(0.01), color='r', ls='-.', label='0.01')
plt.legend()
ax.set_title('Reads Mapping to "tau <=-0.85" Genes / (number of genes / 1e3)')

new_labels = []
for l in ax.get_xticklabels():
    txt = l.get_text()
    clus, c1, c2 = re.match(f"\((\w\w), chr([\w\d]+), chr([\w\d]+)\)", txt).groups()
    new_labels.append(f'{clus} {c1:<2}-{c2:<2}')
ax.set_xticklabels(new_labels, fontsize=8, fontdict=dict(family='Monospace'))

loc = 4.5
for i in range(26):
    ax.axvline(loc, color='k', alpha=0.3)
    loc += 5
    
ax.add_patch(plt.Rectangle([-1, -1], width=10.5, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][0]))
ax.add_patch(plt.Rectangle([9.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][1]))
ax.add_patch(plt.Rectangle([19.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][2]))
ax.add_patch(plt.Rectangle([29.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][3]))
ax.add_patch(plt.Rectangle([39.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][4]))
ax.add_patch(plt.Rectangle([49.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][5]))
ax.add_patch(plt.Rectangle([59.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][6]))
ax.add_patch(plt.Rectangle([69.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][7]))
ax.add_patch(plt.Rectangle([79.5, -1], width=10, height=10, zorder=0, alpha=.4, color=nbconfig.colors['clusters'][8]))

# %%



# %% [markdown]
# ## Example of how to display

# %%
def bootstrap(dat, n_boot=1000, estimator=np.median):
    results = np.empty(n_boot)
    for i in range(n_boot):
        results[i] = estimator(dat.sample(n=dat.shape[0], replace=True))
    return np.percentile(results, [2.5, 97.5])


# %%
chrom_sizes = pd.read_csv('/data/LCDB/lcdb-references/dmel/r6-16/fasta/dmel_r6-16.chromsizes', sep='\t', index_col=0, header=None)
chrom_sizes.index.name = 'chrom'
chrom_sizes.columns = ['chrom_size']
chrom_sizes = chrom_sizes.reindex(nbconfig.chrom_order)
chrom_sizes = chrom_sizes.chrom_size

# %%
def read_chrom_cnt(rep_num):
    chrom_cnt = pd.read_csv(f'../output/scrnaseq-wf/scrnaseq_samples/testis{rep_num}_force/outs/possorted_genome_bam.bam_counts',  sep='\t', index_col=0)
    chrom_cnt.columns = ['chrom', 'UMI']
    chrom_cnt.index = pd.Index([f'rep{rep_num}_{cell_id}' for cell_id in chrom_cnt.index], name='cell_id')
    chrom_cnt.chrom = [f'chr{chrom}' for chrom in chrom_cnt.chrom]
    chrom_cnt_wide = chrom_cnt.set_index('chrom', append=True).unstack().fillna(0)
    chrom_cnt_wide.columns = chrom_cnt_wide.columns.droplevel(0)
    return chrom_cnt_wide

# %%
# get rep 1 chromosome level counts by cell
cnt1 = read_chrom_cnt(1)
cnt1 = cnt1.reindex(clusters.index).dropna()    # Only keep cells that have cluster calls
cnt1_chrom_length_norm = cnt1.div(chrom_sizes / 1e7) 

# %%
chrom_cbns = [
    ('chrX', 'chrA'), 
    ('chr4', 'chrA'), 
    ('chrY', 'chrA'), 
]

# %%
grps = cnt1_chrom_length_norm.join(clusters).groupby('cluster')

results = []
for c, dd in grps:
    for c1, c2 in chrom_cbns:
        chrom1 = dd[c1]
        if c2 == 'chrA':
            chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
            ratio = chrom1 / chrom2
            median_ratio = np.median(ratio)
            lower, upper = bootstrap(ratio)
            log2fold = np.log2(median_ratio)
            pval = permutation_test_chrom1_lt_chrom2(chrom1, chrom2)
            results.append((c, c1, median_ratio, upper, lower, log2fold, pval))

dat = pd.DataFrame(results, columns=['cluster', 'chrom', 'median_ratio', 'upper_ci', 'lower_ci', 'log2fold_change', 'p_value']).set_index(['cluster', 'chrom'])

# %%
datx = dat.query('chrom == "chrX"')
datx.index = datx.index.droplevel('chrom')

# %%
dat4 = dat.query('chrom == "chr4"')
dat4.index = dat4.index.droplevel('chrom')

# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 10))
datx.median_ratio.plot(ax=ax, marker='o', color='k', label='chrX')
dat4.median_ratio.plot(ax=ax, marker='^', color='gray', label='chr4')

for i, (l, dd) in enumerate(datx.iterrows()):
    _color = nbconfig.colors['clusters'][i]
    ax.plot([i, i], [dd.lower_ci, dd.upper_ci], color=_color, zorder=0, lw=1.5)
    ax.plot([i-.1, i+.1], [dd.lower_ci, dd.lower_ci], color=_color, zorder=0, lw=1.5)
    ax.plot([i-.1, i+.1], [dd.upper_ci, dd.upper_ci], color=_color, zorder=0, lw=1.5)
    ax.scatter([i, ], [dd.median_ratio, ], color=_color, s=100, zorder=10)
    if dd.p_value <= 0.001:
        ax.text(i, dd.upper_ci, '***', ha='center', va='bottom')
    elif dd.p_value <= 0.01:
        ax.text(i, dd.upper_ci, '**', ha='center', va='bottom')
    elif dd.p_value <= 0.05:
        ax.text(i, dd.upper_ci, '*', ha='center', va='bottom')
    
for i, (l, dd) in enumerate(dat4.iterrows()):
    _color = nbconfig.colors['clusters'][i]
    ax.plot([i, i], [dd.lower_ci, dd.upper_ci], color=_color, zorder=0, lw=1.5)
    ax.plot([i-.1, i+.1], [dd.lower_ci, dd.lower_ci], color=_color, zorder=0, lw=1.5)
    ax.plot([i-.1, i+.1], [dd.upper_ci, dd.upper_ci], color=_color, zorder=0, lw=1.5)
    ax.scatter([i, ], [dd.median_ratio, ], marker='^', color=_color, s=100, zorder=10)
    if dd.p_value <= 0.001:
        ax.text(i, dd.upper_ci, '***', ha='center', va='bottom')
    elif dd.p_value <= 0.01:
        ax.text(i, dd.upper_ci, '**', ha='center', va='bottom')
    elif dd.p_value <= 0.05:
        ax.text(i, dd.upper_ci, '*', ha='center', va='bottom')
    
ax.set_xticks(range(9))
ax.set_xticklabels(datx.index.tolist())
ax.axhline(1, color='r', alpha=.3, ls='--')
ax.set_ylim(0.25, 1.45)
ax.set_ylabel('X:A Ratio')

plt.legend()
plt.title('X and 4 to Autosome Ratio\n(Mapping to Chromosome)')

# %%



# %%



# %%
fig, ax = plt.subplots(1, 1, figsize=(20, 10))
(-1 * datx.log2fold_change).plot(ax=ax, marker='o', color='k', label='chrX')
(-1 * dat4.log2fold_change).plot(ax=ax, marker='^', color='gray', label='chr4')

ax.set_xticks(range(9))
ax.set_xticklabels(datx.index.tolist())
ax.axhline(0, color='r', alpha=.3, ls='--')
ax.set_ylabel('-log2FoldChange')

plt.legend()

# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



# %%



