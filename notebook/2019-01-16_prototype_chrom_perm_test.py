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
# Replicate 1 seems to work well with this type of permutation test, however replicates 2 and 3 end up calling everything significant.
#
# As a reminder the permutation test that I am doing is as follows.
#
# For each cell I am calculating the median ratio of (X / Autosome). I then permute X and Autosome labels generating random cells and compare the median ratio. I use this to generate a null distribution and calculate an empirical p-value based on the number of permuted median ratios that are more extreme than the observed ratio. 
#
# After talking with Cameron, he thinks that this algorithm is a little strange. Instead he thinks I was on the right track earlier when I was permuting cell type labels and keeping the cells whole. I will explore his suggested algorithm later. 
#
# Here I show that while rep 1 behaves as hypothesized, reps 2 and 3 show everything is significant. 

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
# ## Normalized Chromosome level counts

# %%
fname = '/data/LCDB/lcdb-references/dmel/r6-16/fasta/dmel_r6-16.chromsizes'
chrom_lengths = pd.read_csv(fname, sep='\t', index_col=0, header=None, names=['chrom', 'chrom_size'])
chrom_lengths = chrom_lengths.reindex(nbconfig.chrom_order)
chrom_lengths = chrom_lengths.chrom_size

# %%
def read_chrom_cnt(rep_num, chrom_lengths):
    fname = f'../output/scrnaseq-wf/scrnaseq_samples/testis{rep_num}_force/outs/possorted_genome_bam.bam_counts'
    chrom_cnt = pd.read_csv(fname, sep='\t', index_col=0, header=0, names=['cell_id', 'chrom', 'UMI'])
    # Add `rep#_` and `chr` prefixes
    chrom_cnt.index = pd.Index([f'rep{rep_num}_{cell_id}' for cell_id in chrom_cnt.index], name='cell_id')
    chrom_cnt.chrom = [f'chr{chrom}' for chrom in chrom_cnt.chrom]
    
    chrom_cnt_wide = chrom_cnt.set_index('chrom', append=True).unstack().fillna(0)
    chrom_cnt_wide.columns = chrom_cnt_wide.columns.droplevel(0)
    
    num_reads_per_cell = chrom_cnt_wide.sum(axis=1)
    #chrom_cnt_wide_norm = chrom_cnt_wide.div(num_reads_per_cell / 1e3, axis='index').div(chrom_lengths / 1e7)
    chrom_cnt_wide_norm = (
        chrom_cnt_wide
            .div(num_reads_per_cell / 1e3, axis='index')
            .div(chrom_lengths / 1e7)
    )
    
    return chrom_cnt_wide_norm

# %% [markdown]
# ### Rep 1

# %%
# get rep 1 chromosome level counts by cell
cnt1 = read_chrom_cnt(1, chrom_lengths)
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
ax.set_title('Rep 1 Reads Mapping to Chromosome Arm / (chromosoMe length / 1e7)')

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


# %% [markdown]
# ### Rep 2

# %%
cnt = read_chrom_cnt(2, chrom_lengths)
grp = cnt.join(clusters).groupby('cluster')
dat = grp.get_group('EC')

# %%
ax = np.log10(dat[autosomes + ['chrX']]).plot(kind='kde')
ax.axvline(np.log10(dat['chr2L'].median()), color='blue', ls='--', alpha=.5)
ax.axvline(np.log10(dat['chr2R'].median()), color='orange', ls='--', alpha=.5)
ax.axvline(np.log10(dat['chr3L'].median()), color='green', ls='--', alpha=.5)
ax.axvline(np.log10(dat['chr3R'].median()), color='red', ls='--', alpha=.5)
ax.axvline(np.log10(dat['chrX'].median()), color='purple', ls=':', alpha=.5)
ax.axvline(np.log10(dat[autosomes].median(axis=1).median()), color='k', ls='-.', alpha=.5)

# %%
# get rep 2 chromosome level counts by cell
cnt2 = read_chrom_cnt(2, chrom_lengths)
cnt2 = cnt2.reindex(clusters.index).dropna()    # Only keep cells that have cluster calls
grps = cnt2.join(clusters).groupby('cluster')

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
ax.set_title('Rep 2 Reads Mapping to Chromosome Arm / (chromosome length / 1e7)')

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
# ### Rep 3

# %%
# get rep 3 chromosome level counts by cell
cnt3 = read_chrom_cnt(3)
cnt3 = cnt3.reindex(clusters.index).dropna()    # Only keep cells that have cluster calls
cnt3_chrom_length_norm = cnt3.div(chrom_sizes / 1e7)
grps = cnt3_chrom_length_norm.join(clusters).groupby('cluster')

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
ax.set_title('Rep 3 Reads Mapping to Chromosome Arm / (chromosome length / 1e7)')

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

