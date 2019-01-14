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
# # What permutation tests to do?

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
from collections import defaultdict
from itertools import combinations

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import mannwhitneyu, ks_2samp, friedmanchisquare, kruskal
from statsmodels.stats.multitest import multipletests

# Project level imports
from larval_gonad.notebook import Nb
from larval_gonad.x_to_a import commonly_expressed
from larval_gonad.stats import permutation_sample, permuted_replicates, enrichment_statistic

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
print(cnt1.shape)
cnt1.head()

# %%
grps = cnt1.join(clusters).groupby('cluster')

# %% [markdown]
# ## Does X differ from 2L

# %%
def run(g1):
    x = g1.chrX / 23542271 / 1e7
    a = g1.chr2L / 23513712 / 1e7
    obs = x / a

    cnt = 0
    lt = 0
    for i in range(10_000):
        px, pa = permutation_sample(x, a)
        px = px / 23542271 / 1e7
        pa = pa  / 23513712 / 1e7
        pr = px / pa
        _, pval = mannwhitneyu(obs, pr, alternative='less')
        if pval <= 0.01:
            lt += 1
        cnt += 1

    return 1 - (lt / cnt)

# %%
def run(g1, size=100):
    x = g1.chrX / 23542271 / 1e7
    a = g1.chr2L / 23513712 / 1e7
    obs = np.median(x / a)

    perm_results = np.empty(size)
    for i in range(size):
        px, pa = permutation_sample(x, a)
        px = px / 23542271 / 1e7
        pa = pa  / 23513712 / 1e7
        perm_results[i] = np.median(px / pa)

    return sum(perm_results <= obs) / len(perm_results)

# %%
def run(g1, size=100):
    x = g1.chrX
    a = g1.chr2L
    obs = np.median(x / a)

    perm_results = np.empty(size)
    for i in range(size):
        px, pa = permutation_sample(x, a)
        px = px
        pa = pa
        perm_results[i] = np.median(px / pa)

    return sum(perm_results <= obs) / len(perm_results)

# %%
from larval_gonad.stats import permutation_test_chrom1_lt_chrom2

# %%


# %%
g1 = grps.get_group('SP')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('ES')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('MS')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('LS')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('EC')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('MC')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('LC')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('PC')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%
g1 = grps.get_group('TE')
permutation_test_chrom1_lt_chrom2(g1.chrX, g1.chr2L, alternative='less')

# %%



# %%



# %% [markdown]
# ## Does X:2L differ among SP and spermatocytes

# %%
from scipy.spatial.distance import cosine

# %%
r1 = cnt1.chrX / cnt1.chr2L
r1.name = 'ratio'

# %%
g1_name = 'SP'
g2_name = 'ES'
grps = clusters.join(r1, how='right').groupby('cluster')
g1 = grps.get_group(g1_name).ratio
g2 = grps.get_group(g2_name).ratio

sns.kdeplot(g1, cumulative=False, label=g1_name)
sns.kdeplot(g2, cumulative=False, label=g2_name)

mannwhitneyu(g1, g2)

# %%
def run(g1, g2, size=100):
    _, opval = mannwhitneyu(g1, g2, alternative='two-sided')
    perm_results = np.empty(size)
    for i in range(size):
        pg1, pg2 = permutation_sample(g1, g2)
        _, ppval = mannwhitneyu(pg1, pg2, alternative='two-sided')
        perm_results[i] = pval
    return sum(perm_results <= opval) / len(perm_results)

# %%
def run(g1, g2, size=1_000, plot=False):
    _g1 = g1.median()
    _g2 = g2.median()
    odiff = np.abs(_g1 - _g2)
    perm_results = np.empty(size)
    
    for i in range(size):
        pg1, pg2 = permutation_sample(g1, g2)
        _pg1 = np.median(pg1)
        _pg2 = np.median(pg2)
        pdiff = np.abs(_pg1 - _pg2)
        perm_results[i] = pdiff
        
    if plot:
        sns.kdeplot(g1, color='blue')
        sns.kdeplot(g2, color='orange')
        sns.kdeplot(pg1, color='black', alpha=.2, zorder=0)
        sns.kdeplot(pg2, color='yellow', alpha=.2, zorder=0)
        
        ax = plt.gca()
        ax.axvline(_g1, color='blue', ls='--', alpha=.4)
        ax.axvline(_g2, color='orange', ls='--', alpha=.4)
        ax.axvline(_pg1, color='black', ls='--', alpha=.4)
        ax.axvline(_pg2, color='yellow', ls='--', alpha=.4)
        
    return sum(perm_results >= odiff) / len(perm_results)

# %%
def run(g1, g2, size=100):
    _g1 = g1.sample(100)
    _g2 = g2.sample(100)
    odiff = cosine(_g1, _g2)
    perm_results = np.empty(size)
    for i in range(size):
        pg1, pg2 = permutation_sample(_g1, _g2)
        pdiff = cosine(pg1, pg2)
        perm_results[i] = pdiff
    return sum(perm_results >= odiff) / len(perm_results)

# %%
sp = grps.get_group('SP').ratio
es = grps.get_group('ES').ratio
ms = grps.get_group('MS').ratio
ls = grps.get_group('LS').ratio
ec = grps.get_group('EC').ratio
mc = grps.get_group('MC').ratio
lc = grps.get_group('LC').ratio
pc = grps.get_group('PC').ratio
te = grps.get_group('TE').ratio

# %% [markdown]
# ### SP

# %%
run(sp, es)

# %%
run(sp, ms)

# %%
run(sp, ls)

# %%
run(sp, ec)

# %%
run(sp, mc)

# %%
run(sp, lc)

# %%
run(sp, pc)

# %%
run(sp, te)

# %%



# %% [markdown]
# ### ES

# %%
run(es, sp)

# %%
run(es, ms)

# %%
run(es, ls)

# %%
run(es, ec)

# %%
run(es, mc)

# %%
run(es, lc)

# %%
run(es, pc)

# %%
run(es, te)

# %%



# %% [markdown]
# ### MS

# %%
run(ms, sp)

# %%
run(ms, es)

# %%
run(ms, ls)

# %%
run(ms, ec)

# %%
run(ms, mc)

# %%
run(ms, lc)

# %%
run(ms, pc)

# %%
run(ms, te)

# %%



# %%



# %% [markdown]
# ### LS

# %%
run(ls, sp)

# %%
run(ls, es)

# %%
run(ls, ms)

# %%
run(ls, ec)

# %%
run(ls, mc)

# %%
run(ls, lc)

# %%
run(ls, pc)

# %%
run(ls, te)

# %%



# %%



# %% [markdown]
# ### EC

# %%
run(ec, sp)

# %%
run(ec, es)

# %%
run(ec, ms)

# %%
run(ec, ls)

# %%
run(ec, mc)

# %%
run(ec, lc)

# %%
run(ec, pc)

# %%
run(ec, te)

# %%



# %%



# %% [markdown]
# ### MC

# %%
run(mc, sp)

# %%
run(mc, es)

# %%
run(mc, ms)

# %%
run(mc, ls)

# %%
run(mc, ec)

# %%
run(mc, lc)

# %%
run(mc, pc)

# %%
run(mc, te)

# %%



# %%



# %% [markdown]
# ### LC

# %%
run(lc, sp)

# %%
run(lc, es)

# %%
run(lc, ms)

# %%
run(lc, ls)

# %%
run(lc, ec)

# %%
run(lc, mc)

# %%
run(lc, pc)

# %%
run(lc, te)

# %%



# %%



# %% [markdown]
# ### PC

# %%
run(pc, sp)

# %%
run(pc, es)

# %%
run(pc, ms)

# %%
run(pc, ls)

# %%
run(pc, ec)

# %%
run(pc, mc)

# %%
run(pc, lc)

# %%
run(pc, te)

# %%



# %%



# %% [markdown]
# ### TE

# %%
run(te, sp)

# %%
run(te, es)

# %%
run(te, ms)

# %%
run(te, ls)

# %%
run(te, ec)

# %%
run(te, mc)

# %%
run(te, lc)

# %%
run(pc, te)

# %%



# %%



# %%


