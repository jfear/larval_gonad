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
# # Y Permutation Test

# %%
import os
import sys
import re
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact, cumfreq
from scipy.cluster.hierarchy import linkage, dendrogram
import statsmodels.formula.api as smf

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb
from larval_gonad.x_to_a import commonly_expressed

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %% [markdown]
# ## Data Prep

# %%
clusters = (
    nbconfig.seurat.get_clusters('res.0.6')
    .map(nbconfig.short_cluster_annot)
    .pipe(lambda x: x[x != 'UNK'])
)

clusters.head()

# %%
mapper = {'chrX': 'X', 'chrY': 'Y', 'chr4': '4', 'chr2L': 'A','chr2R': 'A','chr3L': 'A','chr3R': 'A',}
fbgn2chrom = (
    pd.read_csv('../output/fbgn2chrom.tsv', sep='\t', index_col=0)
    .chrom.map(mapper)
    .dropna()
)
fbgn2chrom.head()

# %%
num_genes_per_chrom = fbgn2chrom.value_counts()
num_genes_per_chrom

# %%
stacked = (
    nbconfig.seurat.get_raw()
    .reset_index()
    .join(fbgn2chrom, on='FBgn')
    .melt(id_vars=['FBgn', 'chrom'], var_name='cell_id', value_name='UMI')
    .join(clusters, on='cell_id')
)

stacked.head()

# %%
ratios_by_cell = (
    stacked.groupby(['cell_id', 'cluster', 'chrom'])
    .UMI.sum()
    .div(num_genes_per_chrom, level='chrom')
    .mul(1e3)
    .unstack()
    .assign(ratio_x = lambda df: df.X / df.A)
    .assign(ratio_y = lambda df: df.Y / df.A)
    .assign(ratio_4 = lambda df: df['4'] / df.A)
    .drop(['X', 'Y', 'A', '4'], axis=1)
    .reset_index('cluster')
)

ratios_by_cell.head()

# %%
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=plt.figaspect(1/3))
(
    ratios_by_cell.groupby('cluster')
    .ratio_x.apply(sns.kdeplot, ax=ax1)
)
ax1.set_title('X Ratio')

(
    ratios_by_cell.groupby('cluster')
    .ratio_4.apply(sns.kdeplot, ax=ax2)
)
ax2.set_title('4 Ratio')

(
    ratios_by_cell.groupby('cluster')
    .ratio_y.apply(sns.kdeplot, ax=ax3)
)
ax3.set_title('Y Ratio')

# %%
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=plt.figaspect(1/3), sharey=True, sharex=True)

_defaults = dict(x='cluster', data=ratios_by_cell, order=nbconfig.short_cluster_order, showfliers=False, palette=nbconfig.colors['clusters'])

sns.boxplot(y='ratio_x', ax=ax1, **_defaults)
ax1.set_title('X Ratio')

sns.boxplot(y='ratio_4', ax=ax2, **_defaults)
ax2.set_title('4 Ratio')

sns.boxplot(y='ratio_y', ax=ax3, **_defaults)
ax3.set_title('Y Ratio')

for ax in [ax1, ax2, ax3]:
    ax.axhline(1, ls=':', color='grey')

# %%
median_ratios = ratios_by_cell.groupby('cluster').median().reindex(nbconfig.short_cluster_order)

# %%
pd.concat([
    (median_ratios.loc['SP'] / median_ratios.loc['ES']).rename('SP/ES'), 
    (median_ratios.loc['SP'] / median_ratios.loc['MS']).rename('SP/MS'),
    (median_ratios.loc['SP'] / median_ratios.loc['LS']).rename('SP/LS')
], axis=1, sort=True)

# %%
1 / median_ratios

# %%
median_ratios

# %%


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
fig, axes = plt.subplots(1, 3, figsize=plt.figaspect(1/3), sharey=True)

_defaults = dict(x='cluster', data=ratios_by_cell, order=nbconfig.short_cluster_order, showfliers=False, palette=nbconfig.colors['clusters'])

sns.boxplot(y='ratio_x', ax=axes[0], **_defaults)
plot_pval(ratios_by_cell[['cluster', 'ratio_x']], pvals.pval_x_lt_a, axes[0])
axes[0].text(.01, .95, '* X < A', transform=axes[0].transAxes, fontsize=14)


sns.boxplot(y='ratio_4', ax=axes[1], **_defaults)
plot_pval(ratios_by_cell[['cluster', 'ratio_4']], pvals.pval_4_lt_a, axes[1])
axes[1].text(.01, .95, '* 4 < A', transform=axes[1].transAxes, fontsize=14)


sns.boxplot(y='ratio_y', ax=axes[2], **_defaults)
plot_pval(ratios_by_cell[['cluster', 'ratio_y']], pvals.pval_y_gt_a, axes[2])
axes[2].text(.01, .95, '* A < Y', transform=axes[2].transAxes, fontsize=14)
