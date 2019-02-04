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

# %%
def read_fbgn2chrom():
    mapper = {
        'chrX': 'X',
        'chrY': 'Y',
        'chr4': '4',
        'chr2L': 'A',
        'chr2R': 'A',
        'chr3L': 'A',
        'chr3R': 'A',
    }

    fbgn2chrom = (pd.read_csv('../output/fbgn2chrom.tsv', sep='\t', index_col=0)
                      .query('chrom != "chrM"')
                      .chrom.map(mapper)
                      .astype('category')
                      .cat.as_ordered()
                 )
    
    return fbgn2chrom.cat.reorder_categories(['X', 'A', 'Y', '4'])


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

# %% [markdown]
# ### Permutation Test

# %% [markdown]
# At the experiment level, it is clear that missingness is not random. This maybe due to technical artifacts such as dropout, or maybe related to biological processes (i.e. RNA-content of somatic cells is much smaller than germline). Therefore even a non-parametric test is not appropriate, unless we model the missingness (which is very hard). 
#
# Fortunately, at the cell level missingness appears to be somewhat random in relation to X and A expression. We have proposed using a permutation approach 

# %%
cell_ids = []
flags = []
for cell_id, dd in df.groupby('cell_id'):
    y_data = dd[dd.chrom == "Y"].UMI.values
    a_data = dd[dd.chrom == "A"].UMI.values
    _, p_value = mannwhitneyu(y_data, a_data, alternative='greater')
    
    if p_value <= 0.05:
        flags.append(True)
    else:
        flags.append(False)
        
    cell_ids.append(cell_id)

flag_y_gt_a = pd.Series(flags, index=pd.Index(cell_ids, name='cell_id'), name='flag_y_gt_a')

# %%
flag_y_gt_a_by_cluster = pd.concat([flag_y_gt_a, read_clusters()], axis=1, sort=True)
flag_y_gt_a_by_cluster['rep'] = flag_y_gt_a_by_cluster.index.str.extract('(?P<rep>rep\d)', expand=False)

# %%
obs = flag_y_gt_a_by_cluster.groupby('cluster').mean()

# %%
results = []
for i in range(10_000):
    _df = flag_y_gt_a_by_cluster.copy()
    _df.cluster = _df.cluster.sample(frac=1, replace=False).values
    perm = _df.groupby('cluster').mean()
    perm.columns = [f'permutation_{i}']
    results.append(perm)

perms = pd.concat(results, axis=1).T

# %%
ax = perms.plot(kind='kde')
ax.set_xlabel('Proportion of Cells (Y > A)')
ax.set_title('Permuted Null Distributions by Cluster')

# %%
# Calculate p-value using permutation
results = []
for clus, dd in obs.groupby('cluster'):
    _obs = dd.iloc[0, 0]
    p_val = (perms[clus] >= _obs).mean()
    results.append((clus, p_val))

perm_pval = pd.DataFrame(results, columns=['cluster', 'p_value']).set_index('cluster').p_value

# %%
# calculate bootstrap confidence intervals for plotting
def bootstrap(dat, n_boot=1000, estimator=np.mean):
    results = np.empty(n_boot)
    for i in range(n_boot):
        results[i] = estimator(dat.sample(n=dat.shape[0], replace=True))
    return np.percentile(results, [2.5, 97.5])

prop_flag_y_gt_a = flag_y_gt_a_by_cluster.groupby(['cluster', 'rep']).flag_y_gt_a.mean().to_frame().reset_index()

results = []
for clus, dd in prop_flag_y_gt_a.groupby('cluster'):
    low, high = bootstrap(dd.flag_y_gt_a)
    results.append((clus, low, high))
    
cluster_bootstrap = pd.DataFrame(results, columns=['cluster', 'low', 'high'])
cluster_bootstrap_w_pval = cluster_bootstrap.join(perm_pval, on='cluster')

# %%
fig, ax = plt.subplots(figsize=plt.figaspect(1/2))
ax.plot(prop_flag_y_gt_a.groupby("cluster").flag_y_gt_a.mean().values, color='k', zorder=-10)
sns.pointplot(x='cluster', y='flag_y_gt_a', data=prop_flag_y_gt_a, errwidth=2, capsize=.2, palette=nbconfig.colors['clusters'], zorder=10, ax=ax)
#ax.set_ylim(0, 1.1)
ax.set_ylabel('Prop Cells')
ax.set_title('Proprotion of Cells with Y Enrichment')

for i, row in cluster_bootstrap_w_pval.iterrows():
    if row.p_value <= 0.05:
        ax.text(i, row.high, '*', ha='center', va='bottom')

# %%
prop_missing_by_cell_by_chrom = df.groupby(['cluster', 'cell_id', 'chrom']).missing.mean().unstack()
prop_missing_by_cell_by_chrom.columns = [f'{x}_missingness' for x in prop_missing_by_cell_by_chrom.columns]
prop_missing_by_cell_by_chrom.reset_index(inplace=True)

# %%
dat = prop_missing_by_cell_by_chrom.join(flag_x_lt_a.astype(int), on='cell_id')

# %%
g = sns.FacetGrid(dat, col='cluster', col_wrap=4)
g.map(sns.kdeplot, 'X_missingness', label='X Missingness')
g.map(sns.kdeplot, 'A_missingness', color='orange', label='A Missingness')
g.axes[0].legend()

# %%
results = smf.logit('flag_x_lt_a ~ (X_missingness + A_missingness)/cluster', data=dat).fit()
results.summary2()
