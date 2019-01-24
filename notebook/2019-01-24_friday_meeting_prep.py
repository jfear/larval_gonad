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
# # Friday Meeting Prep

# %% [markdown]
# I want to report on chromosome expression and X:A in this weeks Friday meeting. Here is where I am developing those plots.

# %%
import os
import sys
import re
from pathlib import Path

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact
from scipy.cluster.hierarchy import linkage, dendrogram
import statsmodels.formula.api as smf

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb

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
        .cat.reorder_categories(['X', 'A', 'Y', '4'])
)

# %% [markdown] {"toc-hr-collapsed": true}
# ## Chromosome Expression

# %% [markdown]
# ### Cell level chromosome coverage

# %%
g = sns.FacetGrid(norm_cnts, col='chrom', col_wrap=2, sharey=False)
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

# %% [markdown] {"toc-hr-collapsed": false}
# ### Rep level chromosome coverage

# %%
dat = norm_cnts.groupby(['cluster', 'rep', 'chrom']).norm_cnt.median().to_frame().reset_index()

# %%
g = sns.FacetGrid(dat, col='chrom', col_wrap=2, sharey=False)
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

# %%
del dat

# %% [markdown]
# ### Y Gene Expression

# %%
prop_missing_by_cluster_by_gene = df.groupby(['cluster', 'chrom', 'FBgn']).missing.mean().to_frame().reset_index()

# %%
y_genes = (1 - prop_missing_by_cluster_by_gene.query('chrom == "Y"').set_index(['FBgn', 'cluster']).missing).unstack()

# remove the genes with all 0's
y_genes = y_genes[(y_genes >= 0.05).any(axis=1)]

# %%
tree = dendrogram(linkage(y_genes.values, 'average'), no_plot=True)
leaves = tree['leaves']

fig = plt.figure(figsize=(6, 10))
ax = sns.heatmap(y_genes.iloc[leaves, :], yticklabels=True)
ax.set_title('Proportion Cells with Expression')
labels = []
for label in ax.get_yticklabels():
    labels.append(nbconfig.fbgn2symbol[label.get_text()])
ax.set_yticklabels(labels, fontsize=8);

# %%
del prop_missing_by_cluster_by_gene

# %%
fbgn2symbol = pd.Series(nbconfig.fbgn2symbol, name='gene_symbol')
fbgn2symbol.index.name = 'FBgn'
y_cnts = df.query('chrom == "Y"').groupby(['FBgn', 'cluster']).UMI.sum().to_frame().reset_index().join(fbgn2symbol, on='FBgn')

mask = y_cnts.groupby('FBgn').UMI.sum() > 10

keeps = mask[mask].index.tolist()

g = sns.FacetGrid(y_cnts.query(f'FBgn == {keeps}').sort_values('gene_symbol'), col='gene_symbol', col_wrap=8)
g.map(sns.pointplot, 'cluster', 'UMI', order=nbconfig.short_cluster_order)
g.set_titles('{col_name}')
plt.suptitle('Total Gene Expression sum(UMI)', va='bottom', y=.99)

# %%
del y_cnts

# %%
deg = pd.read_csv('../output/scrnaseq-wf/germcell_soma_deg/germ_vs_cysts.tsv', sep='\t', index_col=0).join(fbgn2chrom)

# %%
deg.query('chrom == "Y"')

# %%
deg = pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0).join(fbgn2chrom)

# %%
deg.query('chrom == "Y"')

# %%
biomarkers = nbconfig.seurat.get_biomarkers('res.0.6').join(fbgn2chrom)

# %%
biomarkers.query('chrom == "Y"')

# %% [markdown]
# ### 4th Expression

# %%
_4_genes = (1 - prop_missing_by_cluster_by_gene.query('chrom == "4"').set_index(['FBgn', 'cluster']).missing).unstack()

# remove the genes with all 0's
_4_genes = _4_genes[(_4_genes >= 0.05).any(axis=1)]

# %%
tree = dendrogram(linkage(_4_genes.values, 'average'), no_plot=True)
leaves = tree['leaves']

fig = plt.figure(figsize=(6, 10))
ax = sns.heatmap(_4_genes.iloc[leaves, :], yticklabels=True)
ax.set_title('Proportion Cells with Expression')
labels = []
for label in ax.get_yticklabels():
    labels.append(nbconfig.fbgn2symbol[label.get_text()])
ax.set_yticklabels(labels, fontsize=8);

# %%
del _4_genes

# %%
fbgn2symbol = pd.Series(nbconfig.fbgn2symbol, name='gene_symbol')
fbgn2symbol.index.name = 'FBgn'
_4_cnts = df.query('chrom == "4"').groupby(['FBgn', 'cluster']).UMI.sum().to_frame().reset_index().join(fbgn2symbol, on='FBgn')

mask = _4_cnts.groupby('FBgn').UMI.sum() > 1e3

keeps = mask[mask].index.tolist()

g = sns.FacetGrid(_4_cnts.query(f'FBgn == {keeps}').sort_values('gene_symbol'), col='gene_symbol', col_wrap=8)
g.map(sns.pointplot, 'cluster', 'UMI', order=nbconfig.short_cluster_order)
g.set_titles('{col_name}')
plt.suptitle('Total Gene Expression sum(UMI)', va='bottom', y=1)

# %%
del _4_cnts

# %%
deg = pd.read_csv('../output/scrnaseq-wf/germcell_soma_deg/germ_vs_cysts.tsv', sep='\t', index_col=0).join(fbgn2chrom)

# %%
deg.query('chrom == "4"')

# %%
deg = pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0).join(fbgn2chrom)

# %%
deg.query('chrom == "4"')

# %%
biomarkers = nbconfig.seurat.get_biomarkers('res.0.6').join(fbgn2chrom)

# %%
biomarkers.query('chrom == "4"')

# %% [markdown] {"toc-hr-collapsed": false}
# ## X:A Testing

# %% [markdown] {"toc-hr-collapsed": true}
# ### Ideal Wolrd

# %% [markdown]
# In and ideal world we could assume missingness was random.
#
# Then I would aggregate gene level counts to the cluster level. I would use the Mann-Whiteney U to test if X expression is less than Autosome expression. 
#
# I would then plot median X:A ratio to show difference of each cluster. 

# %%
# Aggregate gene counts to cluster level
ideal = df.groupby(['cluster', 'rep', 'chrom', 'FBgn']).UMI.sum()

# Run a mannwhitneyU test on X vs A
results = []
for (clus, rep), dd in ideal.groupby(["cluster", 'rep']):
    x_counts = dd.to_frame().query('chrom == "X"').UMI.values
    a_counts = dd.to_frame().query('chrom == "A"').UMI.values
    stat, pval = mannwhitneyu(x_counts, a_counts, alternative='less')
    results.append((clus, rep, pval))
    
ideal_results = pd.DataFrame(results, columns=['cluster', 'rep', 'p_value'])
ideal_results.cluster = ideal_results.cluster.astype('category').cat.as_ordered().cat.reorder_categories(nbconfig.short_cluster_order)
ideal_results.rep = ideal_results.rep.astype('category').cat.as_ordered().cat.reorder_categories(['rep1', 'rep2', 'rep3'])
ideal_results.set_index(['cluster', 'rep'], inplace=True)
ideal_results['significant'] = False
ideal_results.loc[ideal_results.p_value <= 0.05, 'significant'] = True

ideal_results

# %%
# Aggregate gene counts to chromsome level correcting for the number of genes
ideal_agg = ideal.groupby(['cluster', 'rep', 'chrom']).sum().div(num_genes_by_chrom, level='chrom')
ideal_agg.name = 'UMI'
ideal_agg = ideal_agg.to_frame().unstack()
ideal_agg.columns = ideal_agg.columns.droplevel(0)

# Calculate X:A ratio 
xa_ratio = ideal_agg['X'] / (ideal_agg['A'] + 0)
xa_ratio.name = 'xa'
xa_ratio = xa_ratio.to_frame().reset_index()
xa_ratio_means = xa_ratio.groupby(['cluster']).xa.mean().values

# %%
# calculate bootstrap confidence intervals for plotting
def bootstrap(dat, n_boot=1000, estimator=np.mean):
    results = np.empty(n_boot)
    for i in range(n_boot):
        results[i] = estimator(dat.sample(n=dat.shape[0], replace=True))
    return np.percentile(results, [2.5, 97.5])

results = []
for clus, dd in xa_ratio.groupby('cluster'):
    low, high = bootstrap(dd.xa)
    results.append((clus, low, high))
cluster_bootstrap = pd.DataFrame(results, columns=['cluster', 'low', 'high'])

# Merge on significant flag to add '*'
cluster_bootstrap = cluster_bootstrap.join(ideal_results.groupby('cluster').significant.any(), on='cluster')

# %%
# Plot
fig, ax = plt.subplots(figsize=plt.figaspect(1/2))
ax.plot(xa_ratio_means, color='k', zorder=-10, label='X chormosome')
sns.pointplot(x='cluster', y='xa', data=xa_ratio, errwidth=2, capsize=.2, palette=nbconfig.colors['clusters'], zorder=10, ax=ax)
ax.axhline(1, color='gray', ls=':')
ax.set_ylabel('X:A Ratio')
plt.legend(loc=2)

for i, row in cluster_bootstrap.iterrows():
    if row.significant:
        ax.text(i, row.high, '*', ha='center', va='bottom')

# %%
del ideal
del ideal_agg
del ideal_results

# %% [markdown] {"toc-hr-collapsed": false}
# ### Missingness is still problematic

# %% [markdown]
# #### Missingness by cluster

# %%
missing_per_cell = df.groupby(['cell_id', 'cluster']).missing.sum().div(num_genes_by_chrom.sum(), level='chrom')
missing_per_cell.name = 'prop_missing'

# %%
dat = missing_per_cell.reset_index()
ax = sns.boxplot('cluster', 'prop_missing', data=dat, flierprops=dict(alpha=.5), palette=nbconfig.colors['clusters'])
#plt.setp(ax.artists, edgecolor='k', facecolor='w')
#plt.setp(ax.lines, color='k');

# %%
del missing_per_cell

# %% [markdown]
# #### Missingness by cluster by chromosome

# %%
missing_per_cell_per_chrom = df.groupby(['cell_id', 'cluster', 'chrom']).missing.sum().div(num_genes_by_chrom, level='chrom')
missing_per_cell_per_chrom.name = 'prop_missing'

# %%
dat = missing_per_cell_per_chrom.reset_index()

g = sns.FacetGrid(dat, col='cluster', col_wrap=4)
g.map(
    sns.boxplot,
    'chrom',
    'prop_missing',
    order=['X', 'A', 'Y', '4'],
    flierprops=dict(alpha=.5)
)

for ax in g.axes:
    plt.setp(ax.artists, edgecolor='k', facecolor='w')
    plt.setp(ax.lines, color='k')


# %%
del missing_per_cell_per_chrom

# %% [markdown]
# #### Missingness is correlated between X and A

# %%
dat = missing_per_cell_per_chrom.to_frame().query('chrom == "X" | chrom == "A"').unstack()
dat.columns = ['prop_X_missing', 'prop_A_missing']
dat.reset_index(inplace=True)

# %%
def add_rho(color, marker, data):
    cluster = data.cluster.values[0]
    corr = spearmanr(data.prop_X_missing, data.prop_A_missing)[0]
    ax = plt.gca()
    ax.text(0.1, .9, f'r = {np.round(corr, 4)}', fontsize=12)
    
g = sns.lmplot(
    'prop_X_missing', 
    'prop_A_missing', 
    dat, 
    col='cluster', 
    col_wrap=4, 
    size=3, 
    scatter_kws=dict(alpha=.5),
)

g.set(xlim=(0, 1), ylim=(0, 1))
g.map_dataframe(add_rho)
g.set_xlabels('Prop X Missing')
g.set_ylabels('Prop A Missing')

# %%
del dat

# %% [markdown]
# ### Permutation Test

# %% [markdown]
# At the experiment level, it is clear that missingness is not random. This maybe due to technical artifacts such as dropout, or maybe related to biological processes (i.e. RNA-content of somatic cells is much smaller than germline). Therefore even a non-parametric test is not appropriate, unless we model the missingness (which is very hard). 
#
# Fortunately, at the cell level missingness appears to be somewhat random in relation to X and A expression. We have proposed using a permutation approach 

# %%


# %%
cell_ids = []
flags = []
for cell_id, dd in df.groupby('cell_id'):
    x_data = dd[dd.chrom == "X"].UMI.values
    #x_data = x_data[x_data > 0]
    a_data = dd[dd.chrom == "A"].UMI.values
    #a_data = a_data[a_data > 0]
    _, p_value = mannwhitneyu(x_data, a_data, alternative='less')
    
    if p_value <= 0.05:
        flags.append(True)
    else:
        flags.append(False)
        
    cell_ids.append(cell_id)

flag_x_lt_a = pd.Series(flags, index=pd.Index(cell_ids, name='cell_id'), name='flag_x_lt_a')

# %%
flag_x_lt_a_by_cluster = pd.concat([flag_x_lt_a, read_clusters()], axis=1, sort=True)
flag_x_lt_a_by_cluster['rep'] = flag_x_lt_a_by_cluster.index.str.extract('(?P<rep>rep\d)', expand=False)

prop_flag_by_cluster = flag_x_lt_a_by_cluster.groupby(['cluster', 'rep']).flag_x_lt_a.mean()
prop_flag_by_cluster.name = 'prop_cells_x_lt_a'

# %%
means = prop_flag_by_cluster.groupby('cluster').mean().values

# %%
fig, ax = plt.subplots(figsize=plt.figaspect(1/2))
ax.plot(means, color='k', zorder=-10)
sns.pointplot(x='cluster', y='prop_cells_x_lt_a', data=prop_flag_by_cluster.to_frame().reset_index(), errwidth=2, capsize=.2, palette=nbconfig.colors['clusters'], zorder=10, ax=ax)
ax.set_ylim(0, 1)
ax.set_ylabel('Prop Cells')
ax.set_title('Proprotion of Cells with X Depletion')

# %%


# %%


# %%
prop_missing_by_cell = df.groupby('cell_id').missing.mean()
prop_missing_by_cell.name = 'prop_missing_genes'

# %%
dat = pd.concat([prop_missing_by_cell, flag_x_lt_a.astype(int), read_clusters()], axis=1, sort=True)

# %%
dat.head()

# %%
results = smf.logit('flag_x_lt_a ~ prop_missing_genes/cluster', data=dat).fit()
results.summary2()

# %%
results = smf.logit('flag_x_lt_a ~ prop_missing_genes*cluster', data=dat).fit()
results.summary2()

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

