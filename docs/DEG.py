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

# %% [markdown] {"toc-hr-collapsed": false}
# # Differential Expression

# %%
import os
import sys
import re
from pathlib import Path
from yaml import load
import pickle

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, norm, fisher_exact, mannwhitneyu
from scipy.stats.contingency import margins
from statsmodels.stats.multitest import multipletests

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq
from larval_gonad.plotting import format_pval

# %%
sns.set_context('poster')

# %% [markdown]
# ## Data Loading

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
# Get list of gene expressed in the experiment
with open('../output/science_submission/background_fbgns.pkl', 'rb') as fh:
    bg = pickle.load(fh)

# %%
fbgn2symbol = pd.read_pickle('../output/science_submission/fbgn2symbol.pkl')
fbgn2chrom = pd.read_parquet('../output/x-to-a-wf/fbgn2chrom.parquet')

# %%
num_genes = fbgn2chrom.groupby('chrom').size().rename('num_genes')
num_genes.map(lambda x: f'{x:,}').rename('Number of Genes').to_frame()

# %%
# Get mappint of cell_id to short cluster name, remove unknown clusters
clusters = (
    pd.read_parquet('../output/scrnaseq-wf/clusters.parquet')
    .assign(cluster = lambda df: pd.Categorical(df.cluster.map(config['short_cluster_annot']), ordered=True, categories=config['short_cluster_order']))
    .dropna()
)

# %%
autosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R']

# %% [markdown] {"toc-hr-collapsed": true}
# ## Germline DEG

# %% [markdown] {"toc-hr-collapsed": true}
# ### Genes up regulated in Gonia (vs Cytes) are enriched from the X and 4th

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(fbgn2chrom)
    .assign(biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC > 0))
    .loc[:, ['chrom', 'biased']]
    .replace({
        'chrX': 'X',
        'chr2L': 'A',
        'chr2R': 'A',
        'chr3L': 'A',
        'chr3R': 'A',
        'chr4': '4',
    })
)

ct = df.groupby('chrom').biased.value_counts().unstack().loc[['X', 'A', '4']].T

_, pvalx = fisher_exact(ct[['X', 'A']], alternative='two-sided')
_, pval4 = fisher_exact(ct[['4', 'A']], alternative='two-sided')
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = df.groupby('chrom').biased.mean()[['X', 'A', '4']]
ax = dat.plot(kind='bar', width=.9, color=['darkgray', 'w', 'darkgray'], edgecolor='k', lw=1)
format_pval(ax, 0, dat['X'], pvalx)
format_pval(ax, 2, dat['4'], pval4)
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax);

# %% [markdown] {"toc-hr-collapsed": true}
# ### Genes up downregulated in E1° (vs SP) are depleted from the X and 4th

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_early.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(fbgn2chrom)
    .assign(biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC < 0))
    .loc[:, ['chrom', 'biased']]
    .replace({
        'chrX': 'X',
        'chr2L': 'A',
        'chr2R': 'A',
        'chr3L': 'A',
        'chr3R': 'A',
        'chr4': '4',
    })
)

ct = df.groupby('chrom').biased.value_counts().unstack().loc[['X', 'A', '4']].T

_, pvalx = fisher_exact(ct[['X', 'A']], alternative='two-sided')
_, pval4 = fisher_exact(ct[['4', 'A']], alternative='two-sided')
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = df.groupby('chrom').biased.mean()[['X', 'A', '4']]
ax = dat.plot(kind='bar', width=.9, color=['darkgray', 'w', 'darkgray'], edgecolor='k', lw=1)
format_pval(ax, 0, dat['X'], pvalx)
format_pval(ax, 2, dat['4'], pval4)
    
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax);

# %% [markdown] {"toc-hr-collapsed": true}
# ### Genes up downregulated in M1° (vs E1°) are depleted from the X and abscent from the 4th

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/early_vs_mid.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(fbgn2chrom)
    .assign(biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC < 0))
    .loc[:, ['chrom', 'biased']]
    .replace({
        'chrX': 'X',
        'chr2L': 'A',
        'chr2R': 'A',
        'chr3L': 'A',
        'chr3R': 'A',
        'chr4': '4',
    })
)

ct = df.groupby('chrom').biased.value_counts().unstack().loc[['X', 'A', '4']].T.fillna(0)
display(ct)

_, pvalx = fisher_exact(ct[['X', 'A']], alternative='two-sided')
_, pval4 = fisher_exact(ct[['4', 'A']], alternative='two-sided')
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = df.groupby('chrom').biased.mean()[['X', 'A', '4']]
ax = dat.plot(kind='bar', width=.9, color=['darkgray', 'w', 'darkgray'], edgecolor='k', lw=1)
format_pval(ax, 0, dat['X'], pvalx)
format_pval(ax, 2, dat['4'], pval4)
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax);

# %% [markdown] {"toc-hr-collapsed": true}
# ### Genes up downregulated in L1° (vs M1°) are depleted from X and 4th, but gene counts are very low.

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/mid_vs_late.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(fbgn2chrom)
    .assign(biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC < 0))
    .loc[:, ['chrom', 'biased']]
    .replace({
        'chrX': 'X',
        'chr2L': 'A',
        'chr2R': 'A',
        'chr3L': 'A',
        'chr3R': 'A',
        'chr4': '4',
    })
)

ct = df.groupby('chrom').biased.value_counts().unstack().loc[['X', 'A', '4']].T.fillna(0)
display(ct)

_, pvalx = fisher_exact(ct[['X', 'A']], alternative='two-sided')
_, pval4 = fisher_exact(ct[['4', 'A']], alternative='two-sided')
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = df.groupby('chrom').biased.mean()[['X', 'A', '4']]
ax = dat.plot(kind='bar', width=.9, color=['darkgray', 'w', 'darkgray'], edgecolor='k', lw=1)
format_pval(ax, 0, dat['X'], pvalx)
format_pval(ax, 2, dat['4'], pval4)
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax);

# %% [markdown] {"toc-hr-collapsed": true}
# ## Demasculinization of the X chromosome

# %% [markdown]
# ### Bulk RNA-Seq shows demasculinization of the X and the 4th.

# %%
bulk_sig = (
    pd.read_csv('../output/bulk-rnaseq-wf/deseq2_results_tcp_vs_ocp.tsv', sep='\t', index_col=0)
    .dropna()
    .assign(testis_bias = lambda df: (df.log2FoldChange >= 1) & (df.padj <= 0.01))
    .assign(ovary_bias = lambda df: (df.log2FoldChange <= -1) & (df.padj <= 0.01))
    .join(fbgn2chrom, how='outer')
    .fillna({
        'testis_bias': False,
        'ovary_bias': False,
    })
    .replace({
        'chrX': 'X',
        'chr2L': '2L',
        'chr2R': '2R',
        'chr3L': '3L',
        'chr3R': '3R',
        'chr4': '4',
    })
)

bulk_sig.loc[bulk_sig.testis_bias, 'bias'] = 'testis'
bulk_sig.loc[bulk_sig.ovary_bias, 'bias'] = 'ovary'
bulk_sig.bias = bulk_sig.bias.fillna('None')

MALE_BIAS = bulk_sig[bulk_sig.testis_bias].index

# %%
<<<<<<< HEAD
title = 'Bulk'
ct = bulk_sig.groupby('chrom').bias.value_counts().unstack().loc[['X', '2L', '2R', '3L', '3R', '4']].T
res = run_chisq(ct)
display(res.reindex(['ovary', 'None', 'testis'], level=0).loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :])

pvalx = res.loc[('testis', 'fdr q-value'), 'X']
pval4 = res.loc[('testis', 'fdr q-value'), '4']
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = ct.div(ct.sum()).T[['testis', 'None', 'ovary']]
ax = dat.plot(kind='bar', stacked=True, width=.95, color=['w', 'darkgray', 'k'], edgecolor='k', lw=1, figsize=(8, 10))
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome', title=title)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax)
format_pval(ax, 0, dat.loc['X', 'testis'], pvalx)
format_pval(ax, 5, dat.loc['4', 'testis'], pval4)
plt.legend(['Testis Biased', 'Non-Biased', 'Ovary Biased'], loc='upper left', bbox_to_anchor=(1, 1));

# %%
dat
=======
df = bulk_sig.join(fbgn2chrom, how='left').fillna('None').groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%
df = bulk_sig.join(fbgn2chrom, how='outer').fillna('None').groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
>>>>>>> 793cb92f7873d9395f5935b1250878a243e8c9f5

# %% [markdown]
# ### SP testis-biased genes do not show evidence of X or 4th demasculinization.

# %%
title = 'SP-biased'
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC > 0')
    .index.tolist()
)

<<<<<<< HEAD
_bulk_sig = bulk_sig.copy()
_bulk_sig.loc[~_bulk_sig.index.isin(BIASED), 'bias'] = 'None'
ct = _bulk_sig.groupby('chrom').bias.value_counts().unstack().loc[['X', '2L', '2R', '3L', '3R', '4']].T.fillna(0)
=======
# %%
df = bulk_sig.reindex(SP_BIASED).join(fbgn2chrom, how='left').fillna('None').groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%
df = bulk_sig.reindex(SP_BIASED).join(fbgn2chrom, how='outer').fillna('None').groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
>>>>>>> 793cb92f7873d9395f5935b1250878a243e8c9f5

res = run_chisq(ct)
display(res.reindex(['ovary', 'None', 'testis'], level=0).loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :])

pvalx = res.loc[('testis', 'fdr q-value'), 'X']
pval4 = res.loc[('testis', 'fdr q-value'), '4']
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = ct.div(ct.sum()).T[['testis', 'None', 'ovary']]
ax = dat.plot(kind='bar', stacked=True, width=.95, color=['w', 'darkgray', 'k'], edgecolor='k', lw=1, figsize=(8, 10))
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome', title=title)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax)
format_pval(ax, 0, dat.loc['X', 'testis'], pvalx)
format_pval(ax, 5, dat.loc['4', 'testis'], pval4)
plt.legend(['Testis Biased', 'Non-Biased', 'Ovary Biased'], loc='upper left', bbox_to_anchor=(1, 1));

# %% [markdown]
# ### Cyte testis-biased genes show evidence of X demasculinization and tend to have 4th demasculinization

# %%
title = 'Cyte-biased'
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

<<<<<<< HEAD
_bulk_sig = bulk_sig.copy()
_bulk_sig.loc[~_bulk_sig.index.isin(BIASED), 'bias'] = 'None'
ct = _bulk_sig.groupby('chrom').bias.value_counts().unstack().loc[['X', '2L', '2R', '3L', '3R', '4']].T.fillna(0)
=======
# %%
df = (
    bulk_sig.reindex(CYTE_BIASED)
    .join(fbgn2chrom, how='left').fillna('None')
    .assign(bias=lambda x: x.bias.replace({'ovary': "None"}))
    .groupby('chrom').bias.value_counts()
    .unstack()
    .loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
)

res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%
df = (
    bulk_sig.reindex(CYTE_BIASED)
    .join(fbgn2chrom, how='outer').fillna('None')
    .assign(bias=lambda x: x.bias.replace({'ovary': "None"}))
    .groupby('chrom').bias.value_counts()
    .unstack()
    .loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
)
>>>>>>> 793cb92f7873d9395f5935b1250878a243e8c9f5

res = run_chisq(ct)
display(res.reindex(['ovary', 'None', 'testis'], level=0).loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :])

pvalx = res.loc[('testis', 'fdr q-value'), 'X']
pval4 = res.loc[('testis', 'fdr q-value'), '4']
print(f'p-value x: {pvalx:0.4f}, p-value 4: {pval4:0.4f}')

dat = ct.div(ct.sum()).T[['testis', 'None', 'ovary']]
ax = dat.plot(kind='bar', stacked=True, width=.95, color=['w', 'darkgray', 'k'], edgecolor='k', lw=1, figsize=(8, 10))
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Chromosome', title=title)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax)
format_pval(ax, 0, dat.loc['X', 'testis'], pvalx)
format_pval(ax, 5, dat.loc['4', 'testis'], pval4)
plt.legend(['Testis Biased', 'Non-Biased', 'Ovary Biased'], loc='upper left', bbox_to_anchor=(1, 1));

# %%

# %% [markdown] {"toc-hr-collapsed": false}
# ## Gene Movement

# %% [markdown]
# https://genome.cshlp.org/content/19/5/897.long
#
# For DNA-based relocated copies, we observed a higher male-biased expression for X→A movement in comparison to A→ (Table 4), although not significant due to small sample size. However, we found the statistically significant opposite pattern for DNA-based parental copies. X→A parental copies have significant lower male-biased expression compared with A→ cases (Fisher's exact test, P = 0.03; Table 4). The same trend is observed for our RNA-based relocations (Fisher's exact test, P ≤ 0.001; Table 5).
#
# https://genome.cshlp.org/content/19/5/897/T4.expansion.html
#
# https://genome.cshlp.org/content/19/5/897/T5.expansion.html
#
# Note, that it is only significant the parental for both

# %%
assembly = config['assembly']
tag = config['tag']
pth = Path(os.environ['REFERENCES_DIR'], f'{assembly}/{tag}/fb_annotation/{assembly}_{tag}.fb_annotation')

# Create a FBgn sanitizer using secondary IDs
mapper = {}
for record in pd.read_csv(pth, sep='\t').to_records():
    mapper[record.primary_FBgn] = record.primary_FBgn
    
    try:
        for g in record.secondary_FBgn.split(','):
            mapper[g] = record.primary_FBgn
    except AttributeError:
        pass

# %%
movement = (
    pd.read_excel('../data/external/maria/dm6_ver78_genetype.new.xlsx')
    .query('gene_type == ["D", "R", "Dl", "Rl"] and m_type == "M"')
    .assign(child_chrom = lambda df: df.note.str.extract('(chr.*?)-'))
    .assign(parent_chrom = lambda df: df.note.str.extract('-(chr.*?)[:;]'))
    .assign(child_FBgn = lambda df: df.child_id.map(mapper))
    .assign(parent_FBgn = lambda df: df.parent_id.map(mapper))
    .drop(['child_id', 'parent_id', 'note', 'm_type'], axis=1)
) 

movement.loc[(movement.parent_chrom == "chrX") & movement.child_chrom.isin(autosomes), 'movement'] = 'X -> A'
movement.loc[movement.parent_chrom.isin(autosomes) & (movement.child_chrom == "chrX"), 'movement'] = 'A -> X'
movement.loc[movement.parent_chrom.isin(autosomes) & movement.child_chrom.isin(autosomes), 'movement'] = 'A -> A'

movement = movement.query('movement != "A -> X"').dropna().copy()

# %%
movement.groupby(['gene_type', 'movement']).size().unstack()

# %% [markdown]
# ### Gene Movement and Testis-biased Expression

# %% [markdown]
# #### Parent gene location is not associated with testis biased expression

# %%
ct = movement.join(bulk_sig.bias.replace({'ovary': 'Not Biased', 'None': 'Not Biased'}), on='parent_FBgn').groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# #### Child genes that moved from the X are enriched with testis biased expression

# %%
ct = movement.join(bulk_sig.bias.replace({'ovary': 'Not Biased', 'None': 'Not Biased'}), on='child_FBgn').groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# ##### DNA based movements show no association

# %%
ct = movement.join(bulk_sig.bias.replace({'ovary': 'Not Biased', 'None': 'Not Biased'}), on='child_FBgn').query('gene_type == ["D", "Dl"]').groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# ##### RNA based movements tend show testis biased enrichment

# %%
ct = movement.join(bulk_sig.bias.replace({'ovary': 'Not Biased', 'None': 'Not Biased'}), on='child_FBgn').query('gene_type == ["R", "Rl"]').groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# ### Gene Movement and Gonia Biased Expression

# %% [markdown]
# #### There is no association of the parent gene and Gonia biased expression

# %%
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC > 0')
    .index.tolist()
)

_movement = movement[['parent_FBgn', 'child_FBgn', 'movement']].copy()
_movement['bias'] = False
_movement.loc[_movement.parent_FBgn.isin(BIASED), 'bias'] = True

ct = _movement.groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# #### Child genes that moved from the X show an enrichment of gonia biased expression

# %%
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC > 0')
    .index.tolist()
)

_movement = movement[['parent_FBgn', 'child_FBgn', 'movement']].copy()
_movement['bias'] = False
_movement.loc[_movement.child_FBgn.isin(BIASED), 'bias'] = True

ct = _movement.groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# ### Gene Movement and Cyte Biased Expression

# %% [markdown]
# #### There is no association of the parent gene and Cyte biased expression

# %%
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

_movement = movement[['parent_FBgn', 'child_FBgn', 'movement']].copy()
_movement['bias'] = False
_movement.loc[_movement.parent_FBgn.isin(BIASED), 'bias'] = True

ct = _movement.groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %% [markdown]
# #### There is no association of the child gene and Cyte biased expression

# %%
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

_movement = movement[['parent_FBgn', 'child_FBgn', 'movement']].copy()
_movement['bias'] = False
_movement.loc[_movement.child_FBgn.isin(BIASED), 'bias'] = True

ct = _movement.groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %%

# %%

# %%

# %%
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

_movement = movement[['parent_FBgn', 'child_FBgn', 'movement']].copy()
_movement['bias'] = False
_movement.loc[_movement.child_FBgn.isin(BIASED), 'bias'] = True

ct = _movement.groupby('movement').bias.value_counts().unstack().fillna(0)
display(ct)
fisher_exact(ct)[1]

# %%
df = pd.pivot_table(
    (
        pd.read_parquet('../output/scrnaseq-wf/raw_by_cluster.parquet')
        .assign(cluster=lambda df: df.cluster.map(config['short_cluster_annot']))
        .query('cluster == ["SP", "M1º"]')
        .assign(flag_expressed=lambda df: df.UMI > 10)
    ),
    index='FBgn',
    columns='cluster',
    values='flag_expressed'
)


# %%
ct = movement.join(df, on='child_FBgn').loc[:, ['movement', 'SP', 'M1º']].dropna().groupby('movement').sum()
display(ct)
run_chisq(ct)

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

# %%

# %%

# %% [markdown]
# ### Gene movement and missingness

# %%
tpm = (
    pd.read_parquet('../output/scrnaseq-wf/tpm.parquet')
    .assign(cluster=lambda df: df.cluster.map(config['short_cluster_annot']))
    .query('cluster == "M1º"')
    .TPM
    .rename('M1')
)

# %%
BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .assign(gonia_biased=lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC > 0))
    .assign(cyte_biased=lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC < 0))
    .assign(gonia_cyte_bias=lambda df: df[['gonia_biased', 'cyte_biased']].idxmax(axis=1))
    .gonia_cyte_bias
)

# %%
df = movement.join(tpm, on='child_FBgn').join(BIASED, on='child_FBgn')

# %%
df.loc[(df.M1 == 0), 'gonia_cyte_bias'] = 'Not Expressed'
df.loc[(df.M1.isnull()), 'gonia_cyte_bias'] = 'Missing'
df.fillna({'gonia_cyte_bias': 'Not Biased'}, inplace=True)

# %%
ct = df.groupby('movement').gonia_cyte_bias.value_counts().unstack()

# %%
ax = ct.div(ct.sum(axis=1), axis='rows').loc[['X -> A', 'A -> A'], ['gonia_biased', 'cyte_biased', 'Not Biased', 'Not Expressed', 'Missing']].plot(kind='bar', stacked=True, width=.95, figsize=(8, 10), edgecolor='k', lw=1)
plt.legend(['Gonia Biased', 'Cyte Biased', 'Not Biased', 'Not Expressed', 'Missing'], loc='upper left', bbox_to_anchor=[1, 1])
ax.set(ylim=(0, 1), ylabel='Proportion of Genes', xlabel='Movement')
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
sns.despine(ax=ax)

# %%
res = run_chisq(ct)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%

# %%
ct2 = pd.concat([ct[['Not Expressed', 'Missing']].sum(axis=1), ct[['Not Biased', 'cyte_biased', 'gonia_biased']].sum(axis=1)], axis=1).rename({0: 'Off', 1: 'On'}, axis=1)

display(ct2)
pval = fisher_exact(ct2, alternative='two-sided')[1]
print(f"fisher's pval: {pval:0.4f}")

# %%
ct2 = pd.concat([ct['Missing'], ct[['Not Expressed', 'Not Biased', 'cyte_biased', 'gonia_biased']].sum(axis=1)], axis=1).rename({0: 'Not Missing'}, axis=1)

display(ct2)
pval = fisher_exact(ct2, alternative='two-sided')[1]
print(f"fisher's pval: {pval:0.4f}")

# %%

# %%

# %%

# %% [markdown]
# ## Repeat what Maria did.

# %%
maria = (
    pd.read_csv('../output/notebook/2019-04-16_movement_data.csv', index_col=0)
    .assign(movement=lambda df: pd.Categorical(
        df[['moved_x_to_a', 'moved_a_to_x', 'moved_a_to_a']].idxmax(axis=1), 
        ordered=True,
        categories=['moved_x_to_a', 'moved_a_to_a', 'moved_a_to_x']
    ))
    .fillna({
        'biomarker_cluster': "None",
        'bias_gonia_vs_mid_child': "None",
        'bias_gonia_vs_mid_parent': "None",
        'SP_child': 0,
        'M1_child': 0,
        'SP_parent': 0,
        'M1_parent': 0
    })
    .assign(log_M1_child=lambda df: np.log(df.M1_child + 1))
)

# %%
maria.groupby('movement').M1_child.size()

# %%
fig = plt.figure(figsize=(10, 10))
ax = sns.boxplot('movement', 'M1_child', data=maria, color='w')
plt.setp(ax.artists, edgecolor='k', facecolor='w', lw=2)
plt.setp(ax.lines, color='k', lw=2)
ax.set_ylim(None, 1000);

# %%
fig = plt.figure(figsize=(10, 10))
ax = sns.boxplot('movement', 'log_M1_child', data=maria, color='w')
plt.setp(ax.artists, edgecolor='k', facecolor='w', lw=2)
plt.setp(ax.lines, color='k', lw=2)
#ax.set_ylim(None, 1000);
ax.set_title('NA droped')

# %%

# %%

# %%

# %%

# %%

# %%

# %%
mannwhitneyu(maria.query('movement == "moved_x_to_a"').M1_child.dropna(), maria.query('movement == "moved_a_to_a"').M1_child.dropna(), alternative="two-sided")

# %%

# %%
fig = plt.figure(figsize=(10, 10))
ax = sns.boxplot('movement', 'M1_parent', data=maria, color='w')
plt.setp(ax.artists, edgecolor='k', facecolor='w', lw=2)
plt.setp(ax.lines, color='k', lw=2)
ax.set_ylim(None, 1000);

# %%
mannwhitneyu(maria.query('movement == "moved_x_to_a"').M1_parent.dropna(), maria.query('movement == "moved_a_to_a"').M1_parent.dropna(), alternative="two-sided")

# %%

# %%

# %%

# %%
maria.head(20)

# %%
dat = (maria.groupby('movement').bias_gonia_vs_mid_child.value_counts().div( maria.groupby('movement').size()) * 100).rename('prop').to_frame().reset_index()

# %%
fig = plt.figure(figsize=(10, 10))
sns.barplot('movement', 'prop', hue='bias_gonia_vs_mid_child', data=dat)

# %%
fisher_exact(dat.query('movement == ["moved_x_to_a", "moved_a_to_a"] and bias_gonia_vs_mid_child == ["M1", "None"]').set_index(['movement', 'bias_gonia_vs_mid_child']).unstack(), alternative='two-sided')

# %%
dat.query('movement == ["moved_x_to_a", "moved_a_to_a"] and bias_gonia_vs_mid_child == ["M1", "None"]').set_index(['movement', 'bias_gonia_vs_mid_child']).unstack()

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
movement = (
    pd.read_excel('../data/external/maria/dm6_ver78_genetype.new.xlsx')
    .query('gene_type == ["D", "R", "Dl", "Rl"] and m_type == "M"')
    .assign(child_chrom = lambda df: df.note.str.extract('(chr.*?)-'))
    .assign(parent_chrom = lambda df: df.note.str.extract('-(chr.*?)[:;]'))
    .assign(FBgn = lambda df: df.child_id.map(mapper))
    .set_index("FBgn")
    .drop(['child_id', 'parent_id', 'note', 'm_type'], axis=1)
    .dropna()
) 

movement.loc[(movement.parent_chrom == "chrX") & movement.child_chrom.isin(autosomes), 'movement'] = 'X -> A'
movement.loc[movement.parent_chrom.isin(autosomes) & (movement.child_chrom == "chrX"), 'movement'] = 'A -> X'
movement.loc[movement.parent_chrom.isin(autosomes) & movement.child_chrom.isin(autosomes), 'movement'] = 'A -> A'
movement.movement = pd.Categorical(movement.movement, ordered=True, categories=['X -> A', 'A -> X', 'A -> A'])

# %%
m1_tpm = (
    pd.read_parquet('../output/scrnaseq-wf/tpm.parquet')
    .assign(cluster=lambda df: df.cluster.map(config['short_cluster_annot']))
    .query('cluster == "M1º"')
)

# %%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 20))
dat = movement.join(m1_tpm).assign(log_TPM=lambda df: np.log10(df.TPM + 1))
sns.boxplot('movement', 'TPM', data=dat, ax=ax1, showfliers=True)
sns.boxplot('movement', 'log_TPM', data=dat, ax=ax2, showfliers=False)
ax1.set_ylim(None, 1010)

# %%
dat.groupby('movement').TPM.median()

# %%
dat.groupby('movement').TPM.size()

# %%
dat.query('TPM > 0').groupby('movement').TPM.size()

# %%
from scipy.stats import mannwhitneyu

# %%
mannwhitneyu(dat.query('movement == "X -> A"').TPM, dat.query('movement == "A -> A"').TPM)

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

# %%
movement = (
    pd.read_excel('../data/external/maria/dm6_ver78_genetype.new.xlsx')
    .query('gene_type == ["D", "R", "Dl", "Rl"] and m_type == "M"')
    .assign(child_chrom = lambda df: df.note.str.extract('(chr.*?)-'))
    .assign(parent_chrom = lambda df: df.note.str.extract('-(chr.*?)[:;]'))
    .assign(FBgn = lambda df: df.child_id.map(mapper))
    .set_index("FBgn")
    .drop(['child_id', 'parent_id', 'note', 'm_type'], axis=1)
    .dropna()
) 

movement.loc[(movement.parent_chrom == "chrX") & movement.child_chrom.isin(autosomes), 'movement'] = 'X -> A'
movement.loc[movement.parent_chrom.isin(autosomes) & (movement.child_chrom == "chrX"), 'movement'] = 'A -> X'
movement.loc[movement.parent_chrom.isin(autosomes) & movement.child_chrom.isin(autosomes), 'movement'] = 'A -> A'
movement.movement = pd.Categorical(movement.movement, ordered=True, categories=['X -> A', 'A -> X', 'A -> A'])

# %%
movement.head()

# %%
norm = (
    pd.read_csv('../output/scrnaseq-wf/scrnaseq_combine_force/normalized_read_counts.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(movement, how='inner')
    .reset_index()
    .drop(['parent_chrom', 'child_chrom'], axis=1)
    .melt(id_vars=['FBgn', 'movement', 'gene_type'], var_name='cell_id', value_name='norm')
    .set_index('FBgn')
    .assign(rep=lambda df: df.cell_id.str.extract('(rep\d)', expand=False))
    .join(clusters, on='cell_id')
)

# %%
fig = plt.figure(figsize=(20, 10))
ax = sns.boxplot('movement', 'norm', hue='cluster', data=norm.query('norm > 0'), order=['X -> A', 'A -> A'])
plt.legend(loc='upper left', bbox_to_anchor=[1, 1])
ax.set_axisbelow(True)
ax.grid(axis='y')

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

# %%

# %%

# %%

# %%

# %%
