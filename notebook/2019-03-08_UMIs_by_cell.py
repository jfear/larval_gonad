# -*- coding: utf-8 -*-
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

# %%
import os
import sys
import re
from pathlib import Path
from itertools import combinations

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

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
    pd.read_parquet('../output/science_submission/clusters.parquet')
    .query('cluster != "UNK"')
    .assign(
        cluster=lambda df: (
            df.cluster.astype('category')
            .cat.as_ordered()
            .cat.reorder_categories(nbconfig.short_cluster_order)
        )
    )
)

# %%
df = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
    .sum()
    .rename('total_UMI')
    .rename_axis('cell_id')
    .to_frame()
    .assign(rep = lambda df: df.index.str.extract('(rep\d)', expand=False))
    .assign(log_total_UMI=lambda df: np.log10(df.total_UMI))
    .join(clusters)
)

# %%
sns.violinplot('rep', 'log_total_UMI', data=df, scale='count', inner='quartile')
plt.title('Total UMI (log) by Replicate');

# %%
fig = plt.figure(figsize=(plt.figaspect(1/2)))
sns.violinplot('cluster', 'log_total_UMI', data=df, scale='count', inner='quartile', palette=nbconfig.colors['clusters'])
plt.title('Total UMI (log) by Cluster');

# %%
sns.barplot('cluster', 'cnt', data=df.groupby(['rep', 'cluster']).size().rename('cnt').to_frame().reset_index(), hue='rep')
plt.title('Number of Cells Per Replicate by Cluster');

# %%
for c1, c2 in combinations(nbconfig.short_cluster_order, 2):
    _, pval = mannwhitneyu(
        df.query(f'cluster == "{c1}"').total_UMI.values,
        df.query(f'cluster == "{c2}"').total_UMI.values,
        alternative='two-sided'
    )
    print(f'{c1:<5} {c2:<5} {pval < 0.001}')

# %%
df.groupby('cluster').total_UMI.median().map(lambda x: f'{x:,.2f}')

# %%

# %%
df = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
    .pipe(lambda df: df > 0)
    .sum()
    .rename('expressed_genes')
    .rename_axis('cell_id')
    .to_frame()
    .assign(log_expressed_genes=lambda df: np.log10(df.expressed_genes + 1))
    .join(clusters, how='inner')
)

# %%
fig = plt.figure(figsize=plt.figaspect(1/2))
sns.violinplot('cluster', 'expressed_genes', data=df, palette=nbconfig.colors['clusters'], scale='count', inner='quartile')
plt.title('Expressed Genes by Cluster')

# %%
for c1, c2 in combinations(nbconfig.short_cluster_order, 2):
    _, pval = mannwhitneyu(
        df.query(f'cluster == "{c1}"').expressed_genes.values,
        df.query(f'cluster == "{c2}"').expressed_genes.values,
        alternative='two-sided'
    )
    print(f'{c1:<5} {c2:<5} {pval < 0.001}')

# %%
df.groupby('cluster').expressed_genes.mean().map(lambda x: f'{x:,.2f}')

# %%

# %%
biomarkers = (
    pd.read_csv('../output/scrnaseq-wf/scrnaseq_combine_force/biomarkers_res.0.6.tsv', sep='\t', index_col=0, usecols=['primary_FBgn', 'cluster'])
    .cluster
    .rename_axis('FBgn')
    .map(nbconfig.short_cluster_annot)
    .pipe(lambda x: x[x != 'UNK'])
    .astype('category')
    .cat.as_ordered()
    .cat.reorder_categories(nbconfig.short_cluster_order)
    .to_frame()
)

# %%
unique = biomarkers.groupby('FBgn').size().pipe(lambda x: x[x == 1]).index

# %%
l1 = biomarkers[biomarkers.index.isin(unique)].query('cluster == "L1º"').index

# %%
for g in l1:
    print(f'{g:<12} {nbconfig.fbgn2symbol[g]}')

# %%
l1.tolist()

# %%
df = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet').T
    .join(clusters)
    .groupby('cluster')
    .sum()
    .T
)
df

# %%
binned = df.apply(partial(pd.qcut, q=4, labels=['lowest', 'low', 'high', 'highest'], duplicates='drop'), axis=0)

# %%
binned[binned.index.isin(l1)]

# %%
binned[['M1º', 'L1º']]

# %%
g = sns.clustermap(df.corr(method='spearman'), annot=True)
plt.setp(g.ax_heatmap.get_xticklabels(), fontsize=18)
plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=18, rotation=0)

# %%

# %%

# %%
bd = binned[~((binned == 'lowest').all(axis=1) |  (binned == 'low').all(axis=1) | (binned == 'high').all(axis=1) |  (binned == 'highest').all(axis=1))]

# %%
bd

# %%
df.apply(partial(pd.cut, bins=2, labels=['low', 'high'], duplicates='drop'), axis=1)

# %%
df.loc['FBgn0031086', :]

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
X = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet').T
    .join(clusters)
    .query('cluster == ["M1º", "L1º"]')
    .drop('cluster', axis=1)
)

Y = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet').T
    .join(clusters)
    .query('cluster == ["M1º", "L1º"]')
    .cluster
    .values
    .tolist()
)

# %%
cell_ids = X.index
genes = X.columns
X = X.values

# %%
X

# %%
from sklearn.ensemble import RandomForestClassifier

# %%
RF = RandomForestClassifier(n_jobs=10)

# %%
RF.fit(X, Y)

# %%
cand = pd.Series(RF.feature_importances_, index=genes).sort_values(ascending=False).head(50).index
list(map(lambda x: (x, nbconfig.fbgn2symbol[x]), cand))

# %%
df.reindex(cand)

# %%

# %%
(df['M1º'] > df['L1º']).sum()

# %%
(df['M1º'] < df['L1º']).sum()

# %%
df.loc[((df['L1º'] - df['M1º']) > 100), ['M1º', 'L1º']]

# %%
_323.index.map(nbconfig.fbgn2symbol)

# %%
l1.map(nbconfig.fbgn2symbol)

# %%

# %%

# %%
pd.read_csv('../output/scrnaseq-wf/germcell_deg/mid_vs_late.tsv', sep='\t').query('p_val_adj <= 0.01').sort_values('avg_logFC').set_index('primary_FBgn').reindex(male_sterile).dropna()

# %%

# %%

# %%
with open('/home/fearjm/Downloads/male_sterile.txt') as fh:
    male_sterile = fh.read().strip().split('\n')

# %%
