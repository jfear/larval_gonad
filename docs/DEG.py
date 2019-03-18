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

# %% [markdown]
# # Differential Expression

# %%
import os
import sys
import re
from pathlib import Path
from yaml import load

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, norm
from scipy.stats.contingency import margins

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# %% [markdown]
# ## Data

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
fbgn2symbo = pd.read_pickle('../output/science_submission/fbgn2symbol.pkl')
fbgn2chrom = pd.read_parquet('../output/x-to-a-wf/fbgn2chrom.parquet')

# %%
# Get mappint of cell_id to short cluster name, remove unknown clusters
clusters = (
    pd.read_parquet('../output/scrnaseq-wf/clusters.parquet')
    .assign(cluster = lambda df: pd.Categorical(df.cluster.map(config['short_cluster_annot']), ordered=True, categories=config['short_cluster_order']))
    .dropna()
)

# %% [markdown]
# ## Are cluster biomarkers depleted on of X-linked and 4th-linked genes and enriched for Y-linked genes in the germline?

# %%
# Get biomarkers
resolution = config['resolution']
biomarkers = (
    pd.read_csv(f'../output/scrnaseq-wf/scrnaseq_combine_force/biomarkers_{resolution}.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .query('p_val_adj <= 0.05')
    .assign(cluster = lambda df: pd.Categorical(df.cluster.map(config['short_cluster_annot']), ordered=True, categories=config['short_cluster_order']))
    .dropna()
)

# %%
# Number of biomarker genes by cluster and chromosome
df = (
    biomarkers.join(fbgn2chrom)
    .groupby(['cluster', 'chrom']).size()
    .unstack().fillna(0).drop('chrM', axis=1)
    .assign(autosome=lambda df: df[['chr2L', 'chr2R', 'chr3L', 'chr3R']].sum(axis=1))
    .drop(['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrY'], axis=1)
    .reset_index()
    .assign(cluster = lambda df: df.cluster.astype(str))
    .set_index('cluster')
    .T
    .assign(soma = lambda df: df[['EC', 'MC', 'LC', 'TE', 'PC']].sum(axis=1))
    .drop(['EC', 'MC', 'LC', 'TE', 'PC'], axis=1)
)


# %%
def adjusted_residuals(observed, expected):
    resid = (observed - expected) / np.sqrt(expected)
    n = observed.sum().sum()
    rsum, csum = margins(observed)
    v = csum * rsum * (n - rsum) * (n - csum) / n**3
    return (observed - expected) / np.sqrt(v)

# %%
# Chi-squre test with adjusted standardized residuals pos hoc test
stat, pval, degrees, expected = chi2_contingency(df)
print(f'𝛘^2: {stat:,.4f}, p-value: {pval:,.4f}, df: {degrees:,})')
display(df)
resid = adjusted_residuals(df, expected)
display(resid)
print(f'Adjusted Residuals cutoff: N(0, 1)_{{1-.05/15/2}} = {norm.ppf(1 - (.05 / 15) / 2)}')

# %%
# Ylinked bio-marker genes
biomarkers.join(fbgn2chrom).query('chrom == "chrY"')

# %% [markdown]
# ## New genes coming on

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .query('p_val_adj <= 0.05')
    .assign(gonia_biased = lambda df: df.avg_logFC > 0)
    .join(fbgn2chrom)
)
df.gonia_biased.value_counts().rename({True: 'gonia_basied', False: 'cyte_biased'}).rename('Gona vs Cytes').map(lambda x: f'{x:,}').to_frame()

# %%
df2 = (
    df
    .groupby(['gonia_biased', 'chrom']).size()
    .unstack().fillna(0).drop('chrM', axis=1)
    .assign(autosome=lambda x: x[['chr2L', 'chr2R', 'chr3L', 'chr3R']].sum(axis=1))
    .drop(['chr2L', 'chr2R', 'chr3L', 'chr3R'], axis=1)
    .rename({True: 'gonia_basied', False: 'cyte_biased'})
    .rename_axis('Gona vs Cytes')
)

# %%
# Chi-squre test with adjusted standardized residuals pos hoc test
stat, pval, degrees, expected = chi2_contingency(df2)
print(f'𝛘^2: {stat:,.4f}, p-value: {pval:,.4f}, df: {degrees:,})')
display(df2)
resid = adjusted_residuals(df2, expected)
display(resid)
print(f'Adjusted Residuals cutoff: N(0, 1)_{{1-.05/8/2}} = {norm.ppf(1 - (.05 / 8) / 2)}')

# %%
print('\n'.join(df[df.gonia_biased].index.unique().tolist()))

# %%
