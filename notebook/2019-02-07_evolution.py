# -*- coding: utf-8 -*-
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
# # Parse New Gene Table

# %% [markdown]
# **from:** Maria D. Vibranovski
#
# Here attached is a list from Yong Zhang group based on our paper from 2010. But this is a still not published updated version that he shared with me but you can use.
#
# If you need details about the columns, please look at https://genome.cshlp.org/content/suppl/2010/08/27/gr.107334.110.DC1/SupplementalMaterial.pdf  table 2a.
#
# But mainly, what you need to select is the child genes with:
#
# gene_type = D or R or DL or RL
# m_type= M
# note that contains "chrX-"
#
# D and R stands for DNA-based Duplication and RNA-based duplication
# L means that the assignment of the parental genes is less reliable.
# M indicates that is between chromosome movement.
#
# Hope it helps. If you need I can parse for you. please, do not hesitate to ask. But I thought you would prefer a complete list where you can look at subsets.
#
# cheers
#
# Maria
#

# %%
import os
import sys
from pathlib import Path
import re

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from scipy.stats.contingency import margins
import statsmodels.formula.api as smf

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
sys.path.insert(0, '../lib')
from larval_gonad.notebook import Nb
from larval_gonad.plotting import make_figs
from larval_gonad.config import memory

# Setup notebook
nbconfig = Nb.setup_notebook()

# %%
def adjusted_residuals(observed, expected):
    resid = (observed - expected) / np.sqrt(expected)
    n = observed.sum().sum()
    rsum, csum = margins(observed)
    v = csum * rsum * (n - rsum) * (n - csum) / n**3
    return (observed - expected) / np.sqrt(v)

# %% [markdown]
# ## Import data from Maria

# %% [markdown]
# ## FBgn sanitizer

# %% [markdown]
# I don't know where these FBgns are from, so I need to sanitize them to my current annotation.

# %%
assembly = nbconfig.assembly
tag = nbconfig.tag
pth = Path(os.environ['REFERENCES_DIR'], f'{assembly}/{tag}/fb_annotation/{assembly}_{tag}.fb_annotation')

# Create an FBgn 
mapper = {}

for record in pd.read_csv(pth, sep='\t').to_records():
    mapper[record.primary_FBgn] = record.primary_FBgn
    
    try:
        for g in record.secondary_FBgn.split(','):
            mapper[g] = record.primary_FBgn
    except AttributeError:
        pass

# %%
autosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R']

# %%
movement = (
    pd.read_excel('../data/external/maria/dm6_ver78_genetype.new.xlsx')
    .query('gene_type == ["D", "R", "Dl", "Rl"] and m_type == "M"')
    .assign(child_chrom = lambda df: df.note.str.extract('(chr.*?)-'))
    .assign(parent_chrom = lambda df: df.note.str.extract('-(chr.*?)[:;]'))
    .assign(FBgn = lambda df: df.child_id.map(mapper))
    .assign(parent_FBgn = lambda df: df.parent_id.map(mapper))
    .drop(['child_id', 'parent_id', 'note', 'm_type'], axis=1)
    .dropna()
    .set_index('FBgn')
    .assign(x_to_a = lambda df: (df.parent_chrom == 'chrX') & df.child_chrom.isin(autosomes))
    .assign(a_to_a = lambda df: df.parent_chrom.isin(autosomes) & df.child_chrom.isin(autosomes))
    .assign(a_to_x = lambda df: df.parent_chrom.isin(autosomes) & (df.child_chrom == 'chrX'))
    .query('x_to_a | a_to_a | a_to_x')
)

movement.head()

# %%
germline = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t')
    .assign(FBgn = lambda df: df.primary_FBgn)
    .assign(Spermatogonia = lambda df: df.avg_logFC > 0)
    .assign(Spermatocyte = lambda df: df.avg_logFC < 0)
    .set_index('FBgn')
    .iloc[:, -2:]
    .idxmax(axis=1)
    .reindex(background_genes)
    .fillna('Not Different')
)
germline.name = 'Expression Bias'

# %%
ct = pd.crosstab(germline, movement.x_to_a)
display(ct)
_, pval, _, expected = chi2_contingency(ct)
#_, pval_fisher = fisher_exact(ct)
#print(pval, pval_fisher)
pval

# %%
ct = pd.crosstab(germline, movement.a_to_a)
display(ct)
_, pval, _, expected = chi2_contingency(ct)
print(pval)

# %%
ct = pd.crosstab(germline, movement.a_to_x)
display(ct)
_, pval, _, expected = chi2_contingency(ct)
print(pval)

# %%
_movement = movement.iloc[:, -3:].idxmax(axis=1)
_movement = _movement.map({
    'a_to_a': 'A → A',
    'a_to_x': 'A → X',
    'x_to_a': 'X → A',
})

# %%
ct =pd.crosstab(germline, _movement)
ct.columns.name = ''
display(HTML('<h4>Observed Counts</h4>'))
display(ct)
_, pval, _, expected = chi2_contingency(ct)
display(HTML(f'<b>chi<sup>2</sup> p-value</b>: {np.round(pval, 4)}'))
display(HTML('<hr>'))
display(HTML('<h4>Adjusted Residuals</h4>'))
adjusted_residuals(ct, expected)

# %%
_dat = ct.reset_index().melt(id_vars='Expression Bias', var_name='Movement', value_name='Number of Genes')

# %%
sns.barplot(
    x='Movement',
    y='Number of Genes',
    hue='Expression Bias',
    data=_dat,
)

# %%


# %%


# %%


# %%


# %%


# %%
_dat = pd.concat([germline, movement.x_to_a], axis=1, join='inner', sort=True)
_dat.head()
_dat.columns = ['expression', 'x_to_a']
_dat.x_to_a = _dat.x_to_a.astype(int)

# %%
model = smf.logit('x_to_a ~ expression', _dat).fit()

# %%
model.summary2()

# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%

