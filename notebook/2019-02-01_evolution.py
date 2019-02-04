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
# gene_type = D or R or DI or RI
# m_type= M
# note that contains "chrX-"
#
# D and R stands for DNA-based Duplication and RNA-based duplication
# I means that the assignment of the parental genes is less reliable.
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
movement = (
    pd.read_excel('../data/external/maria/dm6_ver78_genetype.new.xlsx')
    .query('gene_type == ["D", "R", "Dl", "Rl"] and m_type == "M"')
    .loc[:, ["child_id", "parent_id", "note"]]
    .assign(child_chrom = lambda df: df.note.str.extract('(chr.*?)-'))
    .assign(parent_chrom = lambda df: df.note.str.extract('-(chr.*?)[:;]'))
    .assign(FBgn = lambda df: df.child_id.map(mapper))
    .assign(parent_FBgn = lambda df: df.parent_id.map(mapper))
    .drop(['child_id', 'parent_id', 'note'], axis=1)
    .dropna()
    .set_index('FBgn')
    .assign(moved_from_x = lambda df: df.parent_chrom == 'chrX')
    .assign(moved_from_2L = lambda df: df.parent_chrom == 'chr2L')
    .assign(moved_from_2R = lambda df: df.parent_chrom == 'chr2R')
    .assign(moved_from_3L = lambda df: df.parent_chrom == 'chr3L')
    .assign(moved_from_3R = lambda df: df.parent_chrom == 'chr3R')
)

movement.head()

# %%
movement.parent_chrom.value_counts()

# %%
movement.child_chrom.value_counts()

# %%
gonia_vs_cytes = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t')
    .assign(FBgn = lambda df: df.primary_FBgn)
    .assign(gonia_bias = lambda df: df.avg_logFC > 0)
    .assign(cyte_bias = lambda df: df.avg_logFC < 0)
    .set_index('FBgn')
    .loc[:, ['gonia_bias', 'cyte_bias']]
)

# %%
gonia_vs_cytes.gonia_bias.value_counts()

# %%
dat = gonia_vs_cytes.join(movement, how='inner')

# %%
ct = pd.crosstab(dat.moved_from_x, dat.cyte_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_2L, dat.cyte_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_2R, dat.cyte_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_3L, dat.cyte_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_3R, dat.cyte_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%


# %%
germ_vs_soma = (
    pd.read_csv('../output/scrnaseq-wf/germcell_soma_deg/germ_vs_soma.tsv', sep='\t')
    .assign(FBgn = lambda df: df.primary_FBgn)
    .assign(germ_bias = lambda df: df.avg_logFC > 0)
    .assign(soma_bias = lambda df: df.avg_logFC < 0)
    .set_index('FBgn')
    .loc[:, ['germ_bias', 'soma_bias']]
)

# %%
germ_vs_soma.germ_bias.value_counts()

# %%
dat = germ_vs_soma.join(movement, how='inner')

# %%
ct = pd.crosstab(dat.moved_from_x, dat.germ_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
_, pval = fisher_exact(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_2L, dat.germ_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_2R, dat.germ_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_3L, dat.germ_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

# %%
ct = pd.crosstab(dat.moved_from_3R, dat.germ_bias)
print('ct')
display(ct)

_, pval, _, exptedted = chi2_contingency(ct)
print('pvale=', pval)

expected = pd.DataFrame(expected, columns=ct.columns, index=ct.index)
print('--------------------------------------------------------------------------------')
print('expected')
display(expected)

print('--------------------------------------------------------------------------------')
print('adjusted residuals')
resid = (ct - expected) / np.sqrt(expected)
n = ct.sum().sum()
rsum, csum = margins(ct)
v = csum * rsum * (n - rsum) * (n - csum) / n**3
(ct - expected) / np.sqrt(v)

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

