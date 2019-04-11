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
import pickle

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, norm, fisher_exact
from scipy.stats.contingency import margins
from statsmodels.stats.multitest import multipletests

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq

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

# %% [markdown]
# ## Genes up regulated in primary spermatocytes are depleted from the X and 4th

# %% [markdown] {"toc-hr-collapsed": true}
# ### Hypothesis
#
# If MSCI is TRUE, then we would expect that genes important for spermiogenesis will be depleted from the X.
#
# ### Method
#
# Using differential expression between spermatogonia and either all primary spermatocytes or M1° primary spermatocytes and look for chromosomal distributions for genes up regulated in spermatocytes.

# %% [markdown]
# #### Gonia vs All Primary Spermatocytes
# **TRUE: Cyte-biased genes are depleted on the X and 4th**

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(fbgn2chrom)
    .assign(SP_biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC > 0))
    .assign(cyte_biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC < 0))
    .assign(no_difference = lambda df: ~df.SP_biased & ~df.cyte_biased)
    .groupby('chrom')[['SP_biased', 'cyte_biased']].sum()
    .T
    .assign(autosomes = lambda df: df[['chr2L', 'chr2R', 'chr3L', 'chr3R']].sum(axis=1))
    .loc[:, ['chrX', 'chr4', 'autosomes']]
)
print("Fisher's exact test: chrX: {}, chr4: {}".format(fisher_exact(df[['chrX', 'autosomes']])[1], fisher_exact(df[['chr4', 'autosomes']])[1]))
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %% [markdown]
# #### Gonia vs M1° Primary spermatocytes
#
# **TRUE: Cyte-biased genes are depleted on the X and 4th**

# %%
df = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_mid.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .join(fbgn2chrom)
    .assign(SP_biased = lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC > 0))
    .assign(**{'M1°_biased': lambda df: (df.p_val_adj <= 0.01) & (df.avg_logFC < 0)})
    .assign(no_difference = lambda df: ~df.SP_biased & ~df['M1°_biased'])
    .groupby('chrom')[['SP_biased', 'M1°_biased']].sum()
    .T
    .assign(autosomes = lambda df: df[['chr2L', 'chr2R', 'chr3L', 'chr3R']].sum(axis=1))
    .loc[:, ['chrX', 'chr4', 'autosomes']]
)
print("Fisher's exact test: chrX: {}, chr4: {}".format(fisher_exact(df[['chrX', 'autosomes']])[1], fisher_exact(df[['chr4', 'autosomes']])[1]))
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%

# %% [markdown]
# ## Demasculinization of the X chromosome
#
# Previous studies have shown that testis biased genes are depleted from the X chromosome. These studies used microarray and bulk RNA-Seq from adult testis.

# %% [markdown]
# ### Verification of X demasculinization in bulk RNA-Seq from larval testes.
#
# #### Hypothesis
#
# Testis bias genes will be depleted from the X.
#
# #### Method
#
# * Compare bulk RNA-Seq from pools of larval testes and larval ovaries. Here we used testes and ovaries after removing fatbody with papain. 
# * We selected genes with an adjusted p-value <= 0.01 and greater than a 2-fold change (i.e., |log2FC| > 1). 
#
# #### Results
#
# * Roughly 60% of the transcriptome is differentially expressed, with nearly 40% of genes showing male-biased expression.
# * Genes with male-baised expression are depleted on the X and 4th.

# %%
bulk_sig = (
    pd.read_csv('../output/bulk-rnaseq-wf/deseq2_results_tcp_vs_ocp.tsv', sep='\t', index_col=0)
    .dropna()
    .assign(testis_bias = lambda df: (df.log2FoldChange >= 1) & (df.padj <= 0.01))
    .assign(ovary_bias = lambda df: (df.log2FoldChange <= -1) & (df.padj <= 0.01))
)

bulk_sig.loc[bulk_sig.testis_bias, 'bias'] = 'testis'
bulk_sig.loc[bulk_sig.ovary_bias, 'bias'] = 'ovary'
bulk_sig.bias = bulk_sig.bias.fillna('None')

MALE_BIAS = bulk_sig[bulk_sig.testis_bias].index

# %%
df = bulk_sig.join(fbgn2chrom, how='outer').fillna('None').groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%

# %% [markdown]
# ### Spermatogonia-testis-biased genes do not show evidence of X demasculinization.
#
# #### Hypothesis
#
# If there is really demasculinization of the X, then genes that are testis-biased and SP-biased would be depleted from the X.
#
# #### Methods
#
# First I selected SP-biased gene as genes up regulated in SP vs Cytes. 
# Next I selected testis-biased genes from bulk RNA-Seq comparison. 
# Testis-biased or ovary-biased genes that were not SP-biased were set to "None". 
# I then look at chromosomal distributions.
#
# #### Results
#
# In intersection of SP-biased and testis-biased there was no evidence of X demasculinization. 
# In fact, the X-linked and 4th-linked genes tended to be enriched

# %%
SP_BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC > 0')
    .index.tolist()
)

# %%
df = bulk_sig.reindex(SP_BIASED).join(fbgn2chrom, how='outer').fillna('None').groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%

# %% [markdown]
# ### Spermatocyte-testis-biased genes do show evidence of X demasculinization.
#
# #### Hypothesis
#
# If there is really demasculinization of the X, then genes that are testis-biased and Cyte-biased would be depleted from the X.
#
# #### Methods
#
# First I selected Cyte-biased gene as genes up regulated in SP vs Cytes. 
# Next I selected testis-biased genes from bulk RNA-Seq comparison. 
# Testis-biased genes that were not SP-biased were set to "None".
# The numbers are very small so I set all ovary-biased genes to "None" as well.
# I then look at chromosomal distributions.
#
# #### Results
#
# In intersection of Cyte-biased and testis-biased show significant evidence of X demasculinization. 

# %%
CYTE_BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

# %%
df = (
    bulk_sig.reindex(CYTE_BIASED)
    .join(fbgn2chrom, how='outer').fillna('None')
    .assign(bias=lambda x: x.bias.replace({'ovary': "None"}))
    .groupby('chrom').bias.value_counts()
    .unstack()
    .loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
)

res = run_chisq(df)
res.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]

# %%

# %%

# %%

# %% [markdown] {"toc-hr-collapsed": true}
# ## Do male-biased genes show movement?

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
    .assign(child_FBgn = lambda df: df.child_id.map(mapper))
    .assign(parent_FBgn = lambda df: df.parent_id.map(mapper))
    .drop(['child_id', 'parent_id', 'note', 'm_type'], axis=1)
    .dropna()
) 

movement.loc[(movement.parent_chrom == "chrX") & movement.child_chrom.isin(autosomes), 'movement'] = 'X -> A'
movement.loc[movement.parent_chrom.isin(autosomes) & (movement.parent_chrom != movement.child_chrom), 'movement'] = 'A -> '
movement.movement = pd.Categorical(movement.movement, ordered=True, categories=['X -> A', 'A -> '])

# %%
movement

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
