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

# %%
autosomes = ['chr2L', 'chr2R', 'chr3L', 'chr3R']

# %% [markdown] {"toc-hr-collapsed": true}
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

# %% [markdown] {"toc-hr-collapsed": true}
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
    .assign(child_testis=lambda df: df.child_FBgn.isin(MALE_BIAS))
    .assign(parent_testis=lambda df: df.parent_FBgn.isin(MALE_BIAS))
    .assign(parent_gonia=lambda df: df.parent_FBgn.isin(SP_BIASED))
    .assign(child_gonia=lambda df: df.child_FBgn.isin(SP_BIASED))
    .assign(parent_cyte=lambda df: df.parent_FBgn.isin(CYTE_BIASED))
    .assign(child_cyte=lambda df: df.child_FBgn.isin(CYTE_BIASED))
    .dropna()
) 

movement.loc[(movement.parent_chrom == "chrX") & movement.child_chrom.isin(autosomes), 'movement'] = 'X -> A'
movement.loc[movement.parent_chrom.isin(autosomes) & (movement.parent_chrom != movement.child_chrom), 'movement'] = 'A -> '
movement.movement = pd.Categorical(movement.movement, ordered=True, categories=['X -> A', 'A -> '])

movement.loc[(movement.parent_chrom == "chrX") & movement.child_chrom.isin(autosomes), 'movement2'] = 'X -> A'
movement.loc[movement.parent_chrom.isin(autosomes) & (movement.child_chrom == "chrX"), 'movement2'] = 'A -> X'
movement.loc[movement.parent_chrom.isin(autosomes) & movement.child_chrom.isin(autosomes), 'movement2'] = 'A -> A'
movement.movement2 = pd.Categorical(movement.movement2, ordered=True, categories=['X -> A', 'A -> X', 'A -> A'])

# %%
# Summary counts to get an idea of sample size
pd.concat([
    movement.query('gene_type == ["D", "Dl"]').movement.value_counts().rename("DNA Movement"),
    movement.query('gene_type == ["R", "Rl"]').movement.value_counts().rename("RNA Movement"),
    movement.movement.value_counts().rename("Total"),
], axis=1)

# %%
# Summary counts to get an idea of sample size
pd.concat([
    movement.query('gene_type == ["D", "Dl"]').movement2.value_counts().rename("DNA Movement"),
    movement.query('gene_type == ["R", "Rl"]').movement2.value_counts().rename("RNA Movement"),
    movement.movement2.value_counts().rename("Total"),
], axis=1)

# %% [markdown]
# ### New genes and their parents do not appear to be enriched for testis-biased expression.
#
# #### Hypothesis
#
# MSCI is thought to be caused evolutionary forces and should lead to driving genes important to spermatogenesis off of the X. Here we would expect an enrichment of male-biased genes that moved from the X to an autosome. New genes can arise by both DNA and RNA based methods.
#
# According to [Vibranovski et al 2009](https://paperpile.com/app/p/538fb8a8-11ad-0848-97fe-b837749f3249) parental genes from DNA-based and RNA-based relocations were enriched for 
#
# #### Method
#
# Using either DNR or RNA derrived new genes I performed a fishers exact test for an association of gene movement (X -> A or A ->) for a gene being male biased.
#
# #### Results
#
# There is no significant association between testis-biased expression and gene movement.

# %%
# DNA based relocations
ctp = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement').parent_testis.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement').child_testis.value_counts()
    .unstack()
)

print('DNA based relocations')
print('Fisher Parent: ', fisher_exact(ctp)[1])
print('Fisher Child: ', fisher_exact(ctc)[1])
pd.concat([ctp, ctc], axis=1)

# %%
# DNA based relocations
ctp = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement2').parent_testis.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement2').child_testis.value_counts()
    .unstack()
)

print('DNA based relocations')
resp = run_chisq(ctp)
resc = run_chisq(ctc)
pd.concat([
    resp.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :],
    resc.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
], axis=1)

# %%
# RNA based relocations
ctp = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement').parent_testis.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement').child_testis.value_counts()
    .unstack()
)

print('RNA based relocations')
print('Fisher Parent: ', fisher_exact(ctp)[1])
print('Fisher Child: ', fisher_exact(ctc)[1])

pd.concat([ctp, ctc], axis=1)

# %%
# RNA based relocations
ctp = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement2').parent_testis.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement2').child_testis.value_counts()
    .unstack()
)

print('RNA based relocations')
resp = run_chisq(ctp)
resc = run_chisq(ctc)
pd.concat([
    resp.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :],
    resc.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
], axis=1)

# %%

# %% [markdown]
# ### New genes are enriched in the spermatogonia.
#
# #### Hypothesis
#
# I would expect an enrichment of new genes in spermatogonia.
#
#
# #### Method
#
# Using either DNR or RNA derrived new genes I performed a fishers exact test for an association of gene movement (X -> A or A ->) for a gene being SP-biased.
#
# #### Results
#
# There is an enrichment of child genes whose parents were on the X. This is true for both DNA based and RNA based new genes.

# %%
# DNA based relocations
ctp = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement').parent_gonia.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement').child_gonia.value_counts()
    .unstack()
)

print('DNA based relocations')
print('Fisher Parent: ', fisher_exact(ctp)[1])
print('Fisher Child: ', fisher_exact(ctc)[1])
pd.concat([ctp, ctc], axis=1)

# %%
# DNA based relocations
ctp = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement2').parent_gonia.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement2').child_gonia.value_counts()
    .unstack()
)

print('DNA based relocations')
resp = run_chisq(ctp)
resc = run_chisq(ctc)
pd.concat([
    resp.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :],
    resc.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
], axis=1)

# %%
# RNA based relocations
ctp = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement').parent_gonia.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement').child_gonia.value_counts()
    .unstack()
    .fillna(0)
)

print('RNA based relocations')
print('Fisher Parent: ', fisher_exact(ctp)[1])
print('Fisher Child: ', fisher_exact(ctc)[1])
pd.concat([ctp, ctc], axis=1)

# %%
# RNA based relocations
ctp = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement2').parent_gonia.value_counts()
    .unstack()
    .fillna(0)
)

ctc = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement2').child_gonia.value_counts()
    .unstack()
    .fillna(0)
)

print('RNA based relocations')
resp = run_chisq(ctp)
resc = run_chisq(ctc)
pd.concat([
    resp.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :],
    resc.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
], axis=1)

# %%

# %% [markdown]
# ### New genes are not enriched in spermatocytes.
#
# #### Hypothesis
#
# I would expect an enrichment of new genes in primary spermatocytes.
#
# #### Method
#
# Using either DNR or RNA derrived new genes I performed a fishers exact test for an association of gene movement (X -> A or A ->) for a gene being Cyte-biased.
#
# #### Results
#
# There is no enrichment.

# %%
# DNA based relocations
ctp = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement').parent_cyte.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement').child_cyte.value_counts()
    .unstack()
)

print('DNA based relocations')
print('Fisher Parent: ', fisher_exact(ctp)[1])
print('Fisher Child: ', fisher_exact(ctc)[1])
pd.concat([ctp, ctc], axis=1)

# %%
# DNA based relocations
ctp = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement2').parent_cyte.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["D", "Dl"]')
    .groupby('movement2').child_cyte.value_counts()
    .unstack()
)

print('DNA based relocations')
resp = run_chisq(ctp)
resc = run_chisq(ctc)
pd.concat([
    resp.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :],
    resc.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
], axis=1)

# %%
# RNA based relocations
ctp = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement').parent_cyte.value_counts()
    .unstack()
)

ctc = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement').child_cyte.value_counts()
    .unstack()
    .fillna(0)
)

print('RNA based relocations')
print('Fisher Parent: ', fisher_exact(ctp)[1])
print('Fisher Child: ', fisher_exact(ctc)[1])
pd.concat([ctp, ctc], axis=1)

# %%
# RNA based relocations
ctp = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement2').parent_cyte.value_counts()
    .unstack()
    .fillna(0)
)

ctc = (
    movement.query('gene_type == ["R", "Rl"]')
    .groupby('movement2').child_cyte.value_counts()
    .unstack()
    .fillna(0)
)

print('RNA based relocations')
resp = run_chisq(ctp)
resc = run_chisq(ctc)
pd.concat([
    resp.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :],
    resc.loc[(slice(None), ['observed', 'adj std residual', 'flag_sig']), :]
], axis=1)

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
norm = pd.read_csv('../output/scrnaseq-wf/scrnaseq_combine_force/normalized_read_counts.tsv', sep='\t', index_col=0)

# %%
norm.head()

# %%
movement_pairs = movement[['parent_FBgn', 'child_FBgn']].values

# %%
results = []
for cell_id, dd in norm.T.iterrows():
    for parent, child in movement_pairs:
        results.append([cell_id, parent, child, dd.get(parent, 0.0) - dd.get(child, 0.0)])

# %%
diffs = (
    pd.DataFrame(results, columns=['cell_id', 'parent_FBgn', 'child_FBgn', 'pmc']).set_index('cell_id')
    .join(clusters)#.query('cluster == ["SP", "E1º", "M1º", "L1º"]')
    .reset_index()
    .merge(movement[['parent_FBgn', 'movement2']], on='parent_FBgn')
)

# %%
df = diffs.query('cell_id == "rep1_AAAGTAGTCCGCATAA"')

# %%
sns.boxplot('movement2', 'pmc', data=df.query('pmc != 0'))

# %%

# %%
fig, ax = plt.subplots(figsize=(20, 10))
sns.boxplot('movement2', 'pmc', hue='cluster', data=diffs.query('pmc != 0'), ax=ax)
ax.axhline(0, ls='--', color='k', lw=3)
plt.legend(loc='upper left', bbox_to_anchor=[1, 1])

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
