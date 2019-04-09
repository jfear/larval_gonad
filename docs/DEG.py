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
# ## General Data

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
fbgn2chrom = pd.read_parquet('../output/x-to-a-wf/fbgn2chrom.parquet').reindex(bg)

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
# ## Pairwise Germline DEG

# %% [markdown] {"toc-hr-collapsed": true}
# ### Are genes coming on during germline development depleted from the X and 4th?

# %% [markdown]
# #### Gonia vs Cytes
#
# * Gonia-biased genes are enriched on the X and 4th
# * Cyte-biased genes are depleted on the X and 4th

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
run_chisq(df)

# %% [markdown]
# #### Gonia vs Mid
#
# * Gonia-biased genes are enriched on the X and 4th
# * M1-biased genes are depleted on the X and 4th

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
run_chisq(df)

# %%
df.sum(axis=1)

# %% [markdown]
# ## Male-biased expression (TCP vs OCP)
#
# Genes with adjusted p-value <= 0.01 and greater than a 2-fold change (i.e., |log2FC| > 1). Roughly 60% of the transcriptome is differentially expressed, with nearly 40% of genes showing male-biased expression.
#
# Genes with male-baised expression are depleted on the X and 4th.

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
df = bulk_sig.join(fbgn2chrom).groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
run_chisq(df)

# %%
df = bulk_sig.join(fbgn2chrom).groupby('chrom').bias.value_counts().unstack()
df = (df.div(df.sum(axis=1), axis='rows') * 100)

fig, ax = plt.subplots(figsize=(2, 4))
df.loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4'], ['testis', 'None', 'ovary']].plot(kind='bar', stacked=True, legend=False, width=.9, color=['C0', 'gray', 'r'], ax=ax)
ax.set_xticklabels([
    l.get_text().replace('chr', '')
    for l in ax.get_xticklabels()
], rotation=0, fontsize=10);

ax.text(0, df.loc['chrX', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')
ax.text(5, df.loc['chr4', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')
ax.margins(0)
sns.despine(ax=ax, left=True)
ax.set_ylabel('% Genes')

# %% [markdown]
# ### Gonia Biased Male expression

# %%
SP_BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC > 0')
    .index.tolist()
)

# %%
df = bulk_sig.reindex(SP_BIASED).join(fbgn2chrom).groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
run_chisq(df)

# %%
df = bulk_sig.reindex(SP_BIASED).join(fbgn2chrom).groupby('chrom').bias.value_counts().unstack()
df = (df.div(df.sum(axis=1), axis='rows') * 100)

# %%
fig, ax = plt.subplots(figsize=(2, 4))
df.loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4'], ['testis', 'None', 'ovary']].plot(kind='bar', stacked=True, legend=False, width=.9, color=['C0', 'gray', 'r'], ax=ax)
ax.set_xticklabels([
    l.get_text().replace('chr', '')
    for l in ax.get_xticklabels()
], rotation=0, fontsize=10);

#ax.text(0, df.loc['chrX', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')
#ax.text(5, df.loc['chr4', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')
ax.margins(0)
sns.despine(ax=ax, left=True)
ax.set_ylabel('% Genes')

# %% [markdown]
# ### Cyte Biased Male expression

# %%
CYTE_BIASED = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

# %%
df = bulk_sig.reindex(CYTE_BIASED).join(fbgn2chrom).groupby('chrom').bias.value_counts().unstack().loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4']].T
run_chisq(df)

# %%
df = bulk_sig.reindex(CYTE_BIASED).join(fbgn2chrom).groupby('chrom').bias.value_counts().unstack()
df = (df.div(df.sum(axis=1), axis='rows') * 100)

# %%
fig, ax = plt.subplots(figsize=(2, 4))
df.loc[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4'], ['testis', 'None', 'ovary']].plot(kind='bar', stacked=True, legend=False, width=.9, color=['C0', 'gray', 'r'], ax=ax)
ax.set_xticklabels([
    l.get_text().replace('chr', '')
    for l in ax.get_xticklabels()
], rotation=0, fontsize=10);

#ax.text(0, df.loc['chrX', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')
#ax.text(5, df.loc['chr4', 'testis'] - 1, '*', color='w', ha='center', va='top', fontsize=12, fontweight='bold')
ax.margins(0)
sns.despine(ax=ax, left=True)
ax.set_ylabel('% Genes')

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

movement = (
    movement
    .merge(bulk_sig.testis_bias.rename('child_testis_bias').to_frame(), left_on='child_FBgn', right_on='FBgn')
    .merge(bulk_sig.testis_bias.rename('parent_testis_bias').to_frame(), left_on='parent_FBgn', right_on='FBgn')
)

# %%
dna = (
    movement.query('gene_type == ["D", "Dl"]')
)


ct = pd.crosstab(dna.movement, dna.parent_testis_bias)
display(ct)
stat, pval = fisher_exact(ct, alternative='two-sided')
print(f"Fisher's exact p-value: {pval:0.4f}\n")

ct = pd.crosstab(dna.movement, dna.child_testis_bias)
display(ct)
stat, pval = fisher_exact(ct, alternative='two-sided')
print(f"Fisher's exact p-value: {pval:0.4f}")

# %%
rna = (
    movement.query('gene_type == ["R", "Rl"]')
)


ct = pd.crosstab(rna.movement, rna.parent_testis_bias)
display(ct)
stat, pval = fisher_exact(ct, alternative='two-sided')
print(f"Fisher's exact p-value: {pval:0.4f}\n")

ct = pd.crosstab(rna.movement, rna.child_testis_bias)
display(ct)
stat, pval = fisher_exact(ct, alternative='two-sided')
print(f"Fisher's exact p-value: {pval:0.4f}")

# %%


# %%
movement.head()

# %%
m1_fbgns = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_mid.tsv', sep='\t', index_col=0)
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index.tolist()
)

# %%
pd.read_parquet('../output/scrnaseq-wf/cpm.parquet').pipe(lambda x: x[x.index.isin(['FBgn0283741'])])

# %%

# %%

# %%
one = pd.read_csv('../output/gene_ts_lengths.tsv', sep='\t', index_col=0).rename(columns={'gene_ts_length': 'one'})
two = pd.read_csv('../output/gene_ts_lengths2.tsv', sep='\t', index_col=0).rename(columns={'gene_ts_length': 'two'})

# %%
df = one.merge(two, on='FBgn', how='outer')

# %%
df[df.isna().any(axis=1)]

# %%
]

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


# %% [markdown] {"toc-hr-collapsed": true}
# ## M1° Gene Expression

# %%
background = (
    pd.read_parquet('../output/scrnaseq-wf/seurat_norm_by_cluster.parquet')
    .join(fbgn2chrom)
    .assign(symbol = lambda df: df.index.map(fbgn2symbol))
    .loc[:, ['symbol', 'chrom']]
)

with open('../output/science_submission/intron_less_genes.pkl', 'rb') as fh:
    intron_less = pickle.load(fh)

background['intron_less'] = False
background.loc[background.index.isin(intron_less), 'intron_less'] = True

with open('../output/science_submission/non_coding_genes.pkl', 'rb') as fh:
    intron_less = pickle.load(fh)

background['non_coding'] = False
background.loc[background.index.isin(intron_less), 'non_coding'] = True

# %%
# intron-less genes are depleted on the X and 3R, while enriched on 2L and 3L
_df = background.groupby('chrom').intron_less.value_counts().unstack(level=0).fillna(0)[config['chrom_order'][:5]]
display(_df)
run_chisq(_df)

# %%
# non-coding genes are depleted on X, 2R, 3R and enriched on 2L, 3L
_df = background.groupby('chrom').non_coding.value_counts().unstack(level=0).fillna(0)[config['chrom_order'][:5]]
display(_df)
run_chisq(_df)

# %%
# intron and non-coding genes are associated
_df = pd.crosstab(background.intron_less, background.non_coding)
display(_df)
run_chisq(_df)

# %%


# %%
pct_cells_expressed = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
    .T
    .join(clusters).dropna()
    .reset_index()
    .melt(id_vars=['cell_id', 'cluster'], var_name='FBgn', value_name='UMI')
    .assign(UMI_on = lambda df: df.UMI > 0)
    .groupby(['cluster', 'FBgn'])
    .UMI_on.mean()
    .unstack(level=0)
    .multiply(100)
)

# %%
print(f"There are {(pct_cells_expressed['M1º'] >= 10).sum():,} expressed in 10% of M1° cells.")

# %%
ax = sns.barplot('cluster', 'num_expressed', data=(pct_cells_expressed >= 10).sum().rename('num_expressed').to_frame().reset_index(), palette=config['colors']['clusters'])
ax.set_ylabel('Number of Genes\nExpressed in >= 10% Cells')

# %%


# %%
m1_coming_on = (
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_mid.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .query('p_val_adj <= 0.01 & avg_logFC < 0')
    .index
)

# %%
print(f"There are {m1_coming_on.shape[0]:,} genes that show M1° biased gene expression (SP vs M1°).")

# %%
background['m1_coming_on'] = False
background.loc[background.index.isin(m1_coming_on), 'm1_coming_on'] = True

# %%
# M1-biased genes are depleted on the X.
df = background.groupby('chrom').m1_coming_on.value_counts().unstack().T[config['chrom_order'][:5]]
run_chisq(df)

# %%
# The 93 of the X-linked escapers are intron-less.
run_chisq(background.query('chrom == "chrX"').groupby('m1_coming_on').intron_less.value_counts().unstack())

# %%
# Only 25 of X-linked escapers are non-coding
run_chisq(background.query('chrom == "chrX"').groupby('m1_coming_on').non_coding.value_counts().unstack())

# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %% [markdown]
# ## scRNA-Seq Expression

# %% [markdown]
# ### Are high expressed genes in the M1° cluster depleted on the X?
#
# There is no association between the number of genes in different expression bins and the chromosome they belong to after correcting for the number of genes on each chromosome (p = 0.4536). 

# %%
# Genes expressed on the X
m1_expression = (
    pd.read_parquet('../output/scrnaseq-wf/raw_by_cluster.parquet')
    .assign(cluster = lambda df: pd.Categorical(df.cluster.map(config['short_cluster_annot']), ordered=True, categories=config['short_cluster_order']))
    .pivot_table(index='FBgn', columns='cluster', values='UMI')
    .loc[:, 'M1º']
)

m1_bins = pd.cut(m1_expression, [-np.inf, 0, 10, 50, 100, 1000, np.inf], labels=['x = 0', 'x < 10', '10 ≤ x < 50', '50 ≤ x < 100', '100 ≤ x < 1,000', '1,000 ≤ x']).rename('bins')
df = pd.concat([m1_bins, fbgn2chrom], join='inner', axis=1).groupby(['bins', 'chrom']).size().unstack().fillna(0)[config['chrom_order'][:5]]
scaled = df.div(num_genes / 1e3, axis='columns').dropna(axis=1)[config['chrom_order'][:5]]
run_chisq(scaled)
