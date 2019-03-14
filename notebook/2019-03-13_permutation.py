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

# %%

# %%

# %%
chrom_counts = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
    .join(pd.read_parquet('../output/x-to-a-wf/fbgn2chrom.parquet'))
    .groupby('chrom')
    .sum()
    .T
    .assign(rep=lambda df: df.index.str.extract('(rep\d)', expand=False))
    .join(pd.read_parquet('../output/scrnaseq-wf/clusters.parquet'))
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .dropna()
)

# %%
chrom_counts

# %%
sns.pointplot('cluster', 'chrY', data=chrom_counts)

# %%
sns.kdeplot(chrom_counts.chrY)

# %%

# %%
expressed = (
    pd.read_parquet('../output/x-to-a-wf/expressed_genes_by_chrom.parquet')
    .join(pd.read_parquet('../output/scrnaseq-wf/clusters.parquet'))
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .dropna()
    .assign(rep=lambda df: df.index.str.extract('(rep\d)', expand=False))
)

# %%
expressed

# %%
cnts = (
    pd.read_parquet('../output/x-to-a-wf/raw_by_chrom.parquet')
    .join(pd.read_parquet('../output/scrnaseq-wf/clusters.parquet'))
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .dropna().
)

# %%
from functools import partial

# %%
sns.pointplot('cluster', 'chrY', data=expressed, estimator=np.max, ci=95)

# %%
fig = plt.figure(figsize=plt.figaspect(1/4))
sns.violinplot('cluster', 'chrY', hue='rep', data=expressed, width=1, scale='count')

# %%
fig = plt.figure(figsize=plt.figaspect(1/4))
sns.boxplot('cluster', 'chrY', data=expressed, notch=True)

# %%
fig = plt.figure(figsize=plt.figaspect(1/4))
sns.violinplot('cluster', 'chrX', data=expressed, width=1, scale='count')

# %%
fig = plt.figure(figsize=plt.figaspect(1/4))
sns.violinplot('cluster', 'chr4', data=expressed, width=1, scale='count')

# %%
df = expressed[['chrY', 'cluster']]

# %%
pd.DataFrame(results, columns=['cluster', 'lower', 'median', 'upper'])

# %%
df.groupby('cluster').chrY.apply(lambda x: (x > 0).mean())

# %%
pd.read_pickle('../output/x-to-a-wf/num_genes_by_chrom.pkl')

# %%
(
    pd.read_csv('../output/scrnaseq-wf/scrnaseq_samples/testis1_force/outs/possorted_genome_bam.bam_counts', sep='\t').query('chromosome == "Y"')
    .assign(cell_id=lambda df: 'rep1_' + df.cell_id)
    .set_index('cell_id')
    .join(pd.read_parquet('../output/scrnaseq-wf/clusters.parquet'), how='inner')
    .groupby('cluster')
    .number_reads.median()
)

# %%
cnts.assign(rep=lambda df: df.index.str.extract('(rep\d)', expand=False)).groupby(['rep', 'cluster']).chrY.sum()

# %%

# %%

# %%
x = np.arange(0, 100, 10)
y = np.arange(0, 100, 10)

# %%
from itertools import product

# %%
xy = np.array(list(product(x, y)))

# %%
sizes = [np.random.randint(1, 1000) for i in range(100)]

# %%
plt.scatter(xy[:, 0], xy[:, 1], s=sizes, c=sizes, vmin=0, vmax=1000, cmap='viridis')
plt.colorbar(orientation='horizontal', aspect=2, fraction=.1, label='colors', ticks=[0, 500, 1000])

# %%

# %%
df.head()


# %%

# %%
def prop_cells(x):
    return (x > 0).mean()



# %%
df = (
    pd.read_parquet('../output/scrnaseq-wf/raw.parquet')
    .join(pd.read_parquet('../output/x-to-a-wf/fbgn2chrom.parquet'))
    .query('chrom == "chrY"')
    .drop('chrom', axis=1)
    .assign(gene_symbol=lambda df: df.index.map(nbconfig.fbgn2symbol))
    .set_index('gene_symbol')
    .T
    .join(pd.read_parquet('../output/scrnaseq-wf/clusters.parquet'))
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .rename_axis('cell_id')
    .reset_index()
    .melt(id_vars=['cell_id', 'cluster'], var_name='gene_symbol', value_name='UMI')
    .groupby(['cluster', 'gene_symbol']).agg({'UMI': ['sum', prop_cells]})
)

df.columns = df.columns.droplevel(0)
df = df.reset_index()

df.gene_symbol = pd.Categorical(df.gene_symbol, ordered=True, categories=sorted([x for x in df.gene_symbol.unique() if not x.startswith('Su(Ste)')], key=lambda x: x.lower()))
df.dropna(inplace=True)

# %%
xvals = df.cluster.cat.codes
xlabels = df.cluster.cat.categories

yvals = df.gene_symbol.cat.codes
ylabels = df.gene_symbol.cat.categories

fig, (ax, cax) = plt.subplots(3, 2, figsize=(2, 2), gridspec_kw=dict(width_ratios=[1, .3], height_ratios=[1, .1, 1]))
sc = ax.scatter(xvals, yvals, c = df['sum'], s=df['prop_cells'] * 200, vmin=1, vmax=1000)
ax.set_axisbelow(True)
ax.grid(axis='y', alpha=.2)
plt.colorbar(sc, cax=cax, orientation='horizontal', ticks=[1, 500, 1000], label='Total Number of Reads')

for sizes in [2, 20, 100, 200]:
    plt.scatter([], [], c='k', alpha=0.3, s=sizes, label=f'{sizes / 2:,.0f} %')
plt.legend(scatterpoints=1, frameon=False, labelspacing=.5, title='Percent of Cells', fontsize=7, loc='upper left', bbox_to_anchor=[1, 1])

ax.set(yticks=range(len(ylabels)), xticks=range(len(xlabels)))
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8, fontstyle='italic');
plt.margins(y=0.01)
sns.despine(ax=ax, left=True, bottom=True);

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
