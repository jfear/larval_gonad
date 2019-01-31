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
# # Compare Gal4 to biomarkers

# %%
import os
import sys
import re
from pathlib import Path
from io import StringIO

from IPython.display import display, HTML, Markdown
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
#nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')
nbconfig = Nb.setup_notebook(seurat_dir='/data/fearjm/local_data_store/larval_gonad/output/scrnaseq-wf/scrnaseq_combine_force/')

# %%
biomarkers = nbconfig.seurat.get_biomarkers('res.0.6')
biomarkers.cluster = biomarkers.cluster.map(nbconfig.short_cluster_annot)

# %%
gonia_genes = biomarkers.query('cluster == ["SP", "ES", "MS", "LS"]').index.unique().values.tolist()
cyte_genes = biomarkers.query('cluster == ["EC", "MC", "LC"]').index.unique().values.tolist()
pc_genes = biomarkers.query('cluster == "PC"').index.unique().values.tolist()
te_genes = biomarkers.query('cluster == "TE"').index.unique().values.tolist()

# %%
expression = (pd.read_parquet('../output/scrnaseq-wf/tpm.parquet', columns=nbconfig.cluster_order[:9])
    # change 0s to NaN so I can set them back to 0 after binning values
    .replace(0.0, np.nan)
    .dropna(how='all')
    .apply(pd.cut, bins=3, labels=[1, 2, 3], axis=1)
    .fillna(0)
    .astype(int)
    .rename(columns=dict(zip(nbconfig.cluster_order[:9], nbconfig.short_cluster_order)))
)

# %%
gal4_str="""
bloomington_stocks		gene_symbol	SP	ES	MS	LS	EC	MC	LC	TE	PC	Hub
63356		CG11658	0	0	0	0	1	1	1	0	0	0
65692		Notum	0	0	0	0	2	0	0	0	0	2
49662		CadN	0	0	0	0	0	0	0	0	0	2
62573		Tsp74F	0	0	0	0	2	2	2	0	0	2
63887		sano	0	0	0	0	2	2	2	0	0	0
65690		rau	0	0	0	0	0	1	1	0	0	0
65516		trn	0	0	0	0	0	1	1	0	0	0
62587		Irk1	0	2	2	2	2	2	2	2	0	2
62607		bnl	0	0	0	0	0	1	1	0	0	2
63731		hng3	0	0	0	0	2	2	2	0	0	0
62609		qjt	0	0	0	0	2	0	0	0	0	0
62708		qin	0	0	0	0	1	0	0	0	0	
63304		mael	0	0	0	0	2	0	0	0	0	0
64689		CG31644	0	0	0	0	1	0	1	0	0	0
28849		svp	0	0	0	0	2	2	2	0	0	0
63510		Meltrin	0	0	0	0	0	1	1	0	0	0
63387		Fili	0	0	0	0	0	1	1	0	0	0
63399		Papss	0	0	0	0	0	1	1	0	0	0
62810		AdamTS-A	0	0	0	0	1	0	0	0	0	0
65506		cdi	0	0	0	0	1	1	1	0	0	0
48881		eya	0	0	0	0	1	1	1	0	0	0
"""

gal4 = (
    pd.read_csv(StringIO(gal4_str), sep='\t')
    .dropna(how='all', axis=1)
    .assign(FBgn=lambda df: df.gene_symbol.map(nbconfig.symbol2fbgn))
    .set_index('FBgn')
    .iloc[:, 2:-1]
    .fillna(0)
    .astype(int)
)

# %%
_title = 'Germ Line (Gal4)'
_clusters = nbconfig.short_cluster_order[:4]
_genes = gal4.index.intersection(gonia_genes).intersection(expression.index)
_dat = (
    pd.concat([gal4.reindex(_genes), expression.reindex(_genes)], keys=['gal4', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
)


ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == ["SP", "ES", "MS", "LS"]', engine='python')

# %%
_title = 'Cytes (Gal4)'
_clusters = nbconfig.short_cluster_order[4:7]
_genes = gal4.index.intersection(cyte_genes).intersection(expression.index)
_dat = (
    pd.concat([gal4.reindex(_genes), expression.reindex(_genes)], keys=['gal4', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
)

ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == ["EC", "MC", "LC"]', engine='python')

# %%
_title = 'TE (Gal4)'
_clusters = nbconfig.short_cluster_order[-2]
_genes = gal4.index.intersection(te_genes).intersection(expression.index)
_dat = (
    pd.concat([gal4.reindex(_genes), expression.reindex(_genes)], keys=['gal4', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
    .to_frame()
)

ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == "TE"', engine='python')

# %%
_title = 'PC (Gal4)'
_clusters = nbconfig.short_cluster_order[-1]
_genes = gal4.index.intersection(pc_genes).intersection(expression.index)
_dat = (
    pd.concat([gal4.reindex(_genes), expression.reindex(_genes)], keys=['gal4', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
    .to_frame()
)

ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == "PC"', engine='python')

# %%
t2a_str = """
bloomington_stock_number		gene_symbol	SP	ES	MS	LS	EC	MC	LC	TE	PC	Hub
76191		vari	0	0	0	0	1	1	1	0	0	0
67472		CG42458	0	0	0	0	0	2	2	1	0	0
67449		Wnt4	0	0	0	0	0	2	2	0	0	0
67509		Eaf	1	1	1	1	1	1	1	0	0	0
76159		Zasp52	0	0	0	0	2	1	1	0	0	0
76164		Su(var)2-10	2	2	2	2	0	0	0	0	0	0
76168		GEFmeso	2	2	2	2	2	2	2	2	0	0
77475		ths	0	0	0	0	1	1	1	0	0	0
76157		pk	0	0	0	0	2	2	2	0	0	0
76757		rdo	0	0	0	0	1	0	0	0	2	0
76193		bin3	0	2	2	2	0	0	0	0	0	0
76181		CG2082	0	0	0	0	1	0	0	0	0	0
76222		Khc-73	2	2	2	2	0	0	0	0	0	0
67448		FER	0	0	0	0	2	2	2	0	0	2
76770		CG31075	0	0	0	0	0	1	1	0	0	0
76739		QC	0	0	0	0	2	2	2	0	0	2
66830		dally	0	0	0	0	2	1	1	0	0	2
66856		CG11658	0	0	0	0	2	2	2	0	0	2
76230		PyK	0	2	2	2	0	0	0	0	0	0
76678		PH4alphaEFB	0	0	0	0	2	2	2	0	0	2
66785		Ino80	2	2	2	2	0	0	0	0	0	0
		RasGAP1	0	0	0	0	2	2	2	2	0	2
"""

# %%
t2a = (
    pd.read_csv(StringIO(t2a_str), sep='\t')
    .dropna(how='all', axis=1)
    .assign(FBgn=lambda df: df.gene_symbol.map(nbconfig.symbol2fbgn))
    .set_index('FBgn')
    .iloc[:, 2:-1]
    .fillna(0)
    .astype(int)
)

# %%
# Plot gonia validation
_title = 'Germ Line (T2A)'
_clusters = nbconfig.short_cluster_order[:4]
_genes = t2a.index.intersection(gonia_genes).intersection(expression.index)
_dat = (
    pd.concat([t2a.reindex(_genes), expression.reindex(_genes)], keys=['t2a', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
)


ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == ["SP", "ES", "MS", "LS"]', engine='python')

# %%


# %%
_title = 'Cytes (t2a)'
_clusters = nbconfig.short_cluster_order[4:7]
_genes = t2a.index.intersection(cyte_genes).intersection(expression.index)
_dat = (
    pd.concat([t2a.reindex(_genes), expression.reindex(_genes)], keys=['t2a', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
)

ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == ["EC", "MC", "LC"]', engine='python')

# %%
_title = 'TE (t2a)'
_clusters = nbconfig.short_cluster_order[-2]
_genes = t2a.index.intersection(te_genes).intersection(expression.index)
_dat = (
    pd.concat([t2a.reindex(_genes), expression.reindex(_genes)], keys=['t2a', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
    .to_frame()
)

ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == "TE"', engine='python')

# %%
_title = 'PC (t2a)'
_clusters = nbconfig.short_cluster_order[-1]
_genes = t2a.index.intersection(pc_genes).intersection(expression.index)
_dat = (
    pd.concat([t2a.reindex(_genes), expression.reindex(_genes)], keys=['t2a', 'scRNASeq'])
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
    .to_frame()
)

ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

loc = 2
for i in range(len(_genes)):
    ax.axhline(loc, color='w')
    loc +=2

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == "PC"', engine='python')

# %%
ptrap_str = """gene_symbol	SP	ES	MS	LS	C	EC	MC	LC	PC	TE	H
ADD1	1	2	3	3	1	0	0	0	0	1	0
Ance	1	2	2	2	1	2	2	2	1	1	1
ATP8A	0	0	0	0	0	0	0	0	0	0	0
bol	1	1	2	3	0	0	ND	ND	0	0	0
CadN	0	0	0	0	0	0	0	0	0	0	0
CadN	0	0	0	0	0	0	0	0	0	0	0
CadN	0	0	0	0	0	0	0	0	0	0	0
CadN_Avg	0	0	0	0	0	0	0	0	0	0	0
CG17349	0	0	0	0	0	0	0	0	0	0	0
CG17646	0	0	0	0	1	1	1	1	1	1	0
CG3277 (receptor protein-tyrosine kinase)	0	0	0	0	0	0	0	0	0	0	0
CG8100	0	0	0	0	0	0	0	0	0	0	0
CG9747	0	1	1	1	0	1	1	1	1	1	0
Cht5	0	0	0	0	0	1	1	2	0	0	0
cindr	1	1	1	1	1	1	1	1	2	2	2
cmpy	0	0	0	0	0	1	1	1	0	0	0
Dek	2	2	2	2	1	1	1	1	2	2	2
Dh31-R	0	0	0	0	0	0	0	0	0	0	0
dpr17	1	1	1	0	0	0	0	0	0	0	0
e(y)3	2	2	2	2	1	1	1	1	2	2	1
Eaat2	0	0	0	0	0	0	0	0	0	0	0
Efa6	2	2	1	1	0	0	0	0	1	2	0
Efa6	1	1	1	1	1	2	2	2	1	1	1
Efa6	2	1	1	1	2	2	2	2	1	1	2
Efa6_Avg	1.666666667	1.333333333	1	1	1	1.333333333	1.333333333	1.333333333	1	1.333333333	1
Fas3	1	1	1	1	1	1	1	1	1	3	3
Fas3	0	0	0	0	0	0	0	0	0	2	3
Fas3_Avg	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	2.5	3
fln	1	1	1	1	2	2	2	2	1	1	ND
foxo	2	2	2	2	1	1	1	1	1	2	1
Fs(2)Ket	2	2	2	2	1	1	1	1	1	1	1
grim	0	0	0	0	0	0	0	0	0	0	0
haf	0	0	0	0	0	0	0	0	0	0	0
kkv, CG14668	1	1	1	1	1	1	2	2	1	1	1
klu	0	0	0	0	0	0	0	0	0	0	0
Mapmodulin	2	2	2	1	1	1	1	1	2	2	1
mbl	1	2	2	2	1	0	0	0	1	3	1
Mi-2	1	1	1	1	2	2	2	2	1	0	ND
Su(Tpl)	1	1	1	1	2	2	2	2	1	0	ND
Mipp1	0	0	0	0	0	0	0	0	0	0	0
Mlc2	0	0	0	0	0	0	0	0	0	0	0
NFAT	0	0	0	0	0	0	0	0	0	0	0
nkd	0	0	0	0	0	0	0	0	0	0	0
Nlg3	1	1	1	1	1	2	2	1	ND	1	1
Nlg3	0	0	0	0	0	1	1	1	0	0	0
Nlg3_Avg	0.5	0.5	0.5	0.5	0.5	1.5	1.5	1	0	0.5	0.5
nord	1	1	1	1	0	0	0	0	0	2	0
Np	0	0	0	0	0	0	0	0	0	0	0
Nrg	2	1	1	1	2	2	2	2	2	2	2
osa	0	0	0	0	0	0	0	0	0	0	0
osa	1	1	1	1	2	2	2	2	2	2	ND
osa_Avg	0.5	0.5	0.5	0.5	1	1	1	1	1	1	0
p53	2	2	1	0	0	0	0	0	ND	1	0
Pdcd4	3	3	3	3	1	2	2	2	3	3	1
Pde11	0	0	0	0	0	0	0	0	0	0	0
Piezo	0	0	0	0	0	0	0	0	0	3	0
ppk19	0	0	0	0	0	0	0	0	0	0	0
ppk30	0	0	0	0	0	0	0	0	0	0	0
rdo	1	1	1	1	2	2	2	2	1	1	1
rdo	1	1	1	1	2	2	3	3	1	1	3
rdo	1	1	1	1	2	2	3	3	1	1	3
rdo_Avg	1	1	1	1	2	2	2.666666667	2.666666667	1	1	2.333333333
RunxA	1	1	1	1	1	1	1	1	1	1	1
Sap-r	1	1	1	1	2	3	3	3	2	3	1
sca	0	0	0	0	0	0	0	0	0	0	0
SNF4gamma	1	1	1	1	1	1	1	1	1	1	1
Snmp1	1	1	1	1	1	1	1	1	1	1	1
sosie	1	1	1	1	1	2	2	2	1	1	1
spir	1	1	1	1	0	0	0	0	0	0	0
SRPK	2	2	2	2	0	0	0	0	1	1	1
stai	3	2	2	2	2	2	2	2	1	3	2
stg	0	0	0	0	0	0	0	0	0	0	0
Syn	1	1	1	1	1	1	1	2	1	1	1
Syn	0	0	0	0	0	0	0	0	0	0	0
Syn	1	1	1	1	1	1	1	1	1	1	0
Syn_Avg	0.6666666667	0.6666666667	0.6666666667	0.6666666667	0.6666666667	0.6666666667	0.6666666667	1	0.6666666667	0.6666666667	0.3333333333
Tep2	0	1	1	1	0	2	2	2	2	0	0
tok	1	1	1	1	1	0	0	0	1	2	1
tutl	1	1	1	1	0	0	0	0	1	1	0
twin	1	1	1	1	0	0	0	0	0	0	0
VGlut	0	0	0	0	0	0	0	0	0	0	0
"""


# %%
gene_mapper = nbconfig.symbol2fbgn.copy()

gene_mapper.update(
    {
        'ATP8A': 'FBgn0259221',
        'CG3277 (receptor protein-tyrosine kinase)': 'FBgn0031518',
        'kkv, CG14668': 'FBgn0037320',
        'SNF4gamma': 'FBgn0264357',
    }
)

ptrap = (
    pd.read_csv(StringIO(ptrap_str), sep='\t')
    .query('not gene_symbol.str.contains("Avg")', engine='python')
    .assign(FBgn=lambda df: df.gene_symbol.map(gene_mapper))
    .set_index('FBgn')
    .drop(columns=['gene_symbol', 'H', 'C'])
    .replace('ND', 0)
    .astype(int)
)

# %%
_title = 'Germ Line (Ptrap)'
_clusters = nbconfig.short_cluster_order[:4]
_genes = ptrap.index.intersection(gonia_genes).intersection(expression.index).unique()
_dat = (
    pd.concat([ptrap.query(f'FBgn.isin({_genes.tolist()})', engine='python'), expression.reindex(_genes)], keys=['ptrap', 'scRNASeq'], sort=True)
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
)

fig = plt.figure(figsize=(8, 15))
ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

previous = ''
for i, gene in enumerate([l.get_text().replace('ptrap-', '').replace('scRNASeq-', '') for l in ax.get_yticklabels()]):
    if gene == previous:
        continue
    ax.axhline(i, color='w')
    previous = gene

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == ["SP", "ES", "MS", "LS"]', engine='python')

# %%
# Plot gonia validation
_title = 'Cytes (Ptrap)'
_clusters = nbconfig.short_cluster_order[4:7]
_genes = ptrap.index.intersection(cyte_genes).intersection(expression.index).unique()
_dat = (
    pd.concat([ptrap.query(f'FBgn.isin({_genes.tolist()})', engine='python'), expression.reindex(_genes)], keys=['ptrap', 'scRNASeq'], sort=True)
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
)

fig = plt.figure(figsize=(8, 15))
ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

previous = ''
for i, gene in enumerate([l.get_text().replace('ptrap-', '').replace('scRNASeq-', '') for l in ax.get_yticklabels()]):
    if gene == previous:
        continue
    ax.axhline(i, color='w', lw=2)
    previous = gene

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == ["EC", "MC", "LC"]', engine='python')

# %%
_title = 'TE (Ptrap)'
_clusters = nbconfig.short_cluster_order[-2]
_genes = ptrap.index.intersection(te_genes).intersection(expression.index).unique()
_dat = (
    pd.concat([ptrap.query(f'FBgn.isin({_genes.tolist()})', engine='python'), expression.reindex(_genes)], keys=['ptrap', 'scRNASeq'], sort=True)
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
    .to_frame()
)

fig = plt.figure(figsize=(8, 15))
ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

previous = ''
for i, gene in enumerate([l.get_text().replace('ptrap-', '').replace('scRNASeq-', '') for l in ax.get_yticklabels()]):
    if gene == previous:
        continue
    ax.axhline(i, color='w', lw=2)
    previous = gene

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == "TE"', engine='python')

# %%
_title = 'PC (Ptrap)'
_clusters = nbconfig.short_cluster_order[-1]
_genes = ptrap.index.intersection(pc_genes).intersection(expression.index).unique()
_dat = (
    pd.concat([ptrap.query(f'FBgn.isin({_genes.tolist()})', engine='python'), expression.reindex(_genes)], keys=['ptrap', 'scRNASeq'], sort=True)
    .sort_index(level=1)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    .loc[:, _clusters]
    .to_frame()
)

fig = plt.figure(figsize=(8, 15))
ax = sns.heatmap(_dat, yticklabels=True, cmap='viridis')
ax.set_title(_title)
ax.set_ylabel('Gene')

previous = ''
for i, gene in enumerate([l.get_text().replace('ptrap-', '').replace('scRNASeq-', '') for l in ax.get_yticklabels()]):
    if gene == previous:
        continue
    ax.axhline(i, color='w', lw=2)
    previous = gene

# %%
biomarkers.query(f'FBgn.isin({_genes.tolist()}) and cluster == "PC"', engine='python')

# %%


