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
# # First Look at Gal4s

# %% [markdown]
# Sharvani has gone through and scored various gal4 lines and I need to figure out what Gal4s are present in the scRNA-Seq data and in what clusters.

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
# create a fbgn2symbol where symbol has '-RNA' appended
fbgn2symbol_RNA = pd.Series([f'{x}-RNA' for x in nbconfig.fbgn2symbol.values()], index=pd.Index(nbconfig.fbgn2symbol.keys(), name='FBgn'), name= 'gene_symbol')

# %%
# Read in zscores
zscore_by_cluster = pd.read_parquet('/data/fearjm/local_data_store/larval_gonad/output/scrnaseq-wf/tpm_zscore.parquet')[nbconfig.sel_cluster_order]
zscore_by_cluster.columns = ['SP',  'ES',  'MS', 'LS',  'EC',  'MC', 'LC',  'TE',  'PC']
zscore_by_cluster = zscore_by_cluster.join(fbgn2symbol_RNA)

# %%
# Read in Gal4
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

gal4 = pd.read_csv(StringIO(gal4_str), sep='\t', usecols=['gene_symbol', 'SP', 'ES', 'MS', 'LS', 'EC', 'MC', 'LC', 'TE', 'PC'])
gal4.index = pd.Index(gal4.gene_symbol.map(nbconfig.symbol2fbgn).values, name='FBgn')
gal4['gene_symbol'] = [f'{x}-gal4' for x in gal4.gene_symbol]
gal4_data = gal4.set_index('gene_symbol', append=True)

# %%
rna_data = zscore_by_cluster.reindex(gal4.index).set_index('gene_symbol', append=True)
data = pd.concat([rna_data, gal4_data])
data.index = data.index.droplevel(0)

# %%
# Plot
fig = plt.figure(figsize=(10, 15))
sns.heatmap(data.sort_index(), cmap='viridis', xticklabels=True, yticklabels=True, vmin=0, vmax=3, cbar=False)
plt.title('Gal4 Expression')

ax = plt.gca()

loc = 2
for i in range(20):
    ax.axhline(loc, color='w', lw=1)
    loc+=2
    
    
ax.axvline(4, color='w', lw=1)
ax.axvline(7, color='w', lw=1)

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

t2a = pd.read_csv(StringIO(t2a_str), sep='\t', usecols=['gene_symbol', 'SP', 'ES', 'MS', 'LS', 'EC', 'MC', 'LC', 'TE', 'PC'])
t2a.index = pd.Index(t2a.gene_symbol.map(nbconfig.symbol2fbgn).values, name='FBgn')
t2a['gene_symbol'] = [f'{x}-t2a' for x in t2a.gene_symbol]
t2a_data = t2a.set_index('gene_symbol', append=True)

# %%
rna_data = zscore_by_cluster.reindex(t2a.index).set_index('gene_symbol', append=True)
data = pd.concat([rna_data, t2a_data])
data.index = data.index.droplevel(0)

# %%
# Plot
fig = plt.figure(figsize=(10, 15))
sns.heatmap(data.sort_index(), cmap='viridis', xticklabels=True, yticklabels=True, vmin=0, vmax=3, cbar=False)
plt.title('T2A Expression')

ax = plt.gca()

loc = 2
for i in range(20):
    ax.axhline(loc, color='w', lw=1)
    loc+=2
    
    
ax.axvline(4, color='w', lw=1)
ax.axvline(7, color='w', lw=1)

# %%
ben_str = """
stock number		gene_symbol	SP	ES	MS	LS	EC	MC	LC	TE	PC	Hub
1946		CG2187	0	0	0	0	2	0	0	0	0	
1808		CrzR	0	0	0	0	2	1	0	0	0	2
1758		Dh31-R	0	0	0	0	1	0	0	0	0	0
#1719		EH										
#1704		EHR	0	0	0	0	1	1	1	0	0	0
1717		Lgr1	0	0	0	0	1	1	1	10	0	0
1660		Pburs	0	0	0	0	0	1	1	0	0	0
"""

ben = pd.read_csv(StringIO(ben_str), sep='\t', usecols=['gene_symbol', 'SP', 'ES', 'MS', 'LS', 'EC', 'MC', 'LC', 'TE', 'PC'], comment='#')
ben.index = pd.Index(ben.gene_symbol.map(nbconfig.symbol2fbgn).values, name='FBgn')
ben['gene_symbol'] = [f'{x}-t2a' for x in ben.gene_symbol]
ben_data = ben.set_index('gene_symbol', append=True)

# %%
rna_data = zscore_by_cluster.reindex(ben.index).set_index('gene_symbol', append=True).dropna()
data = pd.concat([rna_data, ben_data])
data.index = data.index.droplevel(0)

# %%
# Plot
fig = plt.figure(figsize=(10, 15))
sns.heatmap(data.sort_index(), cmap='viridis', xticklabels=True, yticklabels=True, vmin=0, vmax=3, cbar=False)
plt.title('Ben White Expression')

ax = plt.gca()

loc = 2
for i in range(20):
    ax.axhline(loc, color='w', lw=1)
    loc+=2
    
    
ax.axvline(4, color='w', lw=1)
ax.axvline(7, color='w', lw=1)

# %%



# %%



# %%



# %%



# %%


