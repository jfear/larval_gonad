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
from scipy.cluster.hierarchy import linkage, dendrogram

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
    .pipe(lambda df: df[df.sum(axis=1) > 0])
    .astype(int)
)

# %%
_title = 'All Ptraps With Expression'
_genes = ptrap.index.unique()
_dat = (
    pd.concat([ptrap, expression.reindex(_genes)], keys=['ptrap', 'scRNASeq'], sort=True)
    .rename(index=nbconfig.fbgn2symbol, level=1)
    
    # sort by gene name and if it was a ptrap or rnaseq
    .assign(_type = lambda df: df.index.get_level_values(0).str.lower())
    .assign(lower = lambda df: df.index.get_level_values(1).str.lower())
    .sort_values(['lower', '_type'])
    .drop(['lower', '_type'], axis=1)
    .loc[:, nbconfig.short_cluster_order[:9]]
)

# Brian really wants things sorted by expresison and not name
order = expression.reindex(_genes).sort_values(by=['SP', 'ES', 'MS', 'LS', 'EC', 'MC', 'LC', 'TE', 'PC']).index.map(nbconfig.fbgn2symbol)

stack = []
for gene in order:
    stack.append(_dat.query(f'FBgn == "{gene}"'))

_dat = pd.concat(stack)

# %%
fig, ax = plt.subplots(1, 1, figsize=(8, 30))
sns.heatmap(_dat, yticklabels=True, cmap='viridis', ax=ax, cbar=False)
ax.set_title(_title)
ax.set_ylabel('Gene')

previous = ''
for i, gene in enumerate([l.get_text().replace('ptrap-', '').replace('scRNASeq-', '') for l in ax.get_yticklabels()]):
    if gene == previous:
        continue
    ax.axhline(i, color='w', ls='--')
    previous = gene
    
ax.axvline(4, color='w', ls='--')
ax.axvline(7, color='w', ls='--')
ax.axvline(8, color='w', ls='--')
fig.savefig('../output/notebook/2019-02-05_ptrap_scrnaseq_heatmap.pdf', bbox_inches="tight")

# %%
biomarkers = nbconfig.seurat.get_biomarkers('res.0.6').index.unique().tolist()

# %%
ptrap[~ptrap.index.isin(biomarkers)].index.map(nbconfig.fbgn2symbol).unique().tolist()

# %%



# %%



