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
# # Quick Numbers for Paper

# %% [markdown]
# In today's meeting we are going through the paper looking for holes. I am just taking a quick look and filling some of them.

# %%
import os
import sys
import re
from pathlib import Path
from io import StringIO
from yaml import load

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
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
biomakers = nbconfig.seurat.get_biomarkers('res.0.6')
biomarker_genes = biomakers.index.unique()

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
ptrap.index.unique().shape[0],  ptrap.index.unique().intersection(biomarker_genes).shape[0]

# %%
sorted(ptrap.index.unique().intersection(biomarker_genes).map(nbconfig.fbgn2symbol), key=lambda x: x.lower())

# %%


# %%
lit_genes = {
    'germline': ['vas', 'bam', 'Phf7', 'CG11697', 'p53', 'nos', 'bgcn', 'tut', 'Rbp9', 'peb', 'tej', 'Marf',],
    'late_spermatocytes': ['aly', 'nht', 'soti', 'dj', 'ocn', 'can', 'fzo', 'bol', 'mle', 'mia', 'CG3927', 'sunz', 'sowi', 
                           'd-cup', 'c-cup', 'wa-cup', 'p-cup', 'r-cup', 'oys', 'topi', 'sa', 'CG8368',],
    'cyst': ['tj', 'eya', 'zfh1', 'vn', 'foxo', 'apt', 'ImpL2', 'Wnt4', 'Nrt', 'bnb', 'neur', 'robo2', 'EcR', 'gbb', 'spict', 
             'puc', 'sev', 'hui', 'sano', 'glob1', 'Eip93F', 'fax', 'kek1', 'so',],
    'te': ['nord', 'retn', 'abd-A', 'Abd-B', 'Wnt2', 'Six4', 'CG18628', 'MtnA', 'N',],
    'pc': ['vkg', 'Sox100B', 'bw', 'ems',],
}


gene_annot = []
for k, v in lit_genes.items():
    for gene in v:
        fbgn = nbconfig.symbol2fbgn[gene]
        gene_annot.append((fbgn, gene, k))

lit = pd.DataFrame(gene_annot, columns=['FBgn', 'gene_symbol', 'literature']).set_index('FBgn')

lit.literature.value_counts().to_frame()

# %%
cluster_genes = biomakers.cluster.map(nbconfig.short_cluster_annot).to_frame()

(
    cluster_genes.join(lit, how='right')
    .groupby('literature')
    .cluster.value_counts()
    .to_frame()
    .sort_index()
    .rename(columns={'cluster': 'Genes in BioMarker'})
)

# %%
zscore_max = (
    pd.read_parquet('../output/scrnaseq-wf/tpm_zscore.parquet', columns=nbconfig.cluster_order[:9])
    .rename(columns=dict(zip(nbconfig.cluster_order[:9], nbconfig.short_cluster_order)))
    .idxmax(axis=1)
    .rename('best')
)

# %%
(
    lit.join(zscore_max, how='left')
    .groupby('literature')
    .best.value_counts()
    .rename('Highest Expressed Cluster')
    .rename_axis(['literature', 'cluster'], axis=0)
    .to_frame()
)

# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%


# %%

