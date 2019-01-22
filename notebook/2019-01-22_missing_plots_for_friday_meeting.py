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

# %%
import os
import sys
import re
from pathlib import Path

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
nbconfig = Nb.setup_notebook(seurat_dir='/data/fearjm/local_data_store/larval_gonad/output/scrnaseq-wf/scrnaseq_combine_force')

# %%
mapper = {
    'chrX': 'X',
    'chr2L': 'A',
    'chr2R': 'A',
    'chr3L': 'A',
    'chr3R': 'A',
    'chrY': 'Y',
    'chr4': '4',
}

fbgn2chrom = pd.read_csv('/data/fearjm/local_data_store/larval_gonad/output/fbgn2chrom.tsv', sep='\t', index_col=0)
fbgn2chrom = fbgn2chrom.chrom.map(mapper)

# %%
chrom_gene_number = fbgn2chrom.value_counts()
chrom_gene_number

# %%
clusters = pd.read_csv('/data/fearjm/local_data_store/larval_gonad/output/scrnaseq-wf/scrnaseq_rep2_force/metadata.tsv', sep='\t', usecols=['res.0.6']).iloc[:, 0]
clusters.index.name = 'cell_id'
clusters.name = 'cluster'
clusters = clusters[clusters < 9].map(nbconfig.short_cluster_annot)

# %%


# %%
raw = pd.read_csv('/data/fearjm/local_data_store/larval_gonad/output/scrnaseq-wf/scrnaseq_rep2_force/raw.tsv', sep='\t')
raw.index.name = 'FBgn'
raw.reset_index(inplace=True)

# %%
melted = raw.melt(id_vars='FBgn', var_name='cell_id', value_name='UMI')

# %%
df = melted.join(fbgn2chrom, on='FBgn').join(clusters, on='cell_id').set_index(['cluster', 'cell_id', 'chrom', 'FBgn'])

# %%
df.sort_index(inplace=True)

# %%
df.head()

# %%
num_missing = (df == 0).groupby(['cluster', 'cell_id', 'chrom']).sum()

# %%
num_missing.div(chrom_gene_number.T, axis='rows', level='chrom')

# %%


# %%

