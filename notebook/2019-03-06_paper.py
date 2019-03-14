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

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
raw = pd.read_parquet('../output/scrnaseq-wf/raw.parquet')

# %%
print(f'There are {raw.shape[1]:,} cells that we analyzed.')

# %%
avg_num_expressed_genes = (raw > 0).sum().mean()
print(f'The average number of expressed genes per cell is {avg_num_expressed_genes:,.2f}')

# %%
avg_umi = raw.sum().mean()
print(f'Average UMIs per cell is {avg_umi:,.2f}')

# %%
