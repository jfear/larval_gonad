# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.3
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

# %%
# Project level imports
from larval_gonad.notebook import Nb

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
early_X_s2_ratio = [0.35349596, 0.27698098, 0.60134619, 0.75800253, 0.6973117, 0.73992864, 0.9781163, 
                    1.19147551, 0.7181409, 1.02505399, 0.56599837, 0.67134575, 0.54940688, 0.23465409, 
                    1.61441454, 0.64448692, 0.99730443, 1.02293311, 0.29088681, 0.23431985, 0.29259291, 
                    0.29150711, 0.75703393, 0.90944075, 0.47230178, 0.94196939, 0.61131165, 0.35558758, 
                    1.13131282]


early_A_s2_ratio = [.76764632, 1.25233922, 0.82549, 0.87801023, 0.88671539, 0.75408014, 0.86050087, 
                    0.85326033, 0.67235028, 2.28267822, 0.90044429, 0.87139316, 0.87360117, 0.88105006, 
                    1.06206806, 1.12041043, 1.25796181, 0.71964458, 1.09037369, 0.71922405, 0.79095761, 
                    0.88383748, 0.57420116, 1.0713959, 0.7785934, 0.81304242, 1.1874923, 1.22352008, 0.79170292, 
                    0.9853998, 0.91072355, 0.92223768, 0.90594832, 0.91839518, 0.91561202, 0.63535517, 0.42068222, 
                    0.36915103, 0.76466503, 1.0684652, 0.95125716, 0.9380952, 1.64870023, 1.07708215, 1.2716295, 
                    2.21354095, 1.11341159, 0.6598262, 1.12776755, 1.25988547, 0.98705969, 1.1859313, 0.71563593, 
                    1.47365408, 1.14475138, 1.86296908, 0.78344929, 1.06810294, 1.05662953]

# %%
df = pd.concat([pd.Series(early_X_s2_ratio, name='early_X_s2_ratio'), pd.Series(early_A_s2_ratio, name='early_A_s2_ratio')], axis=1)

# %%
df = df.melt(var_name='_type', value_name='relative_intensity')

# %%
sns.boxplot('_type', 'relative_intensity', data=df, showbox=False, showfliers=False)
sns.swarmplot('_type', 'relative_intensity', data=df)

# %%
