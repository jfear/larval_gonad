#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'notebook'))
	print(os.getcwd())
except:
	pass

#%% [markdown]
# # Cell Selection Summary
# In this notebook I am summarizing the cell selection analysis to add to the
# docs.

#%%
SAMPLES = ['testis1', 'testis2', 'testis3', 'testis4']
ORDER = ['cellranger-wf', 'cellranger-force-wf', 'cellranger3-wf', 'droputils']

#%%
from itertools import combinations

import numpy as np
import pandas as pd
from sklearn.metrics import jaccard_similarity_score

from tabulate import tabulate

#%%
def get_calls(sample):
    df = pd.read_feather(f'../output/cellselection-wf/{sample}_combined_cell_calls.feather').set_index('cell_id')
    return df

def jaccard(sample):
    df = get_calls(sample)
    res = []
    for c1, c2 in combinations(df.columns, 2):
        jc = np.round(jaccard_similarity_score(df[c1], df[c2]), 4)
        res.append([c1, c2, jc])
        res.append([c2, c1, jc])

    dfJ = pd.DataFrame(res, columns = ['method 1', 'method 2', 'Jaccard']).set_index(['method 1', 'method 2']).unstack().fillna(1)
    dfJ.columns = dfJ.columns.droplevel(0)
    dfJ = dfJ.loc[ORDER, ORDER]
    return dfJ

#%%
for sample in SAMPLES:
    print(sample)
    print(tabulate(jaccard(sample), headers='keys', tablefmt='github'))
    print('\n\n')

#%% [markdown]
# I am thinking about using a consensus measure. If I include cell ranger v2 with defaults then the consensus will be the same.
#%%
def consensus(sample):
    print(sample)
    df = get_calls(sample)

    # Consensus of all measure: cell ranger defaults, cell ranger force, cell ranger v3, and droplet utils
    flag_all = df.sum(axis=1) == 4

    # Consensus of cell ranger force, cell ranger v3, and droplet utils
    flag_three = df.iloc[:, 1:].sum(axis=1) == 3

    # Consensus of cell ranger v3, and droplet utils
    flag_two = df.iloc[:, 2:].sum(axis=1) == 2

    print('Number of cells with 4-way consensus: ', sum(flag_all))
    print('Number of cells with 3-way consensus: ', sum(flag_three))
    print('Number of cells with 2-way consensus: ', sum(flag_two))
    print('Jaccard cellranger-wf defaults vs full consensus: ', jaccard_similarity_score(df['cellranger-wf'], flag_all))
    print('Jaccard similarity of consensus with 3 vs 2 measures: ', jaccard_similarity_score(flag_three, flag_two))
    print('Number of different calls between consensus with 3 vs 2 measures: ', sum(flag_three != flag_two))
    print('\n\n')

#%%
for sample in SAMPLES:
    consensus(sample)
