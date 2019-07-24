#%%
import os
import pandas as pd

#%%
df = pd.read_feather('output/neox-wf/ortholog_annotation.feather')

#%%

df.to_csv('output/notebook/2019-07-05_orthologs_for_maria.tsv', sep='\t', index=False)

