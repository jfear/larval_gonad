"""Output consensus muller arm calls for Maria.

D. virilis, D. mojavensus, and D. willistoni have their genes mapped to
scaffolds. I used the D. melanogaster annotations to call consensus Muller
arms for these species orthologs. Here I am exporting the table for Maria.

"""
#%%
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'notebook'))
    print(os.getcwd())
except:
    pass

#%%
import pandas as pd

#%%
df = pd.read_feather('../output/neox-wf/muller_arm_consensus.feather')

#%%
df.to_csv('../output/notebook/2019-07-18_consensus_muller_arms_for_maria.tsv', sep='\t', index=False, na_rep='NAN')
