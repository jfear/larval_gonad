#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate
from IPython.display import display

import os

try:
    os.chdir(os.path.join(os.getcwd(), "docs"))
    print(os.getcwd())
except:
    pass

#%%
muller = pd.read_feather("../output/neox-wf/muller_arm_consensus.feather").set_index("FBgn")

#%%

_df = muller.apply(lambda x: x.value_counts(), axis=0).fillna(0).applymap(lambda x: f'{x:,.0f}')
display(_df)
print(tabulate(_df, headers="keys", tablefmt="github"))


#%%
muller_known = muller[muller.dpse != "unknown_scaffold"].dropna(how='all')

#%%
# Genes conserved across all species
flag_conserved = muller_known.apply(lambda x: len(set(x)) == 1, axis=1)
print(f'Conserved Genes: {flag_conserved.sum():,}')
muller_known[flag_conserved].dpse.value_counts().sort_index()

#%%
# Genes recently conserved between D. mel and D. pse (possibly D. wil)
def recent_conserved(x):
    if (x.dmel == x.dpse) & ((x.dmel != x.dmoj) | (x.dmel != x.dvir)):
        return True
    return False

flag_recent_conserved = muller_known.apply(recent_conserved, axis=1)
print(f'Recent Conserved Genes: {flag_recent_conserved.sum():,}')
muller_known[flag_recent_conserved].dpse.value_counts().sort_index()

#%%
# Gene missing in D. pse but present in other species
def lost(x):
    if  x.isna().dpse & (x.dmel == x.dwil) & (x.dmel == x.dmoj) & (x.dmel == x.dvir):
        return True
    return False

flag_lost = muller_known.apply(lost, axis=1)
print(f'Lost Genes: {flag_lost.sum():,}')
muller_known[flag_lost].dmel.value_counts().sort_index()

#%%
# Genes present in D. pse but on a different chromosome than other species
def moved(x):
    if  ~x.isna().dpse & (x.dmel != x.dpse) & (x.dmel == x.dwil) & (x.dmel == x.dmoj) & (x.dmel == x.dvir):
        return True
    return False

flag_moved = muller_known.apply(moved, axis=1)
print(f'Moved Genes: {flag_moved.sum():,}')
muller_known[flag_moved].dmel.value_counts().sort_index()










#%%
dmel_num_genes_per_arm = muller.loc[:, "dmel"].value_counts().sort_index()

#%% [markdown]
# ## Number of orthologs conserved in all species
#
# First I look at how many orthologs are conserved in all examined species:
# * D. melanogaster
# * D. pseudoobscura
# * D. willistoni
# * D. virilis
# * D. mojavensis

#%%
flag_conserved_in_all = muller.fillna("Z").apply(lambda x: len(set(x)) == 1, axis=1)
num_conserved_in_all = muller.loc[flag_conserved_in_all, "dmel"].value_counts().sort_index()

ax = (
    num_conserved_in_all
    .div(dmel_num_genes_per_arm)
    .fillna(0)
    .sort_index(ascending=False)
    .plot(kind="barh")
)
ax.set(
    ylabel="Muller Arm",
    xlabel="Proportion of Genes",
    title="Proportion Conserved Orthologs\n(All Species)",
)
num_conserved_in_all

#%% [markdown]
# ## Number of orthologs conserved among D. mel and D. pse

#%%
flag_conserved_in_dmel_dpse = muller.dmel == muller.dpse
num_conserved_in_dmel_dpse = muller.loc[flag_conserved_in_dmel_dpse, "dmel"].value_counts().sort_index()

ax = (
    num_conserved_in_dmel_dpse
    .div(dmel_num_genes_per_arm)
    .fillna(0)
    .sort_index(ascending=False)
    .plot(kind="barh")
)
ax.set(
    ylabel="Muller Arm",
    xlabel="Proportion of Genes",
    title="Proportion Conserved Orthologs\n(D. mel vs D. pse)",
)
num_conserved_in_dmel_dpse


#%% [markdown]
# ## Conserved in all species but D. pse

#%%
def only_dpse(x):
    if x.dmel == x.dpse:
        return False
    cols = x.index[x.index != "dpse"].values
    if len(set(x[cols])) == 1:
        return True
    else:
        return False

flag_only_dpse = muller.apply(only_dpse, axis=1)
num_only_dpse = muller.loc[flag_only_dpse, "dmel"].value_counts().sort_index()

ax = (
    num_only_dpse
    .div(dmel_num_genes_per_arm)
    .fillna(0)
    .sort_index(ascending=False)
    .plot(kind="barh")
)
ax.set(
    ylabel="Muller Arm",
    xlabel="Proportion of Genes",
    title="Proportion Only Different in D. pse",
)
num_only_dpse

#%%
muller.loc[flag_only_dpse, 'dpse'].value_counts()

#%%
muller.loc[flag_only_dpse, ['dmel', 'dpse']].dropna().query("dpse == 'D'")


#%%
