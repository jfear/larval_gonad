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

# %% {"pycharm": {"is_executing": false, "name": "#%%\n"}}
from yaml import load
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

config = load(open('../config/common.yaml').read())
sns.set_context('poster')


# %% {"pycharm": {"is_executing": false, "metadata": false, "name": "#%%\n"}}
clusters = (
    pd.read_parquet('../output/scrnaseq-wf/clusters.parquet')
    .assign(cluster=lambda df: df.cluster.map(config['short_cluster_annot']))
)


# %% {"pycharm": {"is_executing": false, "metadata": false, "name": "#%%\n"}}
def parse(sample, rep_num):
    cnts = (
        pd.read_csv(f'../output/aagag-wf/{sample}.tsv', sep='\t')
        .assign(cell_id=lambda df: f'rep{rep_num}_' + df.cell_id)
        .groupby('cell_id').size()
        .rename('aagag_cnts')
    )
    
    df = (
        pd.concat([cnts, clusters], axis=1, sort=True)
        .dropna()
        .assign(log_aagag_cnts=lambda df: np.log10(df.aagag_cnts))
    )
    fig = plt.figure(figsize=(10, 8))
    ax = sns.boxplot('cluster', 'log_aagag_cnts', data=df, order=config['short_cluster_order'])
    ax.set_title(sample)


# %% {"pycharm": {"is_executing": false, "metadata": false, "name": "#%%\n"}}
parse('testis1', 1)


# %% {"pycharm": {"is_executing": false, "metadata": false, "name": "#%%\n"}}
parse('testis2', 2)


# %% {"pycharm": {"is_executing": false, "metadata": false, "name": "#%%\n"}}
parse('testis3', 3)


# %% {"pycharm": {"metadata": false, "name": "#%%\n"}}



# %%

# %% [markdown]
# ## Translocations

# %% {"pycharm": {"is_executing": false, "metadata": false, "name": "#%%\n"}}
def parse2(sample):
    clus = (
        pd.read_csv(f'../output/translocations-wf/translocation_{sample}_force/clusters.tsv', sep='\t', index_col=0)
        .loc[:, 'res.0.4']
        .rename_axis('cell_id')
        .rename('cluster')
    )
    
    cnts = (
        pd.read_csv(f'../output/aagag-wf/translocations_{sample}.tsv', sep='\t')
        .assign(cell_id=lambda df: f'rep1_' + df.cell_id)
        .groupby('cell_id').size()
        .rename('aagag_cnts')
    )
    
    df = (
        pd.concat([cnts, clus], axis=1, sort=True)
        .dropna()
        .assign(log_aagag_cnts=lambda df: np.log10(df.aagag_cnts))
    )
    
    fig = plt.figure(figsize=(10, 8))
    ax = sns.boxplot('cluster', 'log_aagag_cnts', data=df)
    ax.set_title(sample)


# %%
parse2('stock')

# %%
parse2('BM5')

# %%
parse2('BM21')

# %%
