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
        pd.concat([cnts, clusters], axis=1)
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


