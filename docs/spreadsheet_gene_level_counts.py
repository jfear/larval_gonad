# -*- coding: utf-8 -*-
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
import yaml

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
config = yaml.load(open('../science_submission/config.yaml'))
headers = yaml.load(open('../science_submission/results_table_config.yaml'))

# %%
writer = pd.ExcelWriter('../output/notebook/spreadsheet_results.xlsx')

# %%
README = writer.book.add_worksheet('README')
README.set_column('A:A', 25)
README.set_column('B:D', 40)


# %%
def iter_list(dat, row=0, col=0):
    for v in dat:
        if isinstance(v, str):
            README.write(row, col, v)
        elif isinstance(v, list):
            row += 1
            col += 1
            row, col = iter_list(v, row, col)
        elif isinstance(v, dict):
            row += 1
            col += 1
            row, col = iter_dict(v, row, col)
        row += 1
    return row, col
            
def iter_dict(dat, row=0, col=0):
    for k, v in dat.items():
        if k == 'fname':
            continue
        README.write(row, col, k)
        col += 1
        if isinstance(v, str):
            README.write(row, col, v)
        elif isinstance(v, list):
            row += 1
            row, col = iter_list(v, row, col)
        elif isinstance(v, dict):
            row += 1
            row, col = iter_dict(v, row, col)
        col -= 1
        row += 1
    return row, col

row, col = iter_dict({'Cluster IDS': config['legend_names']})
row, col = iter_dict(headers, row, col)


# %%
def build_w_rep(fname, title, writer):
    w_rep = (
        pd.read_parquet(fname)
        .assign(gene_symbol=lambda df: df.index.map(nbconfig.fbgn2symbol))
        .assign(chrom=lambda df: df.index.map(nbconfig.fbgn2chrom.iloc[:, 0].to_dict()))
        .set_index(['gene_symbol', 'chrom'], append=True)
        .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
        .assign(rep=lambda df: pd.Categorical(df.rep, ordered=True, categories=['rep1', 'rep2', 'rep3']))
        .sort_values(by=['cluster', 'rep'])
        .set_index(['rep', 'cluster'], append=True)
        .unstack().unstack()
    )

    w_rep.columns = w_rep.columns.droplevel(0)
    w_rep.to_excel(writer, sheet_name=title)


# %%
# Different counts
build_w_rep('../output/scrnaseq-wf/raw_by_cluster_w_rep.parquet', 'Raw Counts', writer)
build_w_rep('../output/scrnaseq-wf/seurat_norm_by_cluster_w_rep.parquet', 'Seurat Normalized Counts', writer)
build_w_rep('../output/scrnaseq-wf/tpm_w_rep.parquet', 'TPM Normalized Counts', writer)
build_w_rep('../output/scrnaseq-wf/rpkm_w_rep.parquet', 'RPKM Normalized Counts', writer)
build_w_rep('../output/scrnaseq-wf/tpm_zscore_w_rep.parquet', 'TPM Zscores', writer)

# %%
# Biomarkers
(
    pd.read_csv('../output/scrnaseq-wf/scrnaseq_combine_force/biomarkers_res.0.6.tsv', sep='\t', index_col=0)
    .rename_axis('FBgn')
    .assign(gene_symbol=lambda df: df.index.map(nbconfig.fbgn2symbol))
    .assign(chrom=lambda df: df.index.map(nbconfig.fbgn2chrom.iloc[:, 0].to_dict()))
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .sort_values(by='cluster')
    .reset_index()
    .melt(id_vars=['FBgn', 'gene_symbol', 'chrom', 'cluster'], var_name='seurat_output', value_name='value')
    .set_index(['FBgn', 'gene_symbol', 'chrom', 'seurat_output', 'cluster'])
    .unstack().unstack()
    .sort_index(level=0)
    .to_excel(writer, sheet_name='Biomarker List')
)

# %%
# tSNE
(
    pd.read_parquet('../output/science_submission/tsne.parquet')
    .to_excel(writer, sheet_name='tSNE')
)

# %%
# X to Autosome Results
(
    pd.read_parquet('../output/x-to-a-wf/expressed_genes_by_chrom.parquet')
    .rename_axis('cell_id')
    .reindex(['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrY'], axis=1)
    .join(
        pd.read_parquet('../output/x-to-a-wf/autosome_ratios_by_cell.parquet')
        .rename_axis('cell_id')
    )
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .assign(rep=lambda df: pd.Categorical(df.rep, ordered=True, categories=['rep1', 'rep2', 'rep3']))
    .set_index(['cluster', 'rep'], append=True)
    .to_excel(writer, sheet_name='Autosome Ratio Results')
)

# %%
# P-values form permutation test
(
    pd.read_parquet('../output/x-to-a-wf/permuted_autosome_ratio_pvalues.parquet')
    .reset_index()
    .assign(cluster=lambda df: pd.Categorical(df.cluster.map(nbconfig.short_cluster_annot), ordered=True, categories=nbconfig.short_cluster_order))
    .set_index('cluster')
    .sort_index()
    .to_excel(writer, 'Autsome Ratio Results')
)

# %%
# SP_vs_Spermatocytes
(
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_cytes.tsv', sep='\t')
    .set_index('primary_FBgn')
    .rename_axis('FBgn')
    .join(nbconfig.fbgn2chrom)
    .set_index(['gene_symbol', 'chrom'], append=True)
    .sort_index(level=0)
    .to_excel(writer, sheet_name='SP vs Spermatocytes')
)

# %%
# SP_vs_early 
(
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/gonia_vs_early.tsv', sep='\t')
    .set_index('primary_FBgn')
    .rename_axis('FBgn')
    .join(nbconfig.fbgn2chrom)
    .set_index(['gene_symbol', 'chrom'], append=True)
    .sort_index(level=0)
    .to_excel(writer, sheet_name='SP vs E1°')
)

# %%
# early_vs_mid 
(
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/early_vs_mid.tsv', sep='\t')
    .set_index('primary_FBgn')
    .rename_axis('FBgn')
    .join(nbconfig.fbgn2chrom)
    .set_index(['gene_symbol', 'chrom'], append=True)
    .sort_index(level=0)
    .to_excel(writer, sheet_name='E1° vs M1°')
)

# %%
# mid_vs_late 
(
    pd.read_csv('../output/scrnaseq-wf/germcell_deg/mid_vs_late.tsv', sep='\t')
    .set_index('primary_FBgn')
    .rename_axis('FBgn')
    .join(nbconfig.fbgn2chrom)
    .set_index(['gene_symbol', 'chrom'], append=True)
    .sort_index(level=0)
    .to_excel(writer, sheet_name='M1° vs L1°')
)

# %%
writer.close()
