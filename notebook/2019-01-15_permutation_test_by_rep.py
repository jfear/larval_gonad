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
from larval_gonad.stats import permutation_test_chrom1_lt_chrom2

# %%
# Setup notebook
nbconfig = Nb.setup_notebook(seurat_dir='../output/scrnaseq-wf/scrnaseq_combine_force')

# %%
# read chromosome level counts
r1 = pd.read_csv('../output/scrnaseq-wf/scrnaseq_samples/testis1_force/outs/possorted_genome_bam.bam_counts', 
                 sep='\t', index_col=0)
r1.index = pd.Index([f'rep1_{x}' for x in r1.index], name='cell_id')
r1.chromosome = [f'chr{x}' for x in r1.chromosome.tolist()]

r2 = pd.read_csv('../output/scrnaseq-wf/scrnaseq_samples/testis2_force/outs/possorted_genome_bam.bam_counts', 
                 sep='\t', index_col=0)
r2.index = pd.Index([f'rep2_{x}' for x in r2.index], name='cell_id')
r2.chromosome = [f'chr{x}' for x in r2.chromosome.tolist()]

r3 = pd.read_csv('../output/scrnaseq-wf/scrnaseq_samples/testis3_force/outs/possorted_genome_bam.bam_counts', 
                 sep='\t', index_col=0)
r3.index = pd.Index([f'rep3_{x}' for x in r3.index], name='cell_id')
r3.chromosome = [f'chr{x}' for x in r3.chromosome.tolist()]

reps = pd.concat([r1, r2, r3])

reps_wide = reps.set_index('chromosome', append=True).unstack().fillna(0)
reps_wide.columns = reps_wide.columns.droplevel(0)

# %%
chrom_sizes = pd.read_csv('/data/LCDB/lcdb-references/dmel/r6-16/fasta/dmel_r6-16.chromsizes', sep='\t', header=None, index_col=0)
chrom_sizes.index.name = 'chromosome'
chrom_sizes.columns = ['chrom_size']
chrom_sizes = chrom_sizes.chrom_size
chrom_sizes = chrom_sizes.reindex(nbconfig.chrom_order).copy()

# %%
norm_reps_wide = reps_wide.div(chrom_sizes / 1e7)

# %%
cluster_annot = {
    0: 'LS',
    1: 'MC',
    2: 'MS',
    3: 'ES',
    4: 'LC',
    5: 'EC',
    6: 'SP',
    7: 'TE',
    8: 'PC',
}

cluster_order = [
    'SP',
    'ES',
    'MS',
    'LS',
    'EC',
    'MC',
    'LC',
    'PC',
    'TE',
]

clusters = nbconfig.seurat.get_clusters('res.0.6')
clusters = clusters[clusters < 9].copy()
clusters = clusters.map(cluster_annot)
clusters = clusters.astype('category')
clusters.cat.as_ordered(inplace=True)
clusters.cat.reorder_categories(cluster_order, inplace=True)

# %%
reps_w_clusters = norm_reps_wide.join(clusters, how='right')

# %%
reps_w_clusters['rep'] = reps_w_clusters.index.str.extract('(rep\d)').values

# %%
reps_w_clusters.head()

# %%
reps_w_clusters.rep.value_counts()

# %%


# %%
results = []
for (clus, rep), dd in reps_w_clusters.groupby(['cluster', 'rep']):
    chrom1 = dd['chrX']
    chrom2 = dd[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)
    p_val = permutation_test_chrom1_lt_chrom2(chrom1, chrom2, size=10_000)
    results.append((clus, rep, 'chrX', 'chrA', p_val))

# %%
all_reps = pd.DataFrame(results, columns=['cluster', 'rep', 'chrom1', 'chrom2', 'p_value']).set_index(['cluster', 'rep', 'chrom1', 'chrom2'])

# %%
all_reps

# %%
chroms = ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']
reps_w_clusters.groupby(['cluster', 'rep']).get_group(('MC', 'rep1'))[chroms].plot(kind='kde')

# %%


# %%


# %%


# %%


# %%
r2_wide = r2.set_index('chromosome', append=True).unstack().fillna(0)
r2_wide.columns = r2_wide.columns.droplevel(0)

# %%
r2_w_clusters = r2_wide.div(chrom_sizes / 1e7).join(clusters, how='inner')

# %%
dd = r2_w_clusters.groupby('cluster').get_group('MS')

# %%
(
    permutation_test_chrom1_lt_chrom2(dd['chrX'], dd['chr2L']), 
    permutation_test_chrom1_lt_chrom2(dd['chrX'], dd['chr2R']), 
    permutation_test_chrom1_lt_chrom2(dd['chrX'], dd['chr3L']), 
    permutation_test_chrom1_lt_chrom2(dd['chrX'], dd['chr3R']), 
)

# %%
bob = dd[['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']].copy()
bob['chrA'] = bob[['chr2L', 'chr2R', 'chr3L', 'chr3R']].median(axis=1)

# %%
bob.plot(kind='kde')

# %%
dd

# %%

