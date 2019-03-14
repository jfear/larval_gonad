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
import pandas as pd

# %%

# %%
pd.read_parquet('../output/scrnaseq-wf/rpkm_zscore.parquet').pivot_table(index='FBgn', columns='cluster', values='rpkm_zscore')

# %%
pd.read_parquet('../output/scrnaseq-wf/rpkm_zscore_w_rep.parquet').pivot_table(index='FBgn', columns=['cluster', 'rep'], values='rpkm_zscore')

# %%

# %%
annotation = nb.short_cluster_annot
cluster_order = nb.short_cluster_order

# %%
biomarkers = (
    pd.read_csv('../output/scrnaseq-wf/scrnaseq_combine_force/biomarkers_res.0.6.tsv', sep='\t', usecols=['primary_FBgn', 'cluster'], index_col=0)
    .cluster
    .map(annotation)
    .pipe(lambda x: x[x != "UNK"])
    .astype('category')
    .cat.as_ordered()
    .cat.reorder_categories(cluster_order)
    .rename_axis('FBgn')
    .to_frame()
)

# %%
unique_fbgns = biomarkers.groupby('FBgn').size().pipe(lambda x: x[x > 1]).index

# %%
df = biomarkers.query(f'FBgn == {unique_fbgns.tolist()}').sort_values(by='cluster')

# %%
multi_clusters = df.groupby('FBgn').apply(lambda df: '|'.join(df.cluster.sort_values().values))

# %%
from itertools import combinations, chain

# %%
germ_combos = ['|'.join(x) for x in chain(
    combinations(cluster_order[:4], 2),
    combinations(cluster_order[:4], 3),
    combinations(cluster_order[:4], 4),
)]

soma_combos = ['|'.join(x) for x in chain(
    combinations(cluster_order[4:], 2),
    combinations(cluster_order[4:], 3),
    combinations(cluster_order[4:], 4),
    combinations(cluster_order[4:], 5),
)]

# %%
flags = pd.Series(index=multi_clusters.index).fillna('GS')
flags[multi_clusters.isin(germ_only)] = 'G'
flags[multi_clusters.isin(soma_only)] = 'S'

# %%
bob = flags.astype('category').cat.as_ordered().cat.reorder_categories(['G', 'GS', 'S']).value_counts().sort_index()

# %%
for i in bob.iteritems():
    break

# %%
flags

# %%
