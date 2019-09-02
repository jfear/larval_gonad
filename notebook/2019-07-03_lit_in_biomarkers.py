import yaml
import pandas as pd
from more_itertools import flatten

import os
os.chdir('notebook')

fbgn2symbol = (
    pd.read_feather('../references/gene_annotation_dmel_r6-26.feather', columns=['FBgn', 'gene_symbol'])
    .set_index('FBgn')
    .to_dict()['gene_symbol']
)

config = yaml.safe_load(open('../config/common.yaml'))
CLUSTER_ANNOT = config['cluster_annot']
CLUSTER_ORDER = config['cluster_order']
lit_genes = yaml.safe_load(open('../config/literature_genes.yaml'))

lit_genes_all = list(flatten([
    v
    for k, v in lit_genes.items()
]))


bm = (
    pd.read_feather('../output/seurat3-cluster-wf/combined_n3_biomarkers.feather', columns=['FBgn', 'gene_symbol', 'cluster', 'p_val_adj', 'pct.1'])
    .query('p_val_adj <= 0.01')
    .drop_duplicates(subset='FBgn', keep=False)
    .reset_index(drop=True)
    .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
    .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
    .set_index('FBgn')
    .groupby('cluster')
)


def get_lit(cluster):
    print(cluster)
    df = bm.get_group(cluster)
    return df.reindex(lit_genes_all).dropna()

get_lit('SP')

get_lit('EPS')

get_lit('PS1')

get_lit('PS2')

get_lit('PS3')

get_lit('ECY')

get_lit("CY1")

get_lit("CY2")

get_lit("TE")

get_lit("PC")
