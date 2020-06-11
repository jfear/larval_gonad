#%%
import os
import yaml

import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.normalization import tpm as lgTPM
from larval_gonad.normalization import rpkm as lgRPKM
os.chdir('notebook')


#%% [markdown]
# ## Making outputs for Brian Galletta
# 
# Brian asked for an update from the latest iteration for his paper. Here I make all 
# of the different outputs that he has asked for. Including:
#
# * Updated tSNE and UMAP plots (SVGs)
# * Update TPM and RPKM aggregated to cluster level (excel workbook)
# * Differential expression results (excel workbook)

#%%
fbgn2symbol = (
    pd.read_feather('../references/gene_annotation_dmel_r6-26.feather', columns=['FBgn', 'gene_symbol'])
    .set_index('FBgn').squeeze()
)

# Brian's Target Gene List
brians_list = pd.read_csv('../data/external/galletta/trail_list_20190705.txt', sep='\t').FBGN.tolist()

#%% [markdown]
# ## Normalized Read Counts

#%%
## TPM By Cluster
def make_tpm():
    fbgn2length = pd.read_feather('../references/gene_annotation_dmel_r6-26.feather', columns=['FBgn', 'length']).set_index('FBgn').squeeze()

    raw = (
        pd.read_feather('../output/paper_submission/raw_by_cluster_rep.feather')
        .groupby(['FBgn', 'cluster']).sum()
        .squeeze()
        .unstack(level='cluster')
    )

    tpm = lgRPKM(raw, fbgn2length).dropna()
    tpm.index = fbgn2symbol.reindex(tpm.index).to_frame().set_index('gene_symbol', append=True).index
    return tpm

tpm = make_tpm()

## RPKM By Cluster
def make_rpkm():
    fbgn2length = pd.read_feather('../references/gene_annotation_dmel_r6-26.feather', columns=['FBgn', 'length']).set_index('FBgn').squeeze()

    raw = (
        pd.read_feather('../output/paper_submission/raw_by_cluster_rep.feather')
        .groupby(['FBgn', 'cluster']).sum()
        .squeeze()
        .unstack(level='cluster')
    )

    rpkm = lgRPKM(raw, fbgn2length).dropna()
    rpkm.index = fbgn2symbol.reindex(rpkm.index).to_frame().set_index('gene_symbol', append=True).index
    return rpkm

rpkm = make_rpkm()

# Save to excel
with pd.ExcelWriter('../output/notebook/2019-07-05_table_for_galletta.xlsx') as writer:
    # Full TPM
    tpm.to_excel(writer, sheet_name='TPM')

    # Brian's Genes TPM
    tpm.query(f'FBgn == {brians_list}').to_excel(writer, sheet_name='TPM (Genes of Interest)')

    # Full TPM
    rpkm.to_excel(writer, sheet_name='RPKM')

    # Brian's Genes TPM
    rpkm.query(f'FBgn == {brians_list}').to_excel(writer, sheet_name='RPKM (Genes of Interest)')

#%% [markdown]
# ## Differential Expression

#%%
## Gonia vs Primary Spermatocytes
gvc = (
    pd.read_feather('../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather')
    .set_index(['FBgn', 'gene_symbol'])
    # .query('p_val_adj <= 0.01')
    .sort_values('avg_logFC')
)
gvc['direction up regulated'] = 'SP'
gvc.loc[gvc['avg_logFC'] < 0, 'direction up regulated'] = 'EPS|PS1|PS2|PS3'
gvc.loc[gvc['p_val_adj'] > 0.01, 'direction up regulated'] = 'NS'

## Gonia vs Early Primary Spermatocytes
gve = (
    pd.read_feather('../output/seurat3-cluster-wf/combined_n3_gonia_vs_eps.feather')
    .set_index(['FBgn', 'gene_symbol'])
    .sort_values('avg_logFC')
)
gve['direction up regulated'] = 'SP'
gve.loc[gve['avg_logFC'] < 0, 'direction up regulated'] = 'EPS'
gve.loc[gve['p_val_adj'] > 0.01, 'direction up regulated'] = 'NS'

## Early Primary Spermatocytes vs Later Primary Spermatocytes
evp = (
    pd.read_feather('../output/seurat3-cluster-wf/combined_n3_eps_vs_ps.feather')
    .set_index(['FBgn', 'gene_symbol'])
    .sort_values('avg_logFC')
)
evp['direction up regulated'] = 'EPS'
evp.loc[evp['avg_logFC'] < 0, 'direction up regulated'] = 'PS1|PS2|PS3'
evp.loc[evp['p_val_adj'] > 0.01, 'direction up regulated'] = 'NS'

# Save to excel
with pd.ExcelWriter('../output/notebook/2019-07-05_table_for_galletta_diff_expression.xlsx') as writer:
    # G vs Cytes
    gvc.query(f'FBgn == {brians_list}').to_excel(writer, sheet_name='SP vs Primary Spermatocytes')

    # G vs Early Primary Spermatocytes
    gve.query(f'FBgn == {brians_list}').to_excel(writer, sheet_name='S vs EPS')

    # Early Primary Spermatocytes vs Later Primary Spermatocytes
    evp.query(f'FBgn == {brians_list}').to_excel(writer, sheet_name='EPS vs PS1|PS2|PS3')

#%% [markdown]
# ## Plots

#%%
# Get colors scheme
config = yaml.safe_load(open('../config/common.yaml'))
cluster_annot = config['cluster_annot']
cluster_order = config['cluster_order']
colors = dict(zip(cluster_order, yaml.full_load(open('../config/colors.yaml'))['clusters']))

# Get cell id to cluster id
cell_annotation = (
    pd.read_feather('../output/seurat3-cluster-wf/combined_n3_metadata.feather', columns=['cell_id', 'cluster']).set_index('cell_id')
    .assign(cluster=lambda df: df.cluster.cat.rename_categories(cluster_annot))
    .assign(cluster=lambda df: df.cluster.cat.reorder_categories(cluster_order))
    .squeeze()
)

#%% [markdown]
# ### UMAP

#%%
umap = (
    pd.read_feather('../output/seurat3-cluster-wf/combined_n3_umap.feather').set_index('cell_id')
    .join(cell_annotation)
    .assign(color=lambda df: df.cluster.map(colors))
)

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    umap.UMAP_1,
    umap.UMAP_2,
    c=umap.color,
    s=3,
    linewidth=0.02,
    edgecolor="k",
    rasterized=True
)

for clus, row in umap.groupby("cluster").agg({"UMAP_1": "mean", "UMAP_2": "mean"}).iterrows():
    ax.text(
        row.UMAP_1,
        row.UMAP_2,
        clus,
        bbox=dict(facecolor=(1, 1, 1, 0.8), edgecolor="none", pad=0.2),
        ha="center",
        va="center",
        fontweight="bold",
    )

    # clean up plot
    plt.setp(
        ax,
        xlabel="UMAP 1",
        ylabel="UMAP 2",
        aspect="equal",
        xmargin=0,
        ymargin=0,
    )

plt.savefig('../output/notebook/2019-07-05_table_for_galletta_umap.svg', bbox_inches='tight')

#%% [markdown]
# ### tSNE

#%% 
tsne = (
    pd.read_feather('../output/notebook/2019-07-05_output_tsne.feather').set_index('cell_id')
    .join(cell_annotation)
    .assign(color=lambda df: df.cluster.map(colors))
)

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(
    tsne.tSNE_1,
    tsne.tSNE_2,
    c=tsne.color,
    s=3,
    linewidth=0.02,
    edgecolor="k",
    rasterized=True
)

for clus, row in tsne.groupby("cluster").agg({"tSNE_1": "mean", "tSNE_2": "mean"}).iterrows():
    ax.text(
        row.tSNE_1,
        row.tSNE_2,
        clus,
        bbox=dict(facecolor=(1, 1, 1, 0.8), edgecolor="none", pad=0.2),
        ha="center",
        va="center",
        fontweight="bold",
    )

    # clean up plot
    plt.setp(
        ax,
        xlabel="tSNE 1",
        ylabel="tSNE 2",
        aspect="equal",
        xmargin=0,
        ymargin=0,
    )

plt.savefig('../output/notebook/2019-07-05_table_for_galletta_tsne.svg', bbox_inches='tight')
