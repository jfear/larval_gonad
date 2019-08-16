#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import pickle_load, shelve_load
from larval_gonad.config import read_config
from larval_gonad.plotting import format_pval

#%%
import os
try:
    os.chdir(os.path.join(os.getcwd(), 'docs'))
    print(os.getcwd())
except:
    pass

#%%
db = shelve_load("../output/x-to-a-wf/db/expressed.dat", "data", "pvalues")
df = db['data']
pvals = db['pvalues']

#%%
config = read_config("../config/common.yaml")
colors = read_config("../config/colors.yaml")['clusters']

#%% [markdown]
# # Autosome Ratios

#%%
g = sns.FacetGrid(data=df, row="ratio_type", sharey=True, sharex=True, height=4, aspect=1.2)
g.map(sns.boxplot, "cluster", "ratio", showfliers=False, palette=colors, order=config['cluster_order'])

#%% [markdown]
# # Gene Expression

#%%
expressed = pickle_load("../output/cellselection-wf/expressed_genes.pkl")
male = pickle_load("../output/expression-atlas-wf/dmel_male_biased_fbgns.pkl")
female = pickle_load("../output/expression-atlas-wf/dmel_female_biased_fbgns.pkl")

#%%
tpm = (
    pd.read_feather("../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather")
    .set_index("FBgn")
    .assign(flag_expressed=lambda x: x.index.isin(expressed))
    .assign(flag_male=lambda x: x.index.isin(male))
    .assign(flag_female=lambda x: x.index.isin(female))
    .assign(log_tpm=lambda x: np.log10(x.TPM + 1))
    .rename({"log_tpm": "Log10(TPM + 1)"}, axis=1)
)

#%%
defaults = dict(
    x="cluster",
    y="Log10(TPM + 1)",
    showfliers=False,
    notch=True,
    palette=colors,
    order=config['cluster_order'],
)
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(4, 12), sharex=True, sharey=True)
sns.boxplot(data=tpm[tpm.flag_expressed], ax=ax1, **defaults)
sns.boxplot(data=tpm[tpm.flag_male], ax=ax2, **defaults)
sns.boxplot(data=tpm[tpm.flag_female], ax=ax3, **defaults)
ax1.set(xlabel="", title="All Genes")
ax2.set(xlabel="", title="Male-Biased Genes")
ax3.set(xlabel="", title="Female-Biased Genes")

#%% [markdown]
# # Regression

#%%
biomarkers = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_biomarkers.feather")
    .set_index("FBgn")
    .loc[:, ['cluster', 'avg_logFC']]
    .rename({"avg_logFC": "scLFC"}, axis=1)
    .assign(cluster=lambda x: x.cluster.map(config['cluster_annot']))
    .query("cluster == ['EPS', 'PS1', 'PS2', 'PS3']")
)

#%%
fbgn2chrom = pickle_load("../output/x-to-a-wf/fbgn2chrom.pkl")

#%% [markdown]
# # Bulk Adult Male v Female

#%%
m_v_f = (
    pd.read_csv("../output/expression-atlas-wf/dmel_male_biased_expression.tsv", sep="\t", index_col=0)
    .assign(bias="NS")
    .join(fbgn2chrom)
)
m_v_f.loc[(m_v_f.padj <= 0.01) & (m_v_f.log2FoldChange > 0), 'bias'] = "Male"
m_v_f.loc[(m_v_f.padj <= 0.01) & (m_v_f.log2FoldChange < 0), 'bias'] = "Female"

#%%
dat = (
    m_v_f.groupby("bias")
    .chrom.value_counts()
    .unstack()
    .div(fbgn2chrom.chrom.value_counts(), axis=1)
    .pipe(lambda x: x.loc[:, ~x.isnull().all(axis=0)])
    .fillna(0)
    .reset_index()
    .melt(id_vars='bias', var_name='chrom', value_name="Pct Differentially Expressed")

)

#%%
g = sns.FacetGrid(data=dat, col="bias", col_order=['Male', 'Female', 'NS'])
g.map(sns.barplot, "chrom", "Pct Differentially Expressed", order=['X', '2L', '2R', '3L', '3R', '4'])
g.set(ylim=(0, 1))

#%% [markdown]
# # Bulk Adult Testis v Ovary

#%%
m_v_f = (
    pd.read_csv("../output/expression-atlas-wf/dmel_gonad_biased_expression.tsv", sep="\t", index_col=0)
    .assign(bias="NS")
    .join(fbgn2chrom)
)
m_v_f.loc[(m_v_f.padj <= 0.01) & (m_v_f.log2FoldChange > 0), 'bias'] = "Male"
m_v_f.loc[(m_v_f.padj <= 0.01) & (m_v_f.log2FoldChange < 0), 'bias'] = "Female"

#%%
dat = (
    m_v_f.groupby("bias")
    .chrom.value_counts()
    .unstack()
    .div(fbgn2chrom.chrom.value_counts(), axis=1)
    .pipe(lambda x: x.loc[:, ~x.isnull().all(axis=0)])
    .fillna(0)
    .reset_index()
    .melt(id_vars='bias', var_name='chrom', value_name="Pct Differentially Expressed")

)

#%%
g = sns.FacetGrid(data=dat, col="bias", col_order=['Male', 'Female', 'NS'])
g.map(sns.barplot, "chrom", "Pct Differentially Expressed", order=['X', '2L', '2R', '3L', '3R', '4'])
g.set(ylim=(0, 1))

#%%

#%% [markdown]
# # Bulk Larval Testis v Ovary

#%%
m_v_f = (
    pd.read_csv("../output/bulk2-rnaseq-wf/deg/bulk_testis_vs_ovary.tsv", sep="\t", index_col=0)
    .assign(bias="NS")
    .join(fbgn2chrom)
)
m_v_f.loc[(m_v_f.padj <= 0.01) & (m_v_f.log2FoldChange > 0), 'bias'] = "Testis"
m_v_f.loc[(m_v_f.padj <= 0.01) & (m_v_f.log2FoldChange < 0), 'bias'] = "Ovary"

#%%
dat = (
    m_v_f.groupby("bias")
    .chrom.value_counts()
    .unstack()
    .div(fbgn2chrom.chrom.value_counts(), axis=1)
    .pipe(lambda x: x.loc[:, ~x.isnull().all(axis=0)])
    .fillna(0)
    .reset_index()
    .melt(id_vars='bias', var_name='chrom', value_name="Pct Differentially Expressed")

)

#%%
g = sns.FacetGrid(data=dat, col="bias", col_order=['Testis', 'Ovary', 'NS'])
g.map(sns.barplot, "chrom", "Pct Differentially Expressed", order=['X', '2L', '2R', '3L', '3R', '4'])
g.set(ylim=(0, 1))

#%% [markdown]
# # SP vs PS

#%%
sp_v_ps = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_gonia_vs_ps.feather")
    .set_index('FBgn')
    .assign(bias="NS")
    .join(fbgn2chrom)
)
sp_v_ps.loc[(sp_v_ps.p_val_adj <= 0.01) & (sp_v_ps.avg_logFC > 0), 'bias'] = "SP"
sp_v_ps.loc[(sp_v_ps.p_val_adj <= 0.01) & (sp_v_ps.avg_logFC < 0), 'bias'] = "PS"

#%%
dat = (
    sp_v_ps.groupby("bias")
    .chrom.value_counts()
    .unstack()
    .div(fbgn2chrom.chrom.value_counts(), axis=1)
    .pipe(lambda x: x.loc[:, ~x.isnull().all(axis=0)])
    .fillna(0)
    .reset_index()
    .melt(id_vars='bias', var_name='chrom', value_name="Pct Differentially Expressed")

)

#%%
g = sns.FacetGrid(data=dat, col="bias", col_order=['SP', 'PS', 'NS'])
g.map(sns.barplot, "chrom", "Pct Differentially Expressed", order=['X', '2L', '2R', '3L', '3R', '4'])
g.set(ylim=(0, .1))

#%% [markdown]
# # SP vs EPS

#%%
sp_v_ps = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_gonia_vs_eps.feather")
    .set_index('FBgn')
    .assign(bias="NS")
    .join(fbgn2chrom)
)
sp_v_ps.loc[(sp_v_ps.p_val_adj <= 0.01) & (sp_v_ps.avg_logFC > 0), 'bias'] = "SP"
sp_v_ps.loc[(sp_v_ps.p_val_adj <= 0.01) & (sp_v_ps.avg_logFC < 0), 'bias'] = "EPS"

#%%
dat = (
    sp_v_ps.groupby("bias")
    .chrom.value_counts()
    .unstack()
    .div(fbgn2chrom.chrom.value_counts(), axis=1)
    .pipe(lambda x: x.loc[:, ~x.isnull().all(axis=0)])
    .fillna(0)
    .reset_index()
    .melt(id_vars='bias', var_name='chrom', value_name="Pct Differentially Expressed")

)

#%%
g = sns.FacetGrid(data=dat, col="bias", col_order=['SP', 'EPS', 'NS'])
g.map(sns.barplot, "chrom", "Pct Differentially Expressed", order=['X', '2L', '2R', '3L', '3R', '4'])
g.set(ylim=(0, .1))

#%% [markdown]
# # EPS vs PS

#%%
sp_v_ps = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_eps_vs_ps.feather")
    .set_index('FBgn')
    .assign(bias="NS")
    .join(fbgn2chrom)
)
sp_v_ps.loc[(sp_v_ps.p_val_adj <= 0.01) & (sp_v_ps.avg_logFC > 0), 'bias'] = "EPS"
sp_v_ps.loc[(sp_v_ps.p_val_adj <= 0.01) & (sp_v_ps.avg_logFC < 0), 'bias'] = "PS"

#%%
dat = (
    sp_v_ps.groupby("bias")
    .chrom.value_counts()
    .unstack()
    .div(fbgn2chrom.chrom.value_counts(), axis=1)
    .pipe(lambda x: x.loc[:, ~x.isnull().all(axis=0)])
    .fillna(0)
    .reset_index()
    .melt(id_vars='bias', var_name='chrom', value_name="Pct Differentially Expressed")

)

#%%
g = sns.FacetGrid(data=dat, col="bias", col_order=['EPS', 'PS', 'NS'])
g.map(sns.barplot, "chrom", "Pct Differentially Expressed", order=['X', '2L', '2R', '3L', '3R', '4'])
g.set(ylim=(0, .1))



#%%
