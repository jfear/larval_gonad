# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import yaml

# %%
def get_genes_on(cluster: str):
    cell_ids = (
        pd.read_feather("../../output/seurat3-cluster-wf/combined_n3_clusters.feather")
        .query("cluster == @cluster")
        .cell_id.to_list()
    )
    genes_on = (
        pd.read_feather(
            "../../output/cellselection-wf/raw.feather", columns=["FBgn"] + cell_ids
        )
        .set_index("FBgn")
        .pipe(lambda df: df > 0)
        .mean(axis=1)
        .pipe(lambda sr: sr[sr > 0.1])
        .index.to_list()
    )
    return set(genes_on)


# %%
genes_on = {
    "G": get_genes_on("G"),
    "EPS": get_genes_on("EPS"),
    "MPS": get_genes_on("MPS"),
    "LPS": get_genes_on("LPS"),
    "C1": get_genes_on("C1"),
    "C2": get_genes_on("C2"),
    "C3": get_genes_on("C3"),
    "C4": get_genes_on("C4"),
    "T": get_genes_on("T"),
    "P": get_genes_on("P"),
}


# %%
tau = (
    pd.read_feather("../../output/expression-atlas-wf/dmel_tau.feather")
    .set_index("FBgn")
    .male_tau
)
tau


# %%
tau_by_cluster_expressed_genes = (
    pd.DataFrame({k: tau[v].rename(k).dropna() for k, v in genes_on.items()})
    .rename_axis("FBgn")
    .reset_index()
    .melt(id_vars="FBgn", var_name="cluster", value_name="tau")
    .dropna()
)
tau_by_cluster_expressed_genes

# %%
colors = yaml.full_load(open("../../config/colors.yaml"))["clusters"]

# %%
fig = plt.figure(figsize=(8, 4))
sns.violinplot(
    "cluster",
    "tau",
    data=tau_by_cluster_expressed_genes,
    scale="area",
    inner="box",
    width=1,
    palette=colors,
)
plt.suptitle("All Expressed")
plt.savefig("../../output/response-to-review-wf/tau_all_expressed_genes.svg")

# %%
widely_expressed = joblib.load(
    "../../output/cellselection-wf/commonly_expressed_genes.pkl"
)
germline_centric = joblib.load(
    "../../output/response-to-review-wf/germline_centric_subset.pkl"
)


# %%
fig = plt.figure(figsize=(8, 4))
sns.violinplot(
    "cluster",
    "tau",
    data=tau_by_cluster_expressed_genes.query("FBgn in @widely_expressed"),
    scale="area",
    inner="box",
    width=1,
    palette=colors,
)
plt.suptitle("Widely Expressed")
plt.savefig("../../output/response-to-review-wf/tau_widely_expressed_genes.svg")

# %%
fig = plt.figure(figsize=(8, 4))
sns.violinplot(
    "cluster",
    "tau",
    data=tau_by_cluster_expressed_genes.query("FBgn in @germline_centric"),
    scale="area",
    inner="box",
    width=1,
    palette=colors,
)
plt.suptitle("Germline Centric")


# %%
