# %%
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

# %%
FBGNS = [
    "FBgn0283451",
    "FBgn0000504",
    "FBgn0010909",
    "FBgn0003415",
    "FBgn0004656",
    "FBgn0262519",
    "FBgn0283657",
    "FBgn0000568",
    "FBgn0263289",
    "FBgn0016917",
    "FBgn0086655",
    "FBgn0026401",
    "FBgn0003371",
    "FBgn0052296",
    "FBgn0053556",
]

# %%
norm = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_normalized.feather")
    .set_index("FBgn")
    .reindex(FBGNS)
    .reset_index()
    .melt(id_vars=["FBgn", "gene_symbol"], var_name="cell_id", value_name="norm")
)

# %%
clusters = pd.read_feather("../output/seurat3-cluster-wf/combined_n3_clusters.feather")

# %%
df = pd.merge(norm, clusters, on="cell_id")

# %%
for fbgn in FBGNS:
    dat = df.query(f"FBgn == '{fbgn}'")
    symbol = dat.gene_symbol.values[0]
    fig = plt.figure(figsize=(10, 8))
    ax = sns.violinplot("cluster", "norm", data=dat, scale="width", inner="box")
    ax.set(ylabel="Normalized Count", xlabel="")
    ax.set_title(symbol, fontstyle="italic", fontsize=24, y=.95)
    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(f"../output/notebook/pradeep/{fbgn}.png")

# %%
