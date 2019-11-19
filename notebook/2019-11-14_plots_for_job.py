#%%
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from larval_gonad.config import read_config
from larval_gonad.plotting.umap import _add_labels, _cleanup_axes
from larval_gonad.plotting.x_to_a import plot_x_to_a

try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

sns.set_style("white")
plt.rcParams['boxplot.showfliers'] = False

#%%
config = read_config("../config/common.yaml")
color_config = read_config("../config/colors.yaml")
colors = color_config["clusters"]

# %%
umap = pd.read_feather("../output/seurat3-cluster-wf/combined_n3_umap.feather").set_index("cell_id")
clusters = pd.read_feather("../output/seurat3-cluster-wf/combined_n3_clusters.feather").set_index(
    "cell_id"
)
df = umap.join(clusters)


# %%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4, 8))
sns.scatterplot(
    x="UMAP_1",
    y="UMAP_2",
    data=df,
    hue="cluster",
    palette=colors,
    s=2,
    linewidth=0.02,
    edgecolor="k",
    rasterized=True,
    legend=False,
    ax=ax1
)
_add_labels(df, ax1)
sns.despine(ax=ax1)
ax1.set_title("Larval Testis\nscRNASeq", fontsize=18)

plot_x_to_a("../output/x-to-a-wf/db/commonly_expressed.bak", colors, config['cluster_order'], ax=ax2)
sns.despine(ax=ax2)
ax2.set_title("X Inactivation", fontsize=18)

plt.savefig("../output/notebook/2019-11-14_plots_for_job.svg")

# %%
