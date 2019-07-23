#%% [markdown]
# # Look at intron-less gene enrichment in Cyte biased expressed genes.

# This is a quick look at if parimary spermatocyte biased genes are enriched in intronless genes.
# Yes this is what we see.

#%%
import os
import pickle
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, contingency
from IPython.display import display, Markdown
import matplotlib.pyplot as plt
import seaborn as sns

try:
    os.chdir(os.path.join(os.getcwd(), "docs"))
    print(os.getcwd())
except:
    pass


#%%
# Get list of intronless FBgns
fbgns_no_intron = pickle.load(open("../output/science_submission/intron_less_genes.pkl", "rb"))
background = pickle.load(open("../output/science_submission/background_fbgns.pkl", "rb"))

#%% [markdown]

#%%
bias = (
    pd.read_feather("../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather")
    .assign(gonia_bias=lambda x: np.where((x.p_val_adj <= 0.01) & (x.avg_logFC > 0), True, False))
    .assign(pct_gonia=lambda x: x["pct.1"])
    .assign(cyte_bias=lambda x: np.where((x.p_val_adj <= 0.01) & (x.avg_logFC < 0), True, False))
    .assign(pct_cyte=lambda x: x["pct.2"])
    .set_index("FBgn")
    .loc[:, ["gonia_bias", "cyte_bias", "pct_gonia", "pct_cyte"]]
    .reindex(background)
    .dropna()
)

#%%
df = bias.copy()
df["intronless"] = np.where(df.index.isin(fbgns_no_intron), True, False)
df["intronless2"] = np.where(df.index.isin(fbgns_no_intron), "intronless", "has_intron")

#%%
# Gonia Biased Enrichment
ct = pd.crosstab(df.gonia_bias, df.intronless)
display(Markdown("### Observed Frequencies"))
display(ct)

display(Markdown("### Expected Frequencies"))
display(pd.DataFrame(contingency.expected_freq(ct), index=ct.index, columns=ct.columns))

display(Markdown(f"**Fisher's exact test:** {fisher_exact(ct)}"))

#%%
# Cyte Biased Enrichment
ct = pd.crosstab(df.cyte_bias, df.intronless)
display(Markdown("### Observed Frequencies"))
display(ct)

display(Markdown("### Expected Frequencies"))
display(pd.DataFrame(contingency.expected_freq(ct), index=ct.index, columns=ct.columns))

display(Markdown(f"**Fisher's exact test:** {fisher_exact(ct)}"))

#%% [markdown]
# * We see a depletion of intronless genes in Spermatogonia biased genes
# * We see an enrichment of intronless genes in Spermatocyte biased genes

#%%
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=plt.figaspect(1 / 2))
sns.boxplot("intronless2", "pct_gonia", data=df, ax=ax1)
ax1.set(title="Gonia Baised", xlabel="", ylabel="Percent of Cells", ylim=(-0.1, 1.1))

sns.boxplot("intronless2", "pct_cyte", data=df, ax=ax2)
ax2.set(title="Cyte Baised", xlabel="", ylabel="")
fig.savefig("../output/docs/2019-07-23_intronless_analysis.svg", bbox_inches="tight")

#%%
zscores = (
    pd.read_feather("../output/science_submission/zscore_by_cluster_rep.feather")
    .set_index(["FBgn", "cluster", "rep"])
    .unstack(level=[-2, -1])
    .reindex(fbgns_no_intron)
    .dropna()
)

#%%
sns.heatmap(zscores, xticklabels=True, yticklabels=False, cmap="viridis", vmin=-3, vmax=3)

#%%
# TODO add X chromosome analysis.