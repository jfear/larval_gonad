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
from statsmodels.api import formula as smf
from tabulate import tabulate

from larval_gonad.io import feather_to_cluster_rep_matrix
from larval_gonad.stats import run_chisq
from larval_gonad.plotting import plot_statsmodels_results

try:
    os.chdir(os.path.join(os.getcwd(), "docs"))
    print(os.getcwd())
except:
    pass


#%%
# Get list of intronless FBgns
fbgns_no_intron = pickle.load(open("../output/science_submission/intron_less_genes.pkl", "rb"))
background = pickle.load(open("../output/science_submission/background_fbgns.pkl", "rb"))

#%%
# Get list of X chromosome genes
fbgn2chrom = (
    pd.read_feather(
        "../references/gene_annotation_dmel_r6-24.feather", columns=["FBgn", "FB_chrom"]
    )
    .set_index("FBgn")
    .squeeze()
)
chrx_fbgns = fbgn2chrom[fbgn2chrom == "X"].index

#%%
# Get gonia biased and cyte biased genes
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
# Munge all into a dataframe
df = bias.copy().join(fbgn2chrom)
df["intronless"] = np.where(df.index.isin(fbgns_no_intron), True, False)
df["X"] = np.where(df.index.isin(chrx_fbgns), True, False)
df["bias"] = "NS"
df.loc[df.gonia_bias, "bias"] = "gonia"
df.loc[df.cyte_bias, "bias"] = "cyte"

#%% [markdown]
# ## How are intronless genes expressed in primary spermatocytes?

#%% [markdown]
# ### Intronless genes are expressed in fewer cells than genes with introns.

#%%
# Plot percent cytes with expression by bias*chrom*intronless
g = sns.FacetGrid(
    df,
    row="bias",
    row_order=["cyte", "gonia", "NS"],
    col="FB_chrom",
    col_order=["X", "2L", "2R", "3L", "3R"],
    sharex=True,
    sharey=True,
    margin_titles=True,
)
g.map(sns.boxplot, "intronless", "pct_cyte", order=[False, True])
g.set_ylabels("% Spermatocyte Cells\nWith Expression")
g.savefig("../output/docs/x_escapers_and_intronless_genes.svg", bbox_inches="tight")

#%% [markdown]
# ### However, intronless genes are enriched in genes with primary spermatocyte biased expression.

#%%
# Cross tab of intronless * bias
ct = pd.crosstab(df.intronless, df.bias)
res = run_chisq(ct).loc[(slice(None), ["observed", "adj std residual", "flag_sig"]), :]
print(tabulate(res.reset_index(), headers="keys", showindex=False, tablefmt="github"))
res

#%%
zscores_intronless = (
    feather_to_cluster_rep_matrix("../output/science_submission/zscore_by_cluster_rep.feather")
    .reindex(fbgns_no_intron)
    .dropna()
)

ax = sns.clustermap(
    zscores_intronless,
    col_cluster=False,
    xticklabels=True,
    yticklabels=False,
    cmap="viridis",
    vmin=-3,
    vmax=3,
    rasterized=True,
)
ax.ax_heatmap.set(xlabel="", ylabel="Intronless Genes")
plt.savefig("../output/docs/x_escapers_and_intronless_genes_heatmap.svg", bbox_inches="tight")


#%% [markdown]
# ## Are intronless genes enriched in X chromosome escapers?

#%% [markdown]
# ### Intronless genes are depleted on the X chromosome.

#%%
# intronless genes across the genome
intronless2chrom = fbgn2chrom.to_frame().query(
    "FB_chrom == ['X', '2L', '2R', '3L', '3R', '4', 'Y']"
)
intronless2chrom["intronless"] = np.where(intronless2chrom.index.isin(fbgns_no_intron), True, False)

ct = pd.crosstab(intronless2chrom.intronless, intronless2chrom.FB_chrom)
res = run_chisq(ct).loc[(slice(None), ["observed", "adj std residual", "flag_sig"]), :]
display(res)

print(tabulate(res.reset_index(), headers="keys", showindex=False, tablefmt="github"))

#%% [markdown]
# ### X chromosome escapers are not enriched for intronless genes.

#%% [markdown]
# #### Main Effects Model Logit(intronless = cyte_biased + X chromosome)

#%%
# Main effects model
model = smf.logit("intronless ~ cyte_bias + X", data=df.replace({True: 1, False: 0}))
results = model.fit()
plot_statsmodels_results(
    "../output/docs/x_escapers_and_intronless_genes_main_effects.png", str(results.summary2())
)
display(results.summary2())

np.exp(results.params).rename("Odds Ratio").to_frame()[results.pvalues <= 0.05]

#%% [markdown]
# #### Full Model Logit(intronless = cyte_biased + X chromosome + cyte_biased * X chromosome)

#%%
# FUll Model
model = smf.logit("intronless ~ cyte_bias * X", data=df.replace({True: 1, False: 0}))
results = model.fit()
plot_statsmodels_results(
    "../output/docs/x_escapers_and_intronless_genes_full.png", str(results.summary2())
)
display(results.summary2())

np.exp(results.params).rename("Odds Ratio").to_frame()[results.pvalues <= 0.05]


#%%
