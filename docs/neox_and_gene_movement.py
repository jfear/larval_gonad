#%% [markdown]
# # Exploring gene movement on the NeoX

#%%
# import matplotlib
# matplotlib.use('Agg')
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact, mannwhitneyu
from tabulate import tabulate
from sklearn.metrics import adjusted_rand_score

from larval_gonad.stats import run_chisq


try:
    os.chdir(os.path.join(os.getcwd(), "docs"))
    print(os.getcwd())
except:
    pass

#%%
# Globals
MULLER = "../output/neox-wf/muller_arm_assignment.feather"
ORDER = ["conserved", "moved_on", "gene_death", "moved_off"]
BOXPLOT_DEFAULTS = {"order": ORDER, "notch": True, "showfliers": False}

#%%
# Helper Functions
def _tabulate(df):
    """Output df as markdown."""
    print(
        tabulate(
            df.applymap(lambda x: f"{x:,}").reset_index(),
            headers="keys",
            showindex=False,
            tablefmt="github",
        )
    )
    print()


def box_plot(x, y, data, fname):
    ax = sns.boxplot(x, y, data=df, **BOXPLOT_DEFAULTS)
    sns.swarmplot(x, y, data=df, order=ORDER, color="k", size=4)
    ax.set(
        xlabel="Evolutionary Status\n(D. pseudoobscura)",
        ylabel="Percent Cells with Expression\n(Spermatocytes)",
        title=f"{x.replace('_', ' ').title()}\nSpermatocyte Expression vs Evolutionary Status",
        ylim=(-0.1, 1.1),
    )
    plt.savefig(fname, bbox_inches="tight")
    return ax


def read_deg(fname: str, labels: list, alpha=0.01):
    deg = pd.read_feather(fname).set_index("FBgn").assign(bias="NS")
    deg.loc[(deg.p_val_adj <= alpha) & (deg.avg_logFC > 0), "bias"] = labels[0]
    deg.loc[(deg.p_val_adj <= alpha) & (deg.avg_logFC < 0), "bias"] = labels[1]
    return deg


def read_muller_arm(fname: str):
    return (
        pd.read_feather(fname)
        .set_index("FBgn")
        .squeeze()
        .replace({"recent_conserved": "conserved"})
    )


def combine_deg_and_muller_arms(fname_deg, *muller_fnames):
    deg = read_deg(fname_deg, ["gonia", "cyte"])
    muller_data = pd.concat(list(map(read_muller_arm, muller_fnames)), axis=1, sort=True)
    return deg.join(muller_data, how="inner", sort=True)


def run_mann(df, col):
    pcts = (
        df.assign(
            muller=lambda x: x[col].replace(
                {
                    "conserved": "Selection Favor",
                    "moved_on": "Selection Favor",
                    "gene_death": "Selection Against",
                    "moved_off": "Selection Against",
                    "other": np.nan,
                }
            )
        )
        .pipe(lambda x: x[x.muller.notna()])
        .pivot(columns="muller", values="pct.2")
    )

    _, pval = mannwhitneyu(
        pcts.loc[:, "Selection Favor"], pcts.loc[:, "Selection Against"], alternative="two-sided"
    )

    return pval


def summarize(df, muller_arm: str, species: str, deg: str):
    # Plot percent of cells with cyte expression vs evolutionary status
    fname = f"../output/docs/neox_analysis_boxplot_{deg}_{species}_{muller_arm}.svg"
    box_plot(muller_arm, "pct.2", df, fname)
    plt.show()

    # Fun Mann Whitney U test for distribution medians (Favor vs Against)
    print(f"Mann-Whitney U (Favor vs Against): {run_mann(df, muller_arm)}\n")

    # Cross tabulation 1: all groups
    ct_full = (
        pd.crosstab(df.bias, df[muller_arm].replace({"other": np.nan}))
        .reindex(ORDER, axis=1)
        .fillna(0)
    )
    _tabulate(ct_full)

    # Cross tabulation 2: collapsed groups
    ct_collapsed = pd.crosstab(
        df.bias.replace({"NS": "Not Biased", "gonia": "Not Biased", "cyte": "Cyte Biased"}),
        df[muller_arm].replace(
            {
                "conserved": "Selection Favor",
                "moved_on": "Selection Favor",
                "gene_death": "Selection Against",
                "moved_off": "Selection Against",
                "other": np.nan,
            }
        ),
    )
    _tabulate(ct_collapsed[["Selection Favor", "Selection Against"]])

    # Fisher's Exact Test: collapsed groups
    print(f"Fisher's Exact Test: {fisher_exact(ct_collapsed, alternative='two-sided')[1]}")


#%% [markdown]
# ## Compare Yang et al and FlyBase

#%%
muller_yang = pd.read_feather("../output/neox-wf/YO_muller_arm_assignment.feather").set_index(
    "FBgn"
)
muller_flybase = pd.read_feather("../output/neox-wf/FB_muller_arm_assignment.feather").set_index(
    "FBgn"
)

muller_yang_and_flybase = muller_yang.join(
    muller_flybase, how="outer", lsuffix="_yo", rsuffix="_fb"
)

#%%
# number of orthologs
print(f"Yang Total Number of Orthologs: {muller_yang.shape[0]:,}")
print(f"FlyBase Total Number of Orthologs: {muller_flybase.shape[0]:,}")

#%%
# Similarity excluding missing
def _compare(species):
    _df = muller_yang_and_flybase[[species + "_yo", species + "_fb"]].copy()
    _df.dropna(inplace=True)
    score = adjusted_rand_score(_df.iloc[:, 0], _df.iloc[:, 1])
    return species, _df.shape[0], score


_df = pd.DataFrame(
    list(map(_compare, muller_flybase.columns)),
    columns=["species", "Num Non Missing", "Similarity of Non Missing"],
).set_index("species")
_tabulate(_df)

#%%
# Similarity of missing
def _compare(species):
    _df = muller_yang_and_flybase[[species + "_yo", species + "_fb"]].isnull().copy()
    yo_missing, fb_missing = _df.sum().values
    score = adjusted_rand_score(_df.iloc[:, 0], _df.iloc[:, 1])
    return species, yo_missing, fb_missing, score


_df = pd.DataFrame(
    list(map(_compare, muller_flybase.columns)),
    columns=["species", "Num Yang Missing", "Num FlyBase Missing", "Similarity of Missing"],
).set_index("species")
_tabulate(_df)

################################################################################
#%% [markdown]
# ## Spermatogonia vs Primary Spermatocytes

#%%
gonia_vs_cytes = "../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather"

#%% [markdown]
# ### *D. pseudoobscura*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dpse_muller_A.feather"
muller_d = "../output/neox-wf/dpse_muller_D.feather"
muller_e = "../output/neox-wf/dpse_muller_E.feather"
df = combine_deg_and_muller_arms(gonia_vs_cytes, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dpse", deg="gonia_vs_cytes")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

#%% [markdown]
# ### *D. willistoni*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dwil_muller_A.feather"
muller_d = "../output/neox-wf/dwil_muller_D.feather"
muller_e = "../output/neox-wf/dwil_muller_E.feather"
df = combine_deg_and_muller_arms(gonia_vs_cytes, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dwil", deg="gonia_vs_cytes")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

################################################################################
#%% [markdown]
# ## Spermatogonia vs Early Primary Spermatocytes

#%%
gonia_vs_eps = "../output/seurat3-cluster-wf/combined_n3_gonia_vs_eps.feather"

#%% [markdown]
# ### *D. pseudoobscura*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dpse_muller_A.feather"
muller_d = "../output/neox-wf/dpse_muller_D.feather"
muller_e = "../output/neox-wf/dpse_muller_E.feather"
df = combine_deg_and_muller_arms(gonia_vs_eps, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dpse", deg="gonia_vs_eps")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

#%% [markdown]
# ### *D. willistoni*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dwil_muller_A.feather"
muller_d = "../output/neox-wf/dwil_muller_D.feather"
muller_e = "../output/neox-wf/dwil_muller_E.feather"
df = combine_deg_and_muller_arms(gonia_vs_eps, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dwil", deg="gonia_vs_eps")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

################################################################################
#%% [markdown]
# ## Spermatogonia vs Later Primary Spermatocytes

#%%
gonia_vs_ps = "../output/seurat3-cluster-wf/combined_n3_gonia_vs_ps.feather"

#%% [markdown]
# ### *D. pseudoobscura*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dpse_muller_A.feather"
muller_d = "../output/neox-wf/dpse_muller_D.feather"
muller_e = "../output/neox-wf/dpse_muller_E.feather"
df = combine_deg_and_muller_arms(gonia_vs_ps, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dpse", deg="gonia_vs_ps")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

#%% [markdown]
# ### *D. willistoni*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dwil_muller_A.feather"
muller_d = "../output/neox-wf/dwil_muller_D.feather"
muller_e = "../output/neox-wf/dwil_muller_E.feather"
df = combine_deg_and_muller_arms(gonia_vs_ps, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dwil", deg="gonia_vs_ps")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

################################################################################
#%% [markdown]
# ## Early Primary Spermatocytes vs Later Primary Spermatocytes

#%%
eps_vs_ps = "../output/seurat3-cluster-wf/combined_n3_eps_vs_ps.feather"

#%% [markdown]
# ### *D. pseudoobscura*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dpse_muller_A.feather"
muller_d = "../output/neox-wf/dpse_muller_D.feather"
muller_e = "../output/neox-wf/dpse_muller_E.feather"
df = combine_deg_and_muller_arms(eps_vs_ps, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dpse", deg="eps_vs_ps")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

#%% [markdown]
# ### *D. willistoni*

#%%
# get conservation calls
muller_a = "../output/neox-wf/dwil_muller_A.feather"
muller_d = "../output/neox-wf/dwil_muller_D.feather"
muller_e = "../output/neox-wf/dwil_muller_E.feather"
df = combine_deg_and_muller_arms(eps_vs_ps, muller_a, muller_d, muller_e)

#%%
# Basic counts
_df = (
    pd.concat(
        [df.muller_A.value_counts(), df.muller_D.value_counts(), df.muller_E.value_counts()],
        axis=1,
        sort=True,
    )
    .reindex(["conserved", "recent_conserved", "moved_on", "gene_death", "moved_off", "other"])
    .fillna(0)
)
_tabulate(_df)

#%%
_common = dict(species="dwil", deg="eps_vs_ps")
summarize(df, "muller_A", **_common)
#%%
summarize(df, "muller_D", **_common)
#%%
summarize(df, "muller_E", **_common)

################################################################################
