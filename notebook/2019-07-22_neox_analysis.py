#%% [markdown]
# # Exploring gene movement on the NeoX in D. psedoobscura

#%%
import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq

from IPython.display import display
import os

try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

#%%
# Globals
GONIA_VS_CYTES = "../output/seurat3-cluster-wf/combined_n3_gonia_vs_cytes.feather"
MULLER = "../output/neox-wf/muller_arm_assignment.feather"
MULLER_A = "../output/neox-wf/muller_A.feather"
MULLER_D = "../output/neox-wf/muller_D.feather"
MULLER_E = "../output/neox-wf/muller_E.feather"

ORDER = ["conserved", "moved_on", "gene_death", "moved_off"]
BOXPLOT_DEFAULTS = {"order": ORDER, "notch": True, "showfliers": False}
#%%
# Helper Functions
def plot(x, y, data, **kwargs):
    ax = sns.boxplot(x, y, data=df, **kwargs)
    sns.swarmplot(x, y, data=df, order=ORDER, color="k", size=4)
    ax.set(
        xlabel="Evolutionary Status\n(D. pseudoobscura)",
        ylabel="Percent Cells with Expression\n(Spermatocytes)",
        title=f"{x.replace('_', ' ').title()}\nSpermatocyte Expression vs Evolutionary Status",
    )
    return ax


def run_mann(df, col):
    selection_pcts = (
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

    return mannwhitneyu(
        selection_pcts.loc[:, "Selection Favor"],
        selection_pcts.loc[:, "Selection Against"],
        alternative="two-sided",
    )


#%%
# Import Spermatogonia vs Spermatocyte Differential Expression
deg = pd.read_feather(GONIA_VS_CYTES).set_index("FBgn").assign(bias="NS")
deg.loc[(deg.p_val_adj <= 0.01) & (deg.avg_logFC > 0), "bias"] = "gonia"
deg.loc[(deg.p_val_adj <= 0.01) & (deg.avg_logFC < 0), "bias"] = "cyte"

#%%
# Merge Muller arm information to DEG.
df = deg.join(
    pd.concat(
        [
            pd.read_feather(MULLER_A)
            .set_index("FBgn")
            .squeeze()
            .replace({"recent_conserved": "conserved"}),
            pd.read_feather(MULLER_D)
            .set_index("FBgn")
            .squeeze()
            .replace({"recent_conserved": "conserved"}),
            pd.read_feather(MULLER_E)
            .set_index("FBgn")
            .squeeze()
            .replace({"recent_conserved": "conserved"}),
        ],
        axis=1,
        sort=True,
    ),
    how="inner",
    sort=True,
)

#%% [markdown]
# ## Muller A

#%%
# Plot percent of cells with cyte expression vs evolutionary status
plot("muller_A", "pct.2", df, **BOXPLOT_DEFAULTS)

#%%
# Cross tabulation 1: all groups
ct_full = pd.crosstab(df.bias, df.muller_A.replace({"other": np.nan}))
ct_full.reindex(ORDER, axis=1).fillna(0)

#%%
# Chi-Square: all groups
run_chisq(ct_full).loc[
    (slice(None), ["observed", "expected", "adj std residual", "flag_sig"]), :
].reindex(ORDER, axis=1).dropna(how="all", axis=1)

#%%
# Cross tabulation 2: collapsed groups
ct_collapsed = pd.crosstab(
    df.bias.replace({"NS": "Not Biased", "gonia": "Not Biased", "cyte": "Cyte Biased"}),
    df.muller_A.replace(
        {
            "conserved": "Selection Favor",
            "moved_on": "Selection Favor",
            "gene_death": "Selection Against",
            "moved_off": "Selection Against",
            "other": np.nan,
        }
    ),
)
ct_collapsed[["Selection Favor", "Selection Against"]]

#%%
# Chi-Square: collapsed groups
run_chisq(ct_collapsed)[["Selection Favor", "Selection Against"]]

#%%
# Fun Mann Whitney U test for distribution medians (Favor vs Against)
run_mann(df, "muller_A")

#%% [markdown]
# ## Muller D

#%%
# Plot percent of cells with cyte expression vs evolutionary status
plot("muller_D", "pct.2", df, **BOXPLOT_DEFAULTS)

#%%
# Cross tabulation 1: all groups
ct_full = pd.crosstab(df.bias, df.muller_D.replace({"other": np.nan}))
ct_full.reindex(ORDER, axis=1).fillna(0)

#%%
# Chi-Square: all groups
run_chisq(ct_full).loc[
    (slice(None), ["observed", "expected", "adj std residual", "flag_sig"]), :
].reindex(ORDER, axis=1).dropna(how="all", axis=1)

#%%
# Cross tabulation 2: collapsed groups
ct_collapsed = pd.crosstab(
    df.bias.replace({"NS": "Not Biased", "gonia": "Not Biased", "cyte": "Cyte Biased"}),
    df.muller_D.replace(
        {
            "conserved": "Selection Favor",
            "moved_on": "Selection Favor",
            "gene_death": "Selection Against",
            "moved_off": "Selection Against",
            "other": np.nan,
        }
    ),
)
ct_collapsed[["Selection Favor", "Selection Against"]]

#%%
# Chi-Square: collapsed groups
run_chisq(ct_collapsed)[["Selection Favor", "Selection Against"]]

#%%
# Fun Mann Whitney U test for distribution medians (Favor vs Against)
run_mann(df, "muller_D")

#%% [markdown]
# ## Muller E

#%%
# Plot percent of cells with cyte expression vs evolutionary status
plot("muller_E", "pct.2", df, **BOXPLOT_DEFAULTS)

#%%
# Cross tabulation 1: all groups
ct_full = pd.crosstab(df.bias, df.muller_E.replace({"other": np.nan}))
ct_full.reindex(ORDER, axis=1).fillna(0)

#%%
# Chi-Square: all groups
run_chisq(ct_full).loc[
    (slice(None), ["observed", "expected", "adj std residual", "flag_sig"]), :
].reindex(ORDER, axis=1).dropna(how="all", axis=1)

#%%
# Cross tabulation 2: collapsed groups
ct_collapsed = pd.crosstab(
    df.bias.replace({"NS": "Not Biased", "gonia": "Not Biased", "cyte": "Cyte Biased"}),
    df.muller_E.replace(
        {
            "conserved": "Selection Favor",
            "moved_on": "Selection Favor",
            "gene_death": "Selection Against",
            "moved_off": "Selection Against",
            "other": np.nan,
        }
    ),
)
ct_collapsed[["Selection Favor", "Selection Against"]]

#%%
# Chi-Square: collapsed groups
run_chisq(ct_collapsed)[["Selection Favor", "Selection Against"]]

#%%
# Fun Mann Whitney U test for distribution medians (Favor vs Against)
run_mann(df, "muller_D")