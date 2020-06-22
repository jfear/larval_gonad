# %% [markdown]
# # Widely Expressed Gene Ontology
#
# This is a GO analysis of the widely expressed genes. Brian wants to summarize
# this to highly level ribbons like displayed on FlyBase.

# %%
import numpy as np
import pandas as pd
import joblib
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from goatools.godag_plot import plot_gos
from goatools.mapslim import mapslim
from larval_gonad.gene_ontology import run_flyslim

# %% [markdown]
# ## Run GO Slim Analysis
# %%
background = joblib.load("../../output/cellselection-wf/expressed_genes.pkl")
widely_expressed = joblib.load(
    "../../output/cellselection-wf/commonly_expressed_genes.pkl"
)

_, res_obj = run_flyslim(widely_expressed, background, cutoff=0.05, return_obj=True)
res = res_obj.run_study(widely_expressed)
res_mapper = {r.GO: r.p_fdr_bh for r in res}

# %% [markdown]
# ## Molecular Function

# %%
def add_other(df: pd.DataFrame, asset: str):
    other_p = 0
    for r in res:
        if r.NS == asset and r.p_fdr_bh <= 0.05 and r.GO not in df.index:
            other_p = max(-np.log10(r.p_fdr_bh), other_p)
    df.loc["GO:999999", "title"] = "other"
    df.loc["GO:999999", "qval"] = other_p
    return df


# %%
fb_mf = {
    "GO:0003824": "enzyme",
    "GO:0098772": "regulator",
    "GO:0038023": "receptor",
    "GO:0005102": "receptor binding",
    "GO:0005215": "transporter",
    "GO:0005198": "structural molecule",
    "GO:0008092": "cytoskeleton binding",
    "GO:0003677": "DNA binding",
    "GO:0003723": "RNA binding",
    "GO:0140110": "transcription factor",
    "GO:0036094": "small molecule binding",
    "GO:0046872": "metal ion binding",
    "GO:0008289": "lipid binding",
    "GO:0030246": "carbohydrate binding",
}

mf_order = [
    "enzyme",
    "regulator",
    "receptor",
    "receptor binding",
    "transporter",
    "structural molecule",
    "cytoskeleton binding",
    "DNA binding",
    "RNA binding",
    "transcription factor",
    "small molecule binding",
    "metal ion binding",
    "lipid binding",
    "carbohydrate binding",
    "other",
]

mf = (
    pd.Series(fb_mf, name="title")
    .rename_axis("GO")
    .to_frame()
    .assign(qval=lambda sr: -np.log10(sr.index.map(res_mapper)))
    .pipe(add_other, asset="MF")
    .groupby("title")
    .qval.max()
    .to_frame()
    .reindex(mf_order)
)

ax = sns.heatmap(
    mf.T,
    vmin=-np.log10(0.05),
    vmax=-np.log10(1e-10),
    square=True,
    cmap="Blues",
    linewidths=0.2,
    linecolor="gray",
    yticklabels=False,
    xticklabels=True,
    cbar=False,
)  # type: plt.Axes

ax.xaxis.set_ticks_position("top")
plt.xticks(rotation=45, ha="left")
ax.set(xlabel="Molecular Function")


# %% [markdown]
# ## Biological Process

fb_bp = {
    "GO:0008283": "cell cycle/proliferation",
    "GO:0007049": "cell cycle/proliferation",
    "GO:0071840": "cellular organization/biogenesis",
    "GO:0051234": "cellular transport/localization",
    "GO:0033036": "cellular transport/localization",
    "GO:0032502": "development",
    "GO:0000003": "reproduction",
    "GO:0002376": "immune system",
    "GO:0050877": "nervous system process",
    "GO:0007610": "behavior",
    "GO:0050896": "response to stimulus",
    "GO:0023052": "signaling",
    "GO:0010467": "gene expression",
    "GO:0019538": "protein metabolism",
    "GO:0006259": "DNA metabolism",
    "GO:0044281": "small molecule metabolism",
}

bp_order = [
    "cell cycle/proliferation",
    "cellular organization/biogenesis",
    "cellular transport/localization",
    "development",
    "reproduction",
    "immune system",
    "nervous system process",
    "behavior",
    "response to stimulus",
    "signaling",
    "gene expression",
    "protein metabolism",
    "DNA metabolism",
    "small molecule metabolism",
    "other",
]

bp = (
    pd.Series(fb_bp, name="title")
    .rename_axis("GO")
    .to_frame()
    .assign(qval=lambda sr: -np.log10(sr.index.map(res_mapper)))
    .pipe(add_other, asset="BP")
    .groupby("title")
    .qval.max()
    .to_frame()
    .reindex(bp_order)
)

ax = sns.heatmap(
    bp.T,
    vmin=-np.log10(0.05),
    vmax=-np.log10(1e-10),
    square=True,
    cmap="Blues",
    linewidths=0.2,
    linecolor="gray",
    yticklabels=False,
    xticklabels=True,
    cbar=False,
)  # type: plt.Axes

ax.xaxis.set_ticks_position("top")
plt.xticks(rotation=45, ha="left")
ax.set(xlabel="Biological Process")


# %% [markdown]
# ## Cellular Component

fb_cc = {
    "GO:0005576": "extracellular",
    "GO:0005829": "cytosol",
    "GO:0005856": "cytoskeleton",
    "GO:0005739": "mitochondrion",
    "GO:0005634": "nucleus",
    "GO:0005694": "chromosome",
    "GO:0016020": "membrane",
    "GO:0071944": "cell periphery",
    "GO:0042995": "cell projection",
    "GO:0005773": "endomembrane system",
    "GO:0012505": "endomembrane system",
    "GO:0045202": "synapse",
    "GO:0032991": "macromolecular complex",
}

cc_order = [
    "extracellular",
    "cytosol",
    "cytoskeleton",
    "mitochondrion",
    "nucleus",
    "chromosome",
    "membrane",
    "cell periphery",
    "cell projection",
    "endomembrane system",
    "synapse",
    "macromolecular complex",
    "other",
]

cc = (
    pd.Series(fb_cc, name="title")
    .rename_axis("GO")
    .to_frame()
    .assign(qval=lambda sr: -np.log10(sr.index.map(res_mapper)))
    .pipe(add_other, asset="CC")
    .groupby("title")
    .qval.max()
    .to_frame()
    .reindex(cc_order)
)

ax = sns.heatmap(
    cc.T,
    vmin=-np.log10(0.05),
    vmax=-np.log10(1e-10),
    square=True,
    cmap="Blues",
    linewidths=0.2,
    linecolor="gray",
    yticklabels=False,
    xticklabels=True,
    cbar=False,
)  # type: plt.Axes

ax.xaxis.set_ticks_position("top")
plt.xticks(rotation=45, ha="left")
ax.set(xlabel="Cellular Component")


# %%
def tweak_ax(ax: plt.Axes, xlabel: str):
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="left")
    ax.set(xlabel=xlabel)


fig = plt.figure(figsize=(4, 10))
gs = GridSpec(nrows=4, ncols=1, height_ratios=[1, 1, 1, 0.1])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[2, 0])
ax4 = fig.add_subplot(gs[3, 0])

defaults = dict(
    vmin=-np.log10(0.05),
    vmax=-np.log10(1e-11),
    square=True,
    cmap="Blues",
    linewidths=0.2,
    linecolor="gray",
    yticklabels=False,
    xticklabels=True,
    cbar=False,
)

sns.heatmap(mf.T, ax=ax1, **defaults)  # type: plt.Axes
sns.heatmap(bp.T, ax=ax2, **defaults)  # type: plt.Axes

defaults.update(
    {
        "cbar": True,
        "cbar_ax": ax4,
        "cbar_kws": {"orientation": "horizontal", "label": "-log10(q-value)"},
    }
)
sns.heatmap(cc.T, ax=ax3, **defaults)  # type: plt.Axes

tweak_ax(ax1, "Molecular Function")
tweak_ax(ax2, "Biological Process")
tweak_ax(ax3, "Cellular Component")

plt.savefig("../../output/response-to-review-wf/widely_expressed_GO_slim_summary.svg")

# %% [markdown]
# ## Save results table

# %%
pd.DataFrame(
    [(r.GO, r.NS, r.goterm.name, r.p_fdr_bh) for r in res],
    columns=["GO", "Asset", "Title", "FDR Adj P-value"],
).sort_values(["Asset", "GO"]).to_csv(
    "../../output/response-to-review-wf/widely_expressed_GO_slim_summary.tsv",
    sep="\t",
    index=False,
)

