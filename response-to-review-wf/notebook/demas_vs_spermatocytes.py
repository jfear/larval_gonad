# %% [markdown]
# # DEG of Bulk RNA-Seq X vs scRNA-Seq


# %%
import numpy as np
import pandas as pd
from IPython.display import display
import joblib

from larval_gonad.config import read_config
from larval_gonad.stats import run_chisq

# %%
def classify(
    sr: pd.Series,
    qval="padj",
    lfc="log2FoldChange",
    up="testis-biased",
    down="ovary-biased",
):
    if sr[qval] > 0.01:
        return "unbiased"

    if sr[lfc] > 0:
        return up

    return down


# %%
# Map FBgn to Chromosome
CHROM = ["X", "2L", "2R", "3L", "3R", "4", "Y"]
fbgn2chrom = (
    pd.read_feather(
        "../../references/gene_annotation_dmel_r6-26.feather",
        columns=["FBgn", "FB_chrom"],
    )
    .set_index("FBgn")
    .squeeze()
    .rename("chrom")
    .pipe(lambda sr: sr[sr.isin(CHROM)])
)

# %%
# Bulk DEG (testis vs ovary) flags with chromosome
bulk = pd.concat(
    [
        (
            pd.read_csv(
                "../../output/bulk-rnaseq-wf/deg/bulk_testis_vs_ovary.tsv", sep="\t"
            )
            .set_index("FBgn")
            .apply(classify, axis=1)
            .rename("deg")
        ),
        fbgn2chrom,
    ],
    axis=1,
    sort=False,
    join="inner",
)

ct = pd.crosstab(bulk.deg, bulk.chrom)[CHROM]
display(ct)
run_chisq(ct)[CHROM]

# %%
# Single Cell DEG (gonia vs cytes) flags with chromosome
sc = pd.concat(
    [
        (
            pd.read_feather("../../output/seurat3-cluster-wf/germline_deg/GvPS.feather")
            .set_index("FBgn")
            .apply(
                classify,
                qval="p_val_adj",
                lfc="avg_logFC",
                up="gonia-biased",
                down="spermatocyte-biased",
                axis=1,
            )
            .rename("deg")
        ),
        fbgn2chrom,
    ],
    axis=1,
    sort=False,
    join="inner",
)
ct = pd.crosstab(sc.deg, sc.chrom)
display(ct)
run_chisq(ct)[CHROM[:-1]]


# %%
# Intersection of Bulk and Single Cell Flags for the X
combined_flags = (
    pd.merge(
        bulk.reset_index(),
        sc.reset_index(),
        suffixes=["_bulk", "_sc"],
        on=["FBgn", "chrom"],
    )
    .set_index(["FBgn"])
    .query("chrom == 'X'")
)
ct = pd.crosstab(combined_flags.deg_bulk, combined_flags.deg_sc)
display(ct)
run_chisq(ct).loc[:, ["unbiased", "gonia-biased", "spermatocyte-biased"]]


# %%
# Are there differences in bulk vs single cell representation.
ct = pd.DataFrame(
    [[756 - 49, 49], [894 - 68, 68], [297 - 11, 11]],
    columns=["bulk-only", "bulk and single cell"],
    index=["testis-biased", "ovary-biased", "unbiased"],
).T
display(ct)
run_chisq(ct)

# %%
expressed_fbgns = set(joblib.load("../../output/cellselection-wf/expressed_genes.pkl"))
cluster_annot = read_config("../../config/common.yaml")["cluster_annot"]
biomarkers = (
    pd.read_feather("../../output/seurat3-cluster-wf/combined_n3_biomarkers.feather")
    .assign(cluster_name=lambda df: df.cluster.map(cluster_annot))
    .drop(["cluster", "gene_symbol"], axis=1)
    .set_index("FBgn")
    .join(fbgn2chrom)
)
biomarkers

# %%
cyst_cells = ["C1", "C2", "C3", "C4"]
cyst_biased = set(biomarkers.query("cluster_name in @cyst_cells").index.unique())

bulk_with_cyst = bulk.reindex(expressed_fbgns).dropna().query("chrom == 'X'").assign(flag="other")
bulk_with_cyst.loc[bulk_with_cyst.index.isin(cyst_biased), "flag"] = "cyst-biased"

ct = pd.crosstab(bulk_with_cyst.flag, bulk_with_cyst.deg)
display(ct)
run_chisq(ct)

# %%
soma_cells = ["P", "T"]
soma_biased = set(biomarkers.query("cluster_name in @soma_cells").index.unique())

bulk_with_soma = bulk.reindex(expressed_fbgns).dropna().query("chrom == 'X'").assign(flag="other")
bulk_with_soma.loc[bulk_with_soma.index.isin(soma_biased), "flag"] = "soma-biased"

ct = pd.crosstab(bulk_with_soma.flag, bulk_with_soma.deg)
display(ct)
run_chisq(ct)


# %%
