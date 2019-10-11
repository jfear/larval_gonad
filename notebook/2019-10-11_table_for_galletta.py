#%%
import os
import yaml

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import feather_to_cluster_matrix

try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

#%% [markdown]
# ## Making outputs for Brian Galletta
#
# Brian asked for an update from the latest iteration for his paper. Here I make all
# of the different outputs that he has asked for. Including:
#
# * Updated tSNE and UMAP plots (SVGs)
# * Update TPM and RPKM aggregated to cluster level (excel workbook)
# * Differential expression results (excel workbook)

#%%
fbgn2symbol = (
    pd.read_feather(
        "../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]
    )
    .set_index("FBgn")
    .squeeze()
)

# Brian's Target Gene List
brians_list = pd.read_csv("../data/external/galletta/trial_list.txt", sep="\t").FBGN.tolist()

#%% [markdown]
# ## Normalized Read Counts

#%%
tpm = feather_to_cluster_matrix("../output/seurat3-cluster-wf/tpm_by_cluster.feather")
rpkm = feather_to_cluster_matrix("../output/seurat3-cluster-wf/rpkm_by_cluster.feather")

# Save to excel
with pd.ExcelWriter("../output/notebook/2019-10-11_table_for_galletta.xlsx") as writer:
    # Full TPM
    tpm.to_excel(writer, sheet_name="TPM")

    # Brian's Genes TPM
    tpm.query(f"FBgn == {brians_list}").to_excel(writer, sheet_name="TPM (Genes of Interest)")

    # Full TPM
    rpkm.to_excel(writer, sheet_name="RPKM")

    # Brian's Genes TPM
    rpkm.query(f"FBgn == {brians_list}").to_excel(writer, sheet_name="RPKM (Genes of Interest)")

#%% [markdown]
# ## Differential Expression

#%%
with pd.ExcelWriter("../output/notebook/2019-10-11_table_for_galletta_diff_expression.xlsx") as writer:
    (
        pd.read_csv("../output/plp-wf/trial_list_deg_gonia_vs_spermatocytes.tsv", sep="\t")
        .set_index(["FBgn", "gene_symbol"])
        .sort_values("avg_logFC")
        .assign(
            **{
                "direction up regulated": lambda x: np.where(
                    (x.p_val_adj < 0.01) & (x.avg_logFC > 0),
                    "G",
                    np.where((x.p_val_adj < 0.01) & (x.avg_logFC < 0), "EPS|MPS|LPS", "NS"),
                )
            }
        )
        .to_excel(writer, sheet_name="G vs EPS|MPS|LPS")
    )

    (
        pd.read_csv("../output/plp-wf/trial_list_deg_gonia_vs_eps.tsv", sep="\t")
        .set_index(["FBgn", "gene_symbol"])
        .sort_values("avg_logFC")
        .assign(
            **{
                "direction up regulated": lambda x: np.where(
                    (x.p_val_adj < 0.01) & (x.avg_logFC > 0),
                    "G",
                    np.where((x.p_val_adj < 0.01) & (x.avg_logFC < 0), "EPS", "NS"),
                )
            }
        )
        .to_excel(writer, sheet_name="G vs EPS")
    )

    (
        pd.read_csv("../output/plp-wf/trial_list_deg_gonia_vs_mid_and_late.tsv", sep="\t")
        .set_index(["FBgn", "gene_symbol"])
        .sort_values("avg_logFC")
        .assign(
            **{
                "direction up regulated": lambda x: np.where(
                    (x.p_val_adj < 0.01) & (x.avg_logFC > 0),
                    "G",
                    np.where((x.p_val_adj < 0.01) & (x.avg_logFC < 0), "MPS|LPS", "NS"),
                )
            }
        )
        .to_excel(writer, sheet_name="G vs MPS|LPS")
    )

    (
        pd.read_csv("../output/plp-wf/trial_list_deg_gonia_vs_mid.tsv", sep="\t")
        .set_index(["FBgn", "gene_symbol"])
        .sort_values("avg_logFC")
        .assign(
            **{
                "direction up regulated": lambda x: np.where(
                    (x.p_val_adj < 0.01) & (x.avg_logFC > 0),
                    "G",
                    np.where((x.p_val_adj < 0.01) & (x.avg_logFC < 0), "MPS", "NS"),
                )
            }
        )
        .to_excel(writer, sheet_name="G vs MPS")
    )

    (
        pd.read_csv("../output/plp-wf/trial_list_deg_gonia_vs_late.tsv", sep="\t")
        .set_index(["FBgn", "gene_symbol"])
        .sort_values("avg_logFC")
        .assign(
            **{
                "direction up regulated": lambda x: np.where(
                    (x.p_val_adj < 0.01) & (x.avg_logFC > 0),
                    "G",
                    np.where((x.p_val_adj < 0.01) & (x.avg_logFC < 0), "LPS", "NS"),
                )
            }
        )
        .to_excel(writer, sheet_name="G vs LPS")
    )
