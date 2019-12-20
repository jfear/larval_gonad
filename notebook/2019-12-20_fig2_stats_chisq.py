# %%
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.stats import run_chisq

# %%
try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

# %% [markdown]
# # Adult Bulk

# %%
adult_ct = (
    pd.read_feather("../output/expression-atlas-wf/w1118_gene_counts.feather")
    .assign(flag_on=lambda x: x.Count >= 5)
    .groupby(["tissue", "chrom"])
    .flag_on.sum()
    .unstack()
    .reindex(index=["testis", "ovary"], columns="X,2L,2R,3L,3R,4,Y".split(","))
)  # type: pd.DataFrame

run_chisq(adult_ct).loc[("testis", ["observed", "adj std residual", "fdr q-value"]),:]

# %% [markdown]
# # Larval Bulk

# %%
larval_ct = (
    pd.read_feather("../output/bulk2-rnaseq-wf/testis_ovary_counts.feather")
    .assign(flag_on=lambda x: x.Count >= 5)
    .groupby(["tissue", "chrom"])
    .flag_on.sum()
    .unstack()
    .reindex(index=["testis", "ovary"], columns="X,2L,2R,3L,3R,4,Y".split(","))
)  # type: pd.DataFrame

run_chisq(larval_ct).loc[("testis", ["observed", "adj std residual", "fdr q-value"]),:]

# %% [markdown]
# # Larval scRNA-Seq

# %%
sc_ct = (
    pd.read_feather("../output/seurat3-cluster-wf/aggegated_gene_counts_by_germ_soma.feather")
    .assign(cell_type=lambda x: x.cell_type.replace({"Germline": "Germline", "Cyst Lineage": "Somatic", "Other Somatic": "Somatic"}))
    .assign(flag_on=lambda x: x.Count >= 5)
    .groupby(["cell_type", "chrom"])
    .flag_on.sum()
    .unstack()
    .reindex(index=["Germline", "Somatic"], columns="X,2L,2R,3L,3R,4,Y".split(","))
) # type: pd.DataFrame

run_chisq(sc_ct).loc[("Germline", ["observed", "adj std residual", "fdr q-value"]),:]
