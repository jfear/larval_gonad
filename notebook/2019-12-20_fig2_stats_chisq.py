# %%
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display

from larval_gonad.stats import run_chisq

pd.options.display.max_rows = 200

# %%
try:
    os.chdir(os.path.join(os.getcwd(), "notebook"))
    print(os.getcwd())
except:
    pass

# %% [markdown]
# # Adult Bulk

# %%
adult_bulk = pd.read_feather("../output/expression-atlas-wf/w1118_gene_counts.feather")
adult_ct = (
    adult_bulk
    .assign(flag_on=lambda x: x.Count >= 5)
    .groupby(["tissue", "chrom"])
    .flag_on.sum()
    .unstack()
    .reindex(index=["testis", "ovary"], columns="X,2L,2R,3L,3R,4,Y".split(","))
)  # type: pd.DataFrame

# run_chisq(adult_ct).loc[("testis", ["observed", "adj std residual", "fdr q-value"]),:]
run_chisq(adult_ct).loc[(slice(None), ["observed", "adj std residual", "fdr q-value"]),:]

# %%
fbgn2symbol = pd.read_feather("../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()
ovary_y_fbgns = adult_bulk.query("tissue == 'ovary' & chrom == 'Y' & Count > 0")
display(ovary_y_fbgns.merge(fbgn2symbol, on="FBgn").set_index(["FBgn", "gene_symbol"]).sort_values("Count", ascending=False))
fbgn2symbol.reindex(ovary_y_fbgns.FBgn.unique()).to_frame()

# %%

# %% [markdown]
# # Larval Bulk

# %%
larval_bulk = pd.read_feather("../output/bulk2-rnaseq-wf/testis_ovary_counts.feather")
larval_ct = (
    larval_bulk
    .assign(flag_on=lambda x: x.Count >= 5)
    .groupby(["tissue", "chrom"])
    .flag_on.sum()
    .unstack()
    .reindex(index=["testis", "ovary"], columns="X,2L,2R,3L,3R,4,Y".split(","))
)  # type: pd.DataFrame

run_chisq(larval_ct).loc[("testis", ["observed", "adj std residual", "fdr q-value"]),:]

# %%
fbgn2symbol = pd.read_feather("../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()
ovary_y_fbgns = larval_bulk.query("tissue == 'ovary' & chrom == 'Y' & Count > 0")
display(ovary_y_fbgns.merge(fbgn2symbol, on="FBgn").set_index(["FBgn", "gene_symbol"]).sort_values("Count", ascending=False))
print(ovary_y_fbgns.Count.sum())
fbgn2symbol.reindex(ovary_y_fbgns.FBgn.unique()).to_frame()

# %% [markdown]
# # Larval scRNA-Seq

# %%
sc = (
    pd.read_feather("../output/seurat3-cluster-wf/aggegated_gene_counts_by_germ_soma.feather")
    .assign(cell_type=lambda x: x.cell_type.replace({"Germline": "Germline", "Cyst Lineage": "Somatic", "Other Somatic": "Somatic"}))
)

sc_ct = (
    sc
    .assign(flag_on=lambda x: x.Count >= 5)
    .groupby(["cell_type", "chrom"])
    .flag_on.sum()
    .unstack()
    .reindex(index=["Germline", "Somatic"], columns="X,2L,2R,3L,3R,4,Y".split(","))
) # type: pd.DataFrame

run_chisq(sc_ct).loc[("Germline", ["observed", "adj std residual", "fdr q-value"]),:]

# %%
fbgn2symbol = pd.read_feather("../references/gene_annotation_dmel_r6-26.feather", columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()
somatic_y_fbgns = sc.query("cell_type == 'Somatic' & chrom == 'Y' & Count > 0")
display(somatic_y_fbgns.merge(fbgn2symbol, on="FBgn").set_index(["FBgn", "gene_symbol"]).sort_values("Count", ascending=False))
print(somatic_y_fbgns.Count.sum())
fbgn2symbol.reindex(somatic_y_fbgns.FBgn.unique()).to_frame()


# %%
