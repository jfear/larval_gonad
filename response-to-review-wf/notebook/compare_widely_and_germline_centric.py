# %%
import joblib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display, Markdown

# %% [markdown]
# ## Overall Counts
# %%
def get_genes_on(cluster: str):
    cell_ids = (
        pd.read_feather("../../output/seurat3-cluster-wf/combined_n3_clusters.feather")
        .query("cluster == @cluster")
        .cell_id.to_list()
    )
    genes_on = (
        pd.read_feather(
            "../../output/cellselection-wf/raw.feather", columns=["FBgn"] + cell_ids
        )
        .set_index("FBgn")
        .pipe(lambda df: df > 0)
        .mean(axis=1)
        .pipe(lambda sr: sr[sr > 0.1])
        .index.to_list()
    )
    return set(genes_on)


# %%
annot = pd.read_feather("../../references/gene_annotation_dmel_r6-26.feather")

expressed = set(joblib.load("../../output/cellselection-wf/expressed_genes.pkl"))
gonia_on = get_genes_on("G")
spermatocyte_on = get_genes_on("EPS") | get_genes_on("MPS") | get_genes_on("LPS")

widely_expressed = set(
    joblib.load("../../output/cellselection-wf/commonly_expressed_genes.pkl")
)
germline_centric = set(
    joblib.load("../../output/response-to-review-wf/germline_centric_subset.pkl")
)
not_in_germline = widely_expressed - germline_centric
not_in_germline_genes = annot.query("FBgn in @not_in_germline").gene_symbol.to_list()

# %%
print(
    f"""
Number Expressed:\t\t{len(expressed):,}
Number Gonia Expressed:\t\t{len(gonia_on):,}
Number Spermatocyte Expressed:\t\t{len(spermatocyte_on):,}
Number Widely Expressed:\t{len(widely_expressed):,}
Number Germline Centric:\t{len(germline_centric):,}
Number Genes In Both:\t\t{len(widely_expressed.intersection(germline_centric)):,}
Genes in Widely not in Centric: {", ".join(not_in_germline_genes)}
"""
)

# %% [markdown]
# ## Proportion Genes On X and 4

# %%
display(Markdown("### All Genes"))
(annot.groupby("FB_chrom").size() / annot.shape[0])[
    ["X", "2L", "2R", "3L", "3R", "4", "Y"]
]

# %%
display(Markdown("### Expressed Genes"))
(
    annot.query("FBgn in @expressed").groupby("FB_chrom").size()
    / annot.query("FBgn in @expressed").shape[0]
)[["X", "2L", "2R", "3L", "3R", "4", "Y"]]

# %%
display(Markdown("### Gonia Expressed Genes"))
(
    annot.query("FBgn in @gonia_on").groupby("FB_chrom").size()
    / annot.query("FBgn in @gonia_on").shape[0]
)[["X", "2L", "2R", "3L", "3R", "4", "Y"]]


# %%
display(Markdown("### Spermatocyte Expressed Genes"))
(
    annot.query("FBgn in @spermatocyte_on").groupby("FB_chrom").size()
    / annot.query("FBgn in @spermatocyte_on").shape[0]
)[["X", "2L", "2R", "3L", "3R", "4", "Y"]]
# %%
display(Markdown("### Widely Expressed Genes"))
(
    annot.query("FBgn in @widely_expressed").groupby("FB_chrom").size()
    / annot.query("FBgn in @widely_expressed").shape[0]
)[["X", "2L", "2R", "3L", "3R", "4", "Y"]]

# %%
display(Markdown("### Germline Centric Genes"))
(
    annot.query("FBgn in @germline_centric").groupby("FB_chrom").size()
    / annot.query("FBgn in @germline_centric").shape[0]
)[["X", "2L", "2R", "3L", "3R", "4", "Y"]]
