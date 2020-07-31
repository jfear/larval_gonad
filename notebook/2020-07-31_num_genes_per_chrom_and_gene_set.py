"""Create table with the number of genes per chromosome for different gene subsets.

rows: major chromosomes
cols: gene sets (FlyBase r6.26, Expressed, Widely Expressed, Tau, TSPS)
values: number of genes
"""
# %%
import pandas as pd
import joblib

from larval_gonad.config import read_config

# %%
# File Names
annot = "../references/gene_annotation_dmel_r6-26.feather"
expressed = "../output/cellselection-wf/expressed_genes.pkl"
widely_expressed = "../output/cellselection-wf/commonly_expressed_genes.pkl"
tau = "../output/expression-atlas-wf/dmel_male_tau_fbgns.pkl"
tsps = "../output/expression-atlas-wf/dmel_male_tsps_fbgns.pkl"

# %%
# Load gene annotation
chrom_order = read_config("../config/common.yaml")["chrom_order"]
genes = (
    pd.read_feather(annot, columns=["FBgn", "FB_chrom"])
    .rename(columns={"FB_chrom": "chrom"})
    .query("chrom == @chrom_order")
    .set_index("FBgn")
)

# %%
def chrom_count(file_name, title, df: pd.DataFrame):
    return (
        df.reindex(joblib.load(file_name))
        .dropna()
        .groupby("chrom")
        .size()
        .rename(title)
    )


# %%
# Create gene count table [chromosome x gene set]
counts_table = (
    pd.concat(
        [
            genes.groupby("chrom").size().rename("FlyBase r6.26"),
            chrom_count(expressed, "Expressed", genes),
            chrom_count(widely_expressed, "Widely Expressed", genes),
            chrom_count(tau, "Tau", genes),
            chrom_count(tsps, "TSPS", genes),
        ],
        axis=1,
    )
    .loc[chrom_order, :]  # Re-order rows to be consistent
    # NaNs from Y chrom turned count to a float. Change back to ints.
    .fillna(0)
    .astype(int)
)

counts_table.index.name = "Chromosome Element"


# %%
counts_table.to_csv(
    "../output/notebook/2020-07-31_num_genes_per_chrom_and_gene_set.tsv", sep="\t"
)
