import os
import pandas as pd
from scipy.stats import spearmanr

from larval_gonad.normalization import tpm
from larval_gonad.constants import L3_BULK, L3_SC


def main():
    df = pd.concat(
        [read_l3_sc(), read_l3_bulk()], sort=True, axis=1, join="inner"
    )  # type: pd.DataFrame

    _, pval = spearmanr(df.values, axis=0)
    pd.DataFrame(pval, index=df.columns, columns=df.columns).rename_axis("Spearman pval").to_csv(
        snakemake.output[0], sep="\t"
    )


def read_l3_sc() -> pd.DataFrame:
    gene_lengths = pd.read_feather(snakemake.input.gene_annot).set_index("FBgn")["length"].squeeze()
    raw = (
        pd.read_feather(snakemake.input.larval_scrnaseq)
        .groupby(["FBgn", "rep"])
        .Count.sum()
        .unstack()
        .rename(columns=L3_SC)
    )
    norm = tpm(raw, gene_lengths).dropna()
    return norm


def read_l3_bulk() -> pd.DataFrame:
    return (
        pd.read_csv(snakemake.input.larval_bulk, sep="\t", index_col=0)
        .rename_axis("FBgn")
        .rename(columns=L3_BULK)
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
                larval_scrnaseq="../../output/seurat3-cluster-wf/aggegated_gene_counts.feather",
                larval_bulk="../../output/bulk-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv",
            )
        )

    main()
