import os
import pandas as pd

from larval_gonad.normalization import tpm


def main():
    df = pd.concat(
        [read_l3_sc(), read_l3_bulk()], sort=True, axis=1, join="inner"
    )  # type: pd.DataFrame

    df.corr(method="spearman").rename_axis("Spearman rho").to_csv(snakemake.output[0], sep="\t")


def read_l3_sc() -> pd.DataFrame:
    gene_lengths = pd.read_feather(snakemake.input.gene_annot).set_index("FBgn")["length"].squeeze()

    raw = (
        pd.read_feather(snakemake.input.larval_scrnaseq)
        .groupby(["FBgn", "rep"])
        .Count.sum()
        .unstack()
    )
    norm = tpm(raw, gene_lengths).dropna()
    norm.columns = [f"l3_scRNAseq_{x}" for x in norm.columns]

    return norm


def read_l3_bulk() -> pd.DataFrame:
    df = pd.read_csv(snakemake.input.larval_bulk, sep="\t", index_col=0).rename_axis("FBgn")
    cols = [f"l3_bulk_{x}" for x in df.columns]
    df.columns = cols
    return df


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
                larval_scrnaseq="../../output/seurat3-cluster-wf/aggegated_gene_counts.feather",
                larval_bulk="../../output/bulk2-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv",
            )
        )

    main()
