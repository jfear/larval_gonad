import os
import pandas as pd

from larval_gonad.normalization import tpm
from larval_gonad.constants import ADULT_BULK, L3_BULK, L3_SC


def main():
    df = pd.concat(
        [read_l3_sc(), read_l3_bulk()], sort=True, axis=1, join="inner"
    )  # type: pd.DataFrame

    df.corr(method="spearman").rename_axis("Spearman rho").to_csv(snakemake.output[0], sep="\t")


def read_l3_sc() -> pd.DataFrame:
    return (
        pd.pivot_table(
            (
                pd.read_feather(snakemake.input.larval_scrnaseq)
                .sort_values(["cluster", "FBgn"])
                .assign(sample_ID=lambda x: x.rep.map(L3_SC))
                .assign(sample_ID=lambda x: x.sample_ID.str.cat(x.cluster, sep="_"))
            ),
            index="FBgn",
            columns=["cluster", "rep", "sample_ID"],
            values="TPM",
            aggfunc="first",
        )
        .sort_values(["cluster", "rep"], axis=1)
        .droplevel([0, 1], axis=1)
    )


def read_l3_bulk() -> pd.DataFrame:
    return (
        pd.read_csv(snakemake.input.larval_bulk, sep="\t", index_col=0)
        .rename_axis("FBgn")
        .rename(columns=L3_BULK)
        .sort_index(axis=1)
    )


def read_adult_bulk() -> pd.DataFrame:
    gene_lens = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    dat = pd.pivot_table(
        pd.read_feather(snakemake.input.adult_bulk),
        index="FBgn",
        columns="sample_ID",
        values="Count",
    ).sort_index(axis=1)

    return tpm(dat, gene_lens).dropna()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from larval_gonad.config import read_config
        from larval_gonad.debug import snakemake_debug

        config = read_config("../../config/common.yaml")

        snakemake = snakemake_debug(
            input=dict(
                gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
                larval_scrnaseq="../../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather",
                larval_bulk="../../output/bulk-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv",
                adult_bulk="../../output/expression-atlas-wf/w1118_gene_counts.feather",
            ),
            params=dict(cluster_order=config["cluster_order"]),
        )

    main()
