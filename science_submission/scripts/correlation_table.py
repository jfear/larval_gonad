import os
import pandas as pd

from larval_gonad.normalization import tpm


def main():
    df = pd.concat(
        [read_l3_sc(), read_l3_bulk(), read_adult_bulk()], sort=True, axis=1, join="inner"
    )  # type: pd.DataFrame

    df.corr(method="spearman").rename_axis("Spearman rho").to_csv(snakemake.output[0], sep="\t")


def read_l3_sc() -> pd.DataFrame:
    return pd.pivot(
        (
            pd.read_feather(snakemake.input.larval_scrnaseq)
            .assign(sample_ID=lambda x: x.cluster.str.cat(x.rep, sep="_"))
            .assign(sample_ID=lambda x: "l3_scRNAseq_" + x.sample_ID)
        ),
        index="FBgn",
        columns="sample_ID",
        values="TPM",
    )


def read_l3_bulk() -> pd.DataFrame:
    df = pd.read_csv(snakemake.input.larval_bulk, sep="\t", index_col=0).rename_axis("FBgn")
    cols = [f"l3_bulk_{x}" for x in df.columns]
    df.columns = cols
    return df


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
    )
    return tpm(dat, gene_lens).dropna()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                larval_scrnaseq="../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather",
                larval_bulk="../output/bulk2-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv",
                adult_bulk="../output/expression-atlas-wf/w1118_gene_counts.feather",
            ),
        )

    main()
