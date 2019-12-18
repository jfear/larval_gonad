"""Annotate Adult Testis Clusters

I want to use our cluster annotation from the larval testis to inform the
adult annotations. I am using a simple correlation metric.

"""
import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.normalization import tpm


def main():
    adult = munge_data(snakemake.input.adult_metadata, snakemake.input.adult)
    larval = munge_data(snakemake.input.larval_metadata, snakemake.input.larval)
    corr = pd.merge(normalize(adult), normalize(larval), on="FBgn").corr()
    sns.clustermap(corr)


def munge_data(metadata, data):
    clusters = pd.read_feather(metadata, columns=["cell_id", "cluster"]).assign(
        cluster=lambda x: x.cluster.astype(str)
    )

    norm = (
        pd.read_feather(data)
        .melt(id_vars="FBgn", var_name="cell_id", value_name="norm")
        .merge(clusters, on="cell_id")
        .groupby(["cluster", "FBgn"])
        .norm.sum()
    )

    return norm.unstack(level=0)


def normalize(df):
    gene_len = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    return tpm(df, gene_len).dropna()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", True):

        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="adult-testis-scRNAseq-wf",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                adult="../output/adult-testis-scRNAseq-wf/raw_Ral517.feather",
                adult_metadata="../output/adult-testis-scRNAseq-wf/Ral517_metadata.feather",
                larval="../output/cellselection-wf/raw.feather",
                larval_metadata="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
            ),
            wildcards=dict(sample="Ral517"),
        )

    main()
