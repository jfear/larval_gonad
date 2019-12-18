"""Annotate Adult Testis Clusters

I want to use our cluster annotation from the larval testis to inform the
adult annotations. I am using a simple correlation metric.

"""
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from more_itertools import flatten

from larval_gonad.config import read_config
from larval_gonad.normalization import tpm


def main():
    corr = plot_corr()

    lit_genes = list(flatten(read_config(snakemake.input.lit_genes).values()))
    biomarkers = pd.read_feather(snakemake.input.adult_biomarkers).query(
        f"FBgn == {lit_genes} & p_val_adj <= 0.01"
    )

    biomarkers.groupby(["FBgn", "gene_symbol"]).apply(
        lambda x: "|".join(x.cluster.sort_values())
    ).rename("clusters").reset_index().set_index("FBgn").reindex(lit_genes).dropna()


def plot_corr():
    adult = munge_data(snakemake.input.adult_metadata, snakemake.input.adult)
    larval = munge_data(snakemake.input.larval_metadata, snakemake.input.larval)
    corr = pd.merge(normalize(adult), normalize(larval), on="FBgn").corr()
    return sns.clustermap(corr)


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
                lit_genes="../config/literature_genes.yaml",
                adult="../output/adult-testis-scRNAseq-wf/raw_Ral517.feather",
                adult_metadata="../output/adult-testis-scRNAseq-wf/Ral517_metadata.feather",
                adult_biomarkers="../output/adult-testis-scRNAseq-wf/Ral517_biomarkers.feather",
                larval="../output/cellselection-wf/raw.feather",
                larval_metadata="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
            ),
            wildcards=dict(sample="Ral517"),
        )

    main()
