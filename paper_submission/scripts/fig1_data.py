""" Data for Figure 2

Calculate the Average TPM per chromosome (expressed genes: TPM > 0)

"""
import os
from typing import List
from collections import ChainMap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.normalization import tpm
from larval_gonad import plotting


def main():
    avg_tpm_per_chrom = (
        pd.concat(
            [
                read_adult_counts(snakemake.input.w1118),
                read_adult_counts(snakemake.input.OreR),
                read_larval_counts(snakemake.input.L3_bulk),
            ],
            sort=False,
        )
        .pipe(add_tpm)
        .groupby(["stage", "strain", "tissue"])
        .apply(calculate_avg_tpm_per_chrom)
        .reset_index()
    )

    avg_tpm_per_chrom.to_feather(snakemake.output[0])


def read_adult_counts(file_name) -> pd.DataFrame:
    return pd.read_feather(file_name)


def read_larval_counts(file_name) -> pd.DataFrame:
    return (
        pd.read_feather(file_name)
        .assign(sex=lambda df: np.where(df.tissue == "ovary", "f", "m"))
        .assign(strain="w1118")
        .replace({"ovary": "L3_GO", "testis": "L3_GO"})
    )


def add_tpm(df: pd.DataFrame) -> pd.DataFrame:
    matrix_ = pd.pivot(df, index="FBgn", columns="sample_ID", values="Count").fillna(0)

    gene_lens = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    tpm_ = (
        tpm(matrix_, gene_lens)
        .dropna()
        .reset_index()
        .melt(id_vars="FBgn", value_name="TPM")
    )
    return df.merge(tpm_, on=["FBgn", "sample_ID"])


def map_chrom(expressed: list) -> pd.Series:
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
        .reindex(expressed)
    )


def calculate_avg_tpm_per_chrom(df: pd.DataFrame) -> pd.DataFrame:
    expressed_fbgns = df.loc[df.Count > 0, "FBgn"].unique().tolist()

    num_expressed_genes_per_chrom = (
        map_chrom(expressed_fbgns)
        .value_counts()
        .rename("num_genes")
        .rename_axis("chrom")
    )

    avg_tpm = (
        df.groupby(["chrom", "rep", "sex"])
        .TPM.sum()
        .div(num_expressed_genes_per_chrom, axis=0)
        .rename("Avg TPM Per Chromosome")
        .to_frame()
    )

    return avg_tpm


if __name__ == "__main__":
    main()
