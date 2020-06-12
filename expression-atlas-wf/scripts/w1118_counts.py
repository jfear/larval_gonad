import os
import numpy as np
import pandas as pd

from larval_gonad.io import pickle_load


def main():
    fbgn2chrom = get_fbgn2chrom()
    yogn2fbgn = pickle_load(snakemake.input.orthologs)
    adult_counts = get_counts(yogn2fbgn)

    df = (
        adult_counts.reset_index()
        .melt(id_vars="FBgn", var_name="sample_ID", value_name="Count")
        .assign(strain="w1118")
        .assign(
            stage=lambda x: x.sample_ID.str.extract(
                r"w1118_(\w+)_\w+_r\d", expand=False
            )
        )
        .assign(
            tissue=lambda x: x.sample_ID.str.extract(
                r"w1118_\w+_(\w+)_r\d", expand=False
            )
        )
        .assign(
            rep=lambda x: x.sample_ID.str.extract(r"w1118_\w+_\w+_r(\d)", expand=False)
        )
        .assign(data_source="RNA-Seq")
        .set_index("FBgn")
        .join(fbgn2chrom, how="inner")
        .reset_index()
    )

    df.to_feather(snakemake.output[0])


def get_fbgn2chrom():
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
    )


def get_counts(yogn2fbgn):
    gene_counts = (
        pd.read_feather(snakemake.input.counts)
        .assign(FBgn=lambda x: x.YOgn.map(lambda y: yogn2fbgn.get(y, np.nan)))
        .dropna()
        .set_index("FBgn", drop=True)
        .drop("YOgn", axis=1)
    )

    return gene_counts


if __name__ == "__main__":
    main()
