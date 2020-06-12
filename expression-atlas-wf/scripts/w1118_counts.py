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
        .assign(stage=lambda x: x.sample_ID.str.extract(r"(\w+)_\w+_r\d", expand=False))
        .assign(tissue=lambda x: x.sample_ID.str.extract(r"\w+_(\w+)_r\d", expand=False))
        .assign(rep=lambda x: x.sample_ID.str.extract(r"\w+_\w+_r(\d)", expand=False))
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
    cols = [
        "w1118_GO_m_r1",
        "w1118_GO_m_r2",
        "w1118_GO_m_r3",
        "w1118_GO_m_r4",
        "w1118_GO_f_r1",
        "w1118_GO_f_r2",
        "w1118_GO_f_r3",
        "w1118_GO_f_r4",
    ]

    gene_counts = (
        pd.read_feather(snakemake.input.counts)
        .set_index("YOgn")
        .reindex(columns=cols)
        .assign(FBgn=lambda x: x.index.map(lambda y: yogn2fbgn.get(y, np.nan)))
        .dropna()
        .set_index("FBgn", drop=True)
    )

    gene_counts.columns = [
        "adult_testis_r1",
        "adult_testis_r2",
        "adult_testis_r3",
        "adult_testis_r4",
        "adult_ovary_r1",
        "adult_ovary_r2",
        "adult_ovary_r3",
        "adult_ovary_r4",
    ]

    return gene_counts


if __name__ == "__main__":
    main()
