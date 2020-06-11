""" Data for Figure 2

prop_reads: is the number of reads (from expressed genes, > 5) mapping to
each chromosome arm scaled by the total number of reads from each sample
divided by the number of expressed genes (> 5) on that chromosome arm divided
by 1e3.

pct_expressed: is the number of expressed genes (> 5) divided by the number
of expressed genes on that chromosome times 100.

"""
import os
from typing import List
from collections import ChainMap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.normalization import tpm
from larval_gonad import plotting

NAMES = ["Average TPM", r"% Gene Expressed"]


def main():
    df = munge_prop_reads(NAMES[0])
    df.to_feather(snakemake.output.prop_reads)

    df = munge_pct_expressed(NAMES[1])
    df.to_feather(snakemake.output.pct_expressed)


def munge_prop_reads(name):
    counts = (
        pd.concat(
            [
                pd.read_feather(snakemake.input.adult_bulk),
                pd.read_feather(snakemake.input.L3_bulk),
                pd.read_feather(snakemake.input.L3_sc_agg),
                pd.read_feather(snakemake.input.L3_sc_by_class),
            ],
            sort=False,
        )
        .fillna("None")
        .assign(
            cell_type=lambda x: x.cell_type.map(
                {
                    "None": "None",
                    "Germline": "Germline",
                    "Other Somatic": "Somatic",
                    "Cyst Lineage": "Somatic",
                }
            )
        )
    )
    counts_tpm = add_tpm(counts)
    expressed = counts.loc[counts.Count > 0, "FBgn"].unique().tolist()

    fbgn2chrom = get_fbgn2chrom(expressed).value_counts().rename("num_genes")
    fbgn2chrom.index.name = "chrom"

    prop_counts = (
        counts_tpm.groupby(["chrom", "stage", "data_source", "cell_type", "rep", "tissue"])
        .TPM.sum()
        .div(fbgn2chrom, axis=0)
        .rename(name)
    )

    return prop_counts.to_frame().reset_index()


def add_tpm(counts):
    counts_matrix = pd.pivot_table(
        counts, index="FBgn", columns="sample_ID", values="Count", aggfunc="first"
    ).fillna(0)

    gene_lens = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
    )

    _tpm = (
        tpm(counts_matrix, gene_lens).dropna().reset_index().melt(id_vars="FBgn", value_name="TPM")
    )
    return counts.merge(_tpm, on=["FBgn", "sample_ID"])


def munge_pct_expressed(name):
    counts = (
        pd.concat(
            [
                pd.read_feather(snakemake.input.adult_bulk),
                pd.read_feather(snakemake.input.L3_bulk),
                pd.read_feather(snakemake.input.L3_sc_agg),
                pd.read_feather(snakemake.input.L3_sc_by_class),
            ],
            sort=False,
        )
        .fillna("None")
        .assign(
            cell_type=lambda x: x.cell_type.map(
                {
                    "None": "None",
                    "Germline": "Germline",
                    "Other Somatic": "Somatic",
                    "Cyst Lineage": "Somatic",
                }
            )
        )
        .groupby(["FBgn", "chrom", "stage", "data_source", "cell_type", "rep", "tissue"])
        .Count.sum()
        .reset_index()
        .assign(flag_on=lambda x: x.Count > 5)
    )

    expressed_FBgns = counts[counts.flag_on].FBgn.unique()

    fbgn2chrom = get_fbgn2chrom(expressed_FBgns).value_counts().rename("num_genes")
    fbgn2chrom.index.name = "chrom"

    pct_expressed = (
        counts.groupby(["chrom", "stage", "data_source", "cell_type", "rep", "tissue"])
        .flag_on.sum()
        .div(fbgn2chrom, axis=0)
        .mul(100)
        .rename(name)
    )

    return pct_expressed.to_frame().reset_index()


def get_fbgn2chrom(expressed):
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
        .reindex(expressed)
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        COLORS = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="paper_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                adult_bulk="../output/expression-atlas-wf/w1118_gene_counts.feather",
                L3_bulk="../output/bulk2-rnaseq-wf/testis_ovary_counts.feather",
                L3_sc_agg="../output/seurat3-cluster-wf/aggegated_gene_counts.feather",
                L3_sc_by_class="../output/seurat3-cluster-wf/aggegated_gene_counts_by_germ_soma.feather",
            ),
            params=dict(colors=COLORS),
        )

    main()
