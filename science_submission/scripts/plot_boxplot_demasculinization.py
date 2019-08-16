"""Demasculinization of the X by cell type.

Look at how and where testis biased expressed genes are expressed in our
different cell-type groups.

"""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns


BULK_DEG = snakemake.input["bulk"]
TPM = snakemake.input["tpm"]
GENE_ANNOTATION = snakemake.input["gene_annotation"]

CLUSTER_COLORS = snakemake.params["colors"]

BULK_DEMASCULINIZATION_FILE = snakemake.output["demas"]
TESTIS_BIASED_FILE = snakemake.output["testis"]
OVARY_BIASED_FILE = snakemake.output["ovary"]

# Debug settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass
# BULK_DEG = '../../output/bulk2-rnaseq-wf/deg/bulk_testis_vs_ovary.tsv'
# TPM = '../../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather'
# GENE_ANNOTATION = '../../references/gene_annotation_dmel_r6-24.feather'
# import yaml
# CLUSTER_COLORS = yaml.full_load(open("../../config/colors.yaml"))['clusters']

BARPLOT_DEFAULTS = dict(order=["X", "2L", "2R", "3L", "3R", "4"], color="C0")

BOXPLOT_DEFAULTS = dict(
    hue="cluster",
    palette=CLUSTER_COLORS,
    showfliers=False,
    notch=True,
    order=["X", "2L", "2R", "3L", "3R", "4"],
    linewidth=0.5,
)


def main():
    plt.style.use("scripts/figure_styles.mplstyle")
    fbgn2chrom = get_fbgn2chrom(GENE_ANNOTATION)

    bulk = pd.read_csv(BULK_DEG, sep="\t").set_index("FBgn").join(fbgn2chrom)
    bulk["testis_biased"] = (bulk.padj < 0.01) & (bulk.log2FoldChange > 0)
    bulk["ovary_biased"] = (bulk.padj < 0.01) & (bulk.log2FoldChange < 0)

    aggregated = bulk.groupby("chrom")[["testis_biased", "ovary_biased"]].mean().reset_index()

    tpm_testis = get_tpm_subset(bulk[bulk.testis_biased].index.tolist(), fbgn2chrom)
    tpm_ovary = get_tpm_subset(bulk[bulk.ovary_biased].index.tolist(), fbgn2chrom)

    fig = plt.figure(figsize=(8, 12))
    gs = GridSpec(3, 4, width_ratios=[1, .2, 1, .2])
    ax_testis_biased = fig.add_subplot(gs[0, 0:2])
    ax_ovary_biased = fig.add_subplot(gs[0, 2:])
    ax_testis_cluster = fig.add_subplot(gs[1, :-1])
    ax_ovary_cluster = fig.add_subplot(gs[2, :-1])

    sns.barplot("chrom", "testis_biased", data=aggregated, ax=ax_testis_biased, **BARPLOT_DEFAULTS)
    ax_testis_biased.set(ylim=(0, 1), ylabel="Percent Genes", xlabel="", title="Bulk\nTestis Biased")

    sns.barplot("chrom", "ovary_biased", data=aggregated, ax=ax_ovary_biased, **BARPLOT_DEFAULTS)
    ax_ovary_biased.set(ylim=(0, 1), ylabel="", xlabel="", yticklabels=[], title="Bulk\nOvary Biased")

    sns.boxplot("chrom", "TPM", data=tpm_testis, ax=ax_testis_cluster, **BOXPLOT_DEFAULTS)
    ax_testis_cluster.set(ylabel="log10(TPM + 1)", xlabel="", title="Testis Biased Genes")
    ax_testis_cluster.legend(loc="upper left", bbox_to_anchor=[1, 1])

    sns.boxplot("chrom", "TPM", data=tpm_ovary, ax=ax_ovary_cluster, **BOXPLOT_DEFAULTS)
    ax_ovary_cluster.set(ylabel="log10(TPM + 1)", xlabel="", title="Ovary Biased Genes")
    ax_ovary_cluster.legend(loc="upper left", bbox_to_anchor=[1, 1])


def get_fbgn2chrom(fname):
    return pd.read_feather(fname).set_index("FBgn").FB_chrom.rename("chrom")


def get_tpm_subset(fbgns, fbgn2chrom):
    return (
        pd.read_feather(TPM)
        .groupby(["FBgn", "cluster"])
        .TPM.median()
        .pipe(lambda x: np.log10(x + 1))
        .to_frame()
        .reset_index()
        .query(f"FBgn == {fbgns}")
        .merge(fbgn2chrom, left_on="FBgn", right_index=True)
    )


if __name__ == "__main__":
    main()
