"""Plot heatmap of all genes.

Plots the tpm normalized zscores of all genes as a heatmap.
"""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns

BULK = snakemake.input.bulk
SCRNASEQ = snakemake.input.scrnaseq

OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass
# BULK = '../../output/bulk2-rnaseq-wf/rnaseq_aggregation/tpm_gene_level_counts.tsv'
# SCRNASEQ = '../../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather'

BULK_MAPPER = {
    "A1_OCP": "Ovary_rep1",
    "A2_OCP": "Ovary_rep2",
    "A3_OCP": "Ovary_rep3",
    "A4_OCP": "Ovary_rep4",
    "A5_TCP": "Testis_rep1",
    "A6_TCP": "Testis_rep2",
    "A7_TCP": "Testis_rep3",
    "A8_TCP": "Testis_rep4",
}


def main():
    df = get_data()

    plt.style.use("scripts/figure_styles.mplstyle")
    fig = plt.figure(figsize=(10, 10))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        df,
        xticklabels=True,
        yticklabels=True,
        vmin=0,
        vmax=1,
        square=True,
        rasterized=False,
        ax=ax,
        annot=True,
        annot_kws=dict(fontsize=5, fontweight="bold"),
        cbar_ax=cax,
        cbar_kws=dict(label="Spearman Correlation", orientation="horizontal", ticks=[0, 0.5, 1]),
    )

    # Clean up X axis
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_tick_params(pad=0, length=2)
    plt.setp(ax.get_xticklabels(), fontsize=7, rotation=90)

    # Clean up Y axis
    ax.set_ylabel("")
    ax.yaxis.set_ticks_position("right")
    ax.yaxis.set_tick_params(pad=0, length=2)
    plt.setp(ax.get_yticklabels(), fontsize=7, rotation=0, ha="left", va="center")

    # Clean up color bar
    plt.setp(cax.xaxis.label, fontsize=6)
    plt.setp(cax.get_xticklabels(), fontsize=5)
    cax.xaxis.set_tick_params(pad=0, length=2)

    fig.savefig(OUTPUT_FILE, bbox_inches="tight")


def get_data():
    bulk = (
        pd.read_csv(BULK, sep="\t")
        .query("Geneid.str.startswith('FBg')", engine="python")
        .set_index("Geneid")
        .rename_axis("FBgn")
        .rename(BULK_MAPPER, axis=1)
    )

    sc = pd.pivot_table(
        pd.read_feather(SCRNASEQ), index="FBgn", columns=["cluster", "rep"], values="TPM"
    )
    sc.columns = sc.columns.map("_".join)

    df = sc.join(bulk, how="inner")

    _corr = df.corr(method="spearman")

    # calculate linkages
    link = linkage(_corr.values, "average")
    tree = dendrogram(link, no_plot=True)
    leaves = tree["leaves"]

    _corr = _corr.iloc[leaves, leaves]
    for i, j in zip(*np.tril_indices(_corr.shape[0], k=-1)):
        _corr.iloc[i, j] = np.nan

    return _corr


if __name__ == "__main__":
    main()
