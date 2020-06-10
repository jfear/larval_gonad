import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import chain

import larval_gonad.plotting


def _debug():
    from larval_gonad.mock import MockSnake
    from larval_gonad.config import read_config

    snakemake = MockSnake(
        input=dict(
            annotation="../../references/gene_annotation_dmel_r6-26.feather",
            zscore="../../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
        ),
        params=dict(
            msl_genes=read_config("../../config/other_lists.yaml")["msl_genes"],
            cluster_order=read_config("../../config/common.yaml")["cluster_order"],
        ),
    )


def main():
    plt.style.use("science_base")
    zscores = load_zscores()
    ax = sns.heatmap(
        zscores,
        vmin=-3,
        vmax=3,
        cmap="viridis",
        xticklabels=True,
        yticklabels=True,
        rasterized=True,
        cbar_kws=dict(label="Z-Score (TPM)", ticks=[-3, 3], orientation="horizontal"),
    )
    _cleanup_xaxis(ax, snakemake.params.cluster_order)
    ax.set_ylabel("")

    for out_format in snakemake.output:
        plt.savefig(out_format)


def _cleanup_xaxis(ax, cluster_order):
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        list(chain.from_iterable([("", x, "") for x in cluster_order])),
        ha="center",
        va="bottom",
    )

    # Add lines separating cell types
    for i in range(1, len(cluster_order)):
        ax.axvline(i * 3, color="w", ls="--", lw=0.5)

    return ax


def load_zscores():
    fbgn2symbol = (
        pd.read_feather(snakemake.input.annotation).set_index("FBgn").gene_symbol
    )
    zscores = (
        pd.read_feather(snakemake.input.zscore)
        .set_index(["FBgn", "cluster", "rep"])
        .unstack(level=[1, 2])
        .query("FBgn in @snakemake.params.msl_genes")
    )
    zscores.columns = zscores.columns.droplevel(0)
    zscores.index = zscores.index.map(fbgn2symbol)
    return zscores


if __name__ == "__main__":
    main()
