import os
import re

import pandas as pd


def main():
    df = pd.concat(
        [cluster_names(), cluster_number(), read_cluster(), pct_cells(), cluster_corr()],
        axis=1,
        sort=False,
    )  # type: pd.DataFrame

    df.reset_index().to_csv(snakemake.output[0], sep="\t", index=False)


def cluster_names():
    names = [re.sub(r"\s\(.*\)", "", x) for x in snakemake.params.legend_names]
    sr = pd.Series(dict(zip(snakemake.params.cluster_order, names))).rename("Description")
    sr.index.name = "cluster"
    return sr


def cluster_number():
    sr = pd.Series({v: k for k, v in snakemake.params.cluster_annot.items()}).rename(
        "Seurat_cluster_number"
    )
    sr.index.name = "cluster"
    return sr.reindex(snakemake.params.cluster_order)


def read_cluster():
    df = (
        pd.read_feather(snakemake.input.clusters)
        .set_index("cell_id")
        .groupby(["cluster", "rep"])
        .size()
        .unstack()
    )

    cols = [f"nCells_{x}" for x in df.columns]
    df.columns = cols

    return df


def pct_cells():
    df = read_cluster()
    pct = df.div(df.sum(axis=1), axis=0) * 100
    cols = [x.replace("nCells", "pctCells") for x in pct.columns]
    pct.columns = cols
    return pct


def cluster_corr():
    tpm = pd.read_feather(snakemake.input.tpm).set_index(["FBgn", "cluster"]).unstack()
    corr = tpm.corr(method="spearman")
    cols = [f"Spearman_rho_{x}" for x in corr.columns.droplevel().astype(str)]
    corr.columns = cols
    corr.index = corr.index.droplevel()
    return corr


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("../../config/common.yaml")

        snakemake = snakemake_debug(
            workdir="paper_submission",
            input=dict(
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
                tpm="../output/seurat3-cluster-wf/tpm_by_cluster.feather",
            ),
            params=dict(
                cluster_annot=config["cluster_annot"],
                cluster_order=config["cluster_order"],
                legend_names=config["legend_names"],
            ),
        )

    main()
