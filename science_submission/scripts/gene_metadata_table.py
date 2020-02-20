import os
import pickle

import pandas as pd

from larval_gonad.io import feather_to_cluster_rep_matrix


def main():
    df = pd.concat(
        [read_gene_annot(), read_biomarkers(), read_zscore(), read_raw(), read_tpm()],
        axis=1,
        sort=False,
    )
    df.rename_axis("FBgn", inplace=True)
    df.to_csv(snakemake.output[0], sep="\t")


def read_gene_annot():
    df = (
        pd.read_feather(snakemake.input.gene_annot)
        .assign(chrom=lambda x: x.FB_chrom)
        .drop(["FB_chrom", "UCSC_chrom"], axis=1)
        .set_index("FBgn")
        .loc[:, ["gene_symbol", "chrom", "start", "end", "length", "strand"]]
    )

    df["tau"] = False
    df.loc[df.index.isin(pickle.load(open(snakemake.input.tau, "rb"))), "tau"] = True

    df["tsps"] = False
    df.loc[df.index.isin(pickle.load(open(snakemake.input.tsps, "rb"))), "tsps"] = True

    return df


def read_biomarkers():
    return (
        pd.read_feather(snakemake.input.biomarkers)
        .assign(
            cluster=lambda x: pd.Categorical(
                x.cluster.map(lambda y: snakemake.params["cluster_annot"][y]),
                ordered=True,
                categories=snakemake.params["cluster_order"],
            )
        )
        .query("p_val_adj <= 0.01")
        .groupby("FBgn")
        .apply(lambda x: "|".join(x.cluster.sort_values().values))
        .rename("marker_for_cell_type")
    )


def read_raw():
    df = feather_to_cluster_rep_matrix(snakemake.input.raw)
    df.columns = [f"raw_{x}_{y}" for x, y in df.columns.to_flat_index()]
    return df


def read_tpm():
    df = feather_to_cluster_rep_matrix(snakemake.input.tpm)
    df.columns = [f"tpm_{x}_{y}" for x, y in df.columns.to_flat_index()]
    return df


def read_zscore():
    df = feather_to_cluster_rep_matrix(snakemake.input.zscore)
    df.columns = [f"zscore_{x}_{y}" for x, y in df.columns.to_flat_index()]
    return df


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("../../config/common.yaml")

        snakemake = snakemake_debug(
            input=dict(
                gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
                biomarkers="../../output/seurat3-cluster-wf/combined_n3_biomarkers.feather",
                raw="../../output/seurat3-cluster-wf/raw_by_cluster_rep.feather",
                tpm="../../output/seurat3-cluster-wf/tpm_by_cluster_rep.feather",
                zscore="../../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
                tau="../../output/expression-atlas-wf/tau_housekeeping/male_fbgns.pkl",
                tsps="../../output/expression-atlas-wf/tsps_housekeeping/male_fbgns.pkl",
            ),
            params=dict(
                cluster_annot=config["cluster_annot"], cluster_order=config["cluster_order"]
            ),
        )

    main()
