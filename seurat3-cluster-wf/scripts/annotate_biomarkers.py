"""Convert cluster to ordered cluster annotations."""
import os
import pandas as pd


def main():
    (
        pd.read_feather(snakemake.input[0])
        .assign(
            cluster=lambda df: pd.Categorical(
                df.cluster.map(dict(snakemake.params.cluster_annot)),
                categories=snakemake.params.cluster_order,
                ordered=True,
            )
        )
        .to_feather(snakemake.output[0])
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input="../output/seurat3-cluster-wf/combined_n3_biomarkers.feather",
            params=dict(
                color="viridis",
                cluster_annot=config["cluster_annot"],
                cluster_order=config["cluster_order"],
            ),
        )

    main()
