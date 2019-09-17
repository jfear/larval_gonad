"""Create UMAP panel of individual reps and combined"""
import os

import matplotlib.pyplot as plt

from larval_gonad.config import read_config
from larval_gonad.plotting.umap import plot_umap_panel


def main():

    plt.style.use("science_base")
    _, axes = plt.subplots(1, 4, figsize=(3.5, 1))

    plot_umap_panel(
        umap_data=snakemake.input.umap,
        cluster_data=snakemake.input.clusters,
        colors=snakemake.params.colors,
        axes=axes,
    )

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        color_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
            ),
            params=dict(colors=color_config["clusters"]),
        )

    main()
