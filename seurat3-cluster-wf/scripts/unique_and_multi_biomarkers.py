import os

import matplotlib.pyplot as plt

from larval_gonad.config import read_config
from larval_gonad.plotting.biomarkers import plot_all_biomarkers


def main():
    lit_genes = read_config(snakemake.input.lit_genes)

    plt.style.use("science_base")
    fig = plt.figure(figsize=(6, 8))

    plot_all_biomarkers(
        snakemake.input.biomarkers, snakemake.input.zscores, lit_genes, fig=fig
    )

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                biomarkers="../output/seurat3-cluster-wf/combined_n3_biomarkers.feather",
                zscore="../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
                lit_genes="../config/literature_genes.yaml",
            ),
        )

    main()
