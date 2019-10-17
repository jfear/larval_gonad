import os

import matplotlib.pyplot as plt

from larval_gonad.config import read_config
from larval_gonad.plotting.literature import plot_lit_evidence_zscore_soma_profile

def main():
    plt.style.use("science_base")

    _, axes = plt.subplots(1, 2, figsize=(1.6, 3.4))

    plot_lit_evidence_zscore_soma_profile(
        gene_metadata=snakemake.input.gene_metadata,
        lit_evidence=snakemake.input.lit_evidence,
        zscore_by_cluster_rep=snakemake.input.zscore_by_cluster_rep,
        soma_clusters=snakemake.config['soma'],
        axes=axes
    )

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        config = read_config("config/common.yaml")
        TAG = config["tag"]

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_metadata=f"../references/gene_annotation_dmel_{TAG}.feather",
                lit_evidence="../data/external/miriam/lit_gene_dummy_vars.tsv",
                tpm_by_cluster="../output/seurat3-cluster-wf/tpm_by_cluster.feather"
            ),
            config=config
        )

    main()
