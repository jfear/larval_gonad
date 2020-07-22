import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
colors = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        tpm="../../output/seurat3-cluster-wf/tpm_by_cluster.feather",
    ),
    params=dict(
        cluster_order=config["cluster_order"],
        cluster_names=config["cluster_names"],
        chrom_order=config["chrom_order"],
        cluster_colors=colors["clusters"],
    ),
)
