import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
colors = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        clusters="../../output/seurat3-cluster-wf/combined_n3_clusters.feather",
        all_genes="../../output/chromosome-expression-wf/cell_level_chromosome_proportions_all_expressed_genes.feather",
        widely_expressed_genes="../../output/chromosome-expression-wf/cell_level_chromosome_proportions_widely_expressed_genes.feather",
        tau_genes="../../output/chromosome-expression-wf/cell_level_chromosome_proportions_tau.feather",
        tsps_genes="../../output/chromosome-expression-wf/cell_level_chromosome_proportions_tsps.feather",
    ),
    params=dict(
        chrom_order=config["chrom_order"],
        cluster_order=config["cluster_order"],
        cluster_names=config["cluster_names"],
        cluster_colors=colors["clusters"],
    ),
)
