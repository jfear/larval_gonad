import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
colors_config = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        raw="../../output/seurat3-cluster-wf/raw_by_cluster.feather",
        zscores="../../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
    ),
    params=dict(
        color=colors_config["heatmap"],
        chrom_order=config["chrom_order"]
    )
)
