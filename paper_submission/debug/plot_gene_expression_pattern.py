import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

builtins.snakemake = MockSnake(
    input=dict(
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        clusters="../../output/seurat3-cluster-wf/combined_n3_clusters.feather",
        raw="../../output/cellselection-wf/raw.feather",
    ),
    params=dict(colors=read_config("../../config/colors.yaml")["clusters"]),
)
