import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        deg="../../output/seurat3-cluster-wf/germline_deg/GvLPS.feather",
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
    ),
    params=dict(
        chrom_order=config["chrom_order"]
    )
)
