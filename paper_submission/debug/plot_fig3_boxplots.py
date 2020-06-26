import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
color_config = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        expressed_fbgns="../../output/cellselection-wf/expressed_genes.pkl",
        widely_expressed_fbgns="../../output/cellselection-wf/commonly_expressed_genes.pkl",
        tpm="../../output/seurat3-cluster-wf/tpm_by_cluster.feather",
        expressed_ratios="../../output/x-to-a-wf/db/expressed.dat",
        widely_expressed_ratios="../../output/x-to-a-wf/db/commonly_expressed.dat",
    ),
    params=dict(
        cluster_color=color_config["clusters"], cluster_order=config["cluster_order"]
    ),
)
