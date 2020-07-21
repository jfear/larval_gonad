import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
colors = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        expressed_fbgns="../../output/cellselection-wf/expressed_genes.pkl",
        widely_expressed_fbgns="../../output/cellselection-wf/commonly_expressed_genes.pkl",
        tau="../../output/expression-atlas-wf/dmel_tau.feather",
        tsps="../../output/expression-atlas-wf/dmel_tsps.feather",
        clusters="../../output/seurat3-cluster-wf/raw_by_cluster.feather"
    ),
    params=dict(
        order=config["cluster_order"],
        colors=colors["clusters"],
        names=config["cluster_names"]
    )
)
