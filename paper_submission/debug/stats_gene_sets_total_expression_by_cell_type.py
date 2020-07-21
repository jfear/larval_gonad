import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        tpm="../../output/seurat3-cluster-wf/tpm_by_cluster.feather",
        expressed_fbgns="../../output/cellselection-wf/expressed_genes.pkl",
        widely_expressed_fbgns="../../output/cellselection-wf/commonly_expressed_genes.pkl",
        tau_fbgns="../../output/expression-atlas-wf/dmel_male_tau_fbgns.pkl",
        tsps_fbgns="../../output/expression-atlas-wf/dmel_male_tsps_fbgns.pkl",
    ),
    params=dict(
        order=config["cluster_order"]
    )
)