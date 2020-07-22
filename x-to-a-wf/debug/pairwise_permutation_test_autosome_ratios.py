import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")

builtins.snakemake = MockSnake(
    input="../../output/x-to-a-wf/autosome_ratios_expressed_by_cell.feather",
    params=dict(
        cluster_order=config["cluster_order"]
    ),
    threads=20
)
