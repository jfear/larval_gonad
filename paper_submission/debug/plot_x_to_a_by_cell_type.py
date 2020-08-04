import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
color_config = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        tau="../../output/x-to-a-wf/autosome_ratios_tau_by_cell.feather",
        tsps="../../output/x-to-a-wf/autosome_ratios_tsps_by_cell.feather",
        widely_expressed="../../output/x-to-a-wf/autosome_ratios_commonly_expressed_by_cell.feather",
        tau_pvals="../../output/x-to-a-wf/db/tau.dat",
        tsps_pvals="../../output/x-to-a-wf/db/tsps.dat",
        widely_expressed_pvals="../../output/x-to-a-wf/db/commonly_expressed.dat",
    ),
    params=dict(
        colors=color_config["clusters"],
        order=config["cluster_order"],
        names=config["cluster_names"],
    ),
)
