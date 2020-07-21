import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")
colors = read_config("../../config/colors.yaml")

builtins.snakemake = MockSnake(
    input="../../output/paper_submission/gsea_by_cell_type.feather",
    params=dict(
        colors=colors["clusters"],
        order=config["cluster_order"],
        names=config["cluster_names"],
    ),
)
