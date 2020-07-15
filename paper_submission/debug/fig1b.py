import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

builtins.snakemake = MockSnake(
    input=dict(
        data="../../output/paper_submission/fig1b_sex_biased_expression.feather",
        stats="../../output/paper_submission/fig1b_stats.tsv",
    ),
    params=dict(colors=read_config("../../config/colors.yaml")),
)
