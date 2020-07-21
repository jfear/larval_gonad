import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

config = read_config("../../config/common.yaml")

builtins.snakemake = MockSnake(
    input=dict(
        data="../../output/paper_submission/fig1b_sex_biased_expression.feather",
        stats="../../output/paper_submission/fig1b_stats.tsv",
    ),
    params=dict(chrom_order=config["chrom_order"]),
)
