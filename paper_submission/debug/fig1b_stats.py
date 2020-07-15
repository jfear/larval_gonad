import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input="../../output/paper_submission/fig1b_sex_biased_expression.feather",
)