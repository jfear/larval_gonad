import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input="../../output/paper_submission/gsea_by_cell_type.feather",
)
