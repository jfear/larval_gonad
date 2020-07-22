import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input="../../data/external/miriam/oligopaint_volumes.csv",
    params=dict(colors=["purple", "yellow"]),
)
