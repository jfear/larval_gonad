import builtins

from larval_gonad.mock import MockSnake

builtins. snakemake = MockSnake(
    input="../../data/external/camila/fourth_distance_20200228.tsv",
    params=dict(colors=["red", "cyan"]),
)
