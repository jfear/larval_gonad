import builtins

from snakemake.io import expand

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        deg=expand(
            "../../output/expression-atlas-wf/sex_biased_expression/{species}_{tissue}.tsv",
            species=["w1118", "orgR"],
            tissue=["AC", "DG", "GE", "GO", "HD", "RE", "TX", "WB"],
        ),
        yogn2fbgn="../../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dmel.pkl"
    )
)
