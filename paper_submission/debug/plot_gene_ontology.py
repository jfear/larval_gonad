import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        background="../../output/cellselection-wf/expressed_genes.pkl",
        target_genes="../../output/cellselection-wf/commonly_expressed_genes.pkl",
    )
)
