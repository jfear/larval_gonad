import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        raw="../../output/cellselection-wf/raw.feather",
        gene_list="../../output/cellselection-wf/commonly_expressed_genes.pkl",
        gene_annotation="../../references/gene_annotation_dmel_r6-26.feather",
    )
)
