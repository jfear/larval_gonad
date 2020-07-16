import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        adult="../../output/expression-atlas-wf/dmel_sex_biased_expression.feather",
        larval="../../output/bulk-rnaseq-wf/deg/bulk_testis_vs_ovary.tsv",
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
    )
)
