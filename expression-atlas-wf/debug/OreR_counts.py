import builtins
from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        orthologs="../../output/expression-atlas-wf/YOgn_to_dmel_ortholog/dmel.pkl",
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        counts="../../output/expression-atlas-wf/aggregated_counts_table/orgR.feather",
    ),
)
