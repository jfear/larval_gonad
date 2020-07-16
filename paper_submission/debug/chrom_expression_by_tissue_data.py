import builtins

from larval_gonad.mock import MockSnake
from larval_gonad.config import read_config

builtins.snakemake = MockSnake(
    input=dict(
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        w1118="../../output/expression-atlas-wf/w1118_gene_counts.feather",
        OreR="../../output/expression-atlas-wf/OreR_gene_counts.feather",
        L3_bulk="../../output/bulk2-rnaseq-wf/testis_ovary_counts.feather",
    ),
)
