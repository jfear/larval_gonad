import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        bulk="../../output/bulk-rnaseq-wf/rnaseq_aggregation/gene_level_counts.tsv",
        sc="../../output/cellselection-wf/raw.feather",
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        male_biased="../../output/bulk-rnaseq-wf/deg/larval_testis_biased_fbgns.pkl",
    )
)
