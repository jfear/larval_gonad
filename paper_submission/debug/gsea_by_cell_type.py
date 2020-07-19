import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        sc="../../output/cellselection-wf/raw.feather",
        clusters="../../output/seurat3-cluster-wf/combined_n3_clusters.feather",
        expressed="../../output/cellselection-wf/expressed_genes.pkl",
        male_biased="../../output/bulk-rnaseq-wf/deg/larval_testis_biased_fbgns.pkl",
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
    ),
    threads=20,
)
