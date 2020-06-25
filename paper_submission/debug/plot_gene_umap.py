import builtins

from larval_gonad.mock import MockSnake

builtins.snakemake = MockSnake(
    input=dict(
        gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
        zscores="../../output/seurat3-cluster-wf/zscore_by_cell.feather",
        umap="../../output/seurat3-cluster-wf/combined_n3_umap.feather",
    ),
)
