"""Save Zscores for Literature genes.

Save a version of the zscore table with the Literature genes pulled out.

1. Get list of literature gene FBgns and gene_symbols
2. Get zscores and pull out literature genes.
3. Map FBgns to gene symbol.
4. Get blocks to seperate lit genes for each cell type.

"""
from more_itertools import flatten
import numpy as np
import pandas as pd

from larval_gonad.io import feather_to_cluster_rep_matrix, shelve_dump


def main(snake):
    lit_genes = get_lit_genes(snake["gene_annot"], list(flatten(snake["lit_genes"].values())))
    zscores = get_zscores(snake["zscores"], lit_genes)
    blocks = get_blocks(snake["lit_genes"])
    shelve_dump(snake["output_file"], data=zscores, blocks=blocks)


def get_lit_genes(input_file: str, fbgns: list) -> pd.Series:
    """Create a series of Y genes."""
    return (
        pd.read_feather(input_file).set_index("FBgn").reindex(fbgns).loc[:, "gene_symbol"].squeeze()
    )


def get_zscores(input_file: str, target_genes: list) -> pd.DataFrame:
    """Import Zscores and pull out Y genes."""
    df = feather_to_cluster_rep_matrix(input_file).reindex(target_genes.index.tolist()).dropna()
    df.index = df.index.map(target_genes)
    return df


def get_blocks(lit_genes: dict) -> list:
    """Count the number of genes.

    Count the number of genes describing each cell type. Then take the
    cumulative sum to get locations to draw lines on a heatmap.

    """
    return np.cumsum([len(x) for x in lit_genes.values()])


if __name__ == "__main__":
    SNAKE = dict(
        zscores=snakemake.input.zscores,
        gene_annot=snakemake.input.gene_annot,
        output_file=snakemake.output[0],
        lit_genes=snakemake.params[0],
    )

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
    #     print(os.getcwd())
    # except:
    #     pass
    # from larval_gonad.config import read_config
    # snake = dict(
    #     zscores="../../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
    #     gene_annot="../../references/gene_annotation_dmel_r6-24.feather",
    #     output_file="",
    #     lit_genes=read_config("../../config/literature_genes.yaml")
    # )

    main(SNAKE)
