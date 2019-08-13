"""Save Zscores for Y genes.

Save a version of the zscore table with the Y genes pulled out.

1. Get list of Y gene FBgns and gene_symbols
2. Get zscores and pull out Y genes.
3. Map FBgns to gene symbol.
4. Sort by gene symbol

"""
import pandas as pd

from larval_gonad.io import feather_to_cluster_rep_matrix, shelve_dump


def main(snake):
    y_genes = get_y_genes(snake["gene_annot"])
    zscores = get_zscores(snake["zscores"], y_genes)
    shelve_dump(snake["output_file"], data=zscores)


def get_y_genes(input_file: str) -> pd.Series:
    """Create a series of Y genes."""
    return (
        pd.read_feather(input_file)
        .query("FB_chrom == 'Y'")
        .set_index("FBgn")
        .loc[:, "gene_symbol"]
        .squeeze()
    )


def get_zscores(input_file: str, target_genes: list) -> pd.DataFrame:
    """Import Zscores and pull out Y genes."""
    df = feather_to_cluster_rep_matrix(input_file).reindex(target_genes.index.tolist()).dropna()
    df.index = df.index.map(target_genes)
    return df.iloc[df.index.str.lower().argsort(), :]


if __name__ == "__main__":
    SNAKE = dict(
        zscores=snakemake.input.zscores,
        gene_annot=snakemake.input.gene_annot,
        output_file=snakemake.output[0],
    )

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
    #     print(os.getcwd())
    # except:
    #     pass
    # snake = dict(
    #     zscores="../../output/seurat3-cluster-wf/zscore_by_cluster_rep.feather",
    #     gene_annot="../../references/gene_annotation_dmel_r6-24.feather",
    #     output_file="",
    # )

    main(SNAKE)
