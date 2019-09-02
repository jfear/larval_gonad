"""Save zscores for houskeeping genes (tau).

Save a version of the zscore table with the housekeeping genes pulled out.

1. Get list of housekeeping FBgns and gene_symbols
2. Get zscores and pull out houskeeping genes.
3. Map FBgns to gene symbol.
4. Group genes by hierarchical clustering.

"""
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram

from larval_gonad.io import feather_to_cluster_rep_matrix, shelve_dump, pickle_load


def main(snake):
    fbgns = get_genes(snake['gene_annot'], pickle_load(snake['fbgns']))
    zscores = get_zscores(snake["zscores"], fbgns)
    zscores_sorted = hierarchical_cluster_zscores(zscores)

    shelve_dump(snake["output_file"], data=zscores_sorted)


def get_genes(input_file: str, fbgns: list) -> pd.Series:
    """Create a series of target genes."""
    return (
        pd.read_feather(input_file).set_index("FBgn").reindex(fbgns).loc[:, "gene_symbol"].squeeze()
    )


def get_zscores(input_file: str, target_genes: list) -> pd.DataFrame:
    """Import zscores and pull out target genes."""
    df = feather_to_cluster_rep_matrix(input_file).reindex(target_genes.index.tolist()).dropna()
    df.index = df.index.map(target_genes)
    return df


def hierarchical_cluster_zscores(df: pd.DataFrame) -> pd.DataFrame:
    links = linkage(df, method='average')
    tree = dendrogram(links, no_plot=True)
    leaves = tree['leaves']
    return df.iloc[leaves, :]


if __name__ == "__main__":
    SNAKE = dict(
        zscores=snakemake.input.zscores,
        gene_annot=snakemake.input.gene_annot,
        fbgns=snakemake.input.fbgns,
        output_file=snakemake.output[0],
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
    #     gene_annot="../../references/gene_annotation_dmel_r6-26.feather",
    #     fbgns="../../output/expression-atlas-wf/dmel_male_tau_fbgns.pkl",
    #     output_file="",
    # )

    main(SNAKE)
