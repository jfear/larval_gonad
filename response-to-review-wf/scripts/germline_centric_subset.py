import pandas as pd
import joblib

"""Reviewer 1 proposed:

A better approach may be to systematically analyze genes that are expressed
in both the germ cells of interest and one or more somatic cell types

"""

GERMLINE = ["G", "EPS", "MPS", "LPS"]
SOMA = ["C1", "C2", "C3", "C4", "T", "P"]


def _debug():
    from larval_gonad.mock import MockSnake

    snakemake = MockSnake(
        input=dict(
            raw="../../output/cellselection-wf/raw.feather",
            clusters="../../output/seurat3-cluster-wf/combined_n3_clusters.feather",
        )
    )


def main():
    clusters = pd.read_feather(snakemake.input.clusters)
    prop_cells_on_by_cluster = (
        read_raw(snakemake.input.raw, clusters.cell_id.to_list())
        .merge(clusters, on="cell_id")
        .groupby(["FBgn", "cluster"])
        .gene_on.mean()
    )

    germline_fbgns = genes_expressed_in_at_least_two_germline(prop_cells_on_by_cluster)
    soma_fbgns = genes_expressed_in_at_least_one_soma(prop_cells_on_by_cluster)
    target_fbgns = list(germline_fbgns.intersection(soma_fbgns))

    joblib.dump(target_fbgns, snakemake.output[0])


def read_raw(file_name: str, cells: list) -> pd.DataFrame:
    return (
        pd.read_feather(file_name, columns=["FBgn"] + cells)
        .set_index("FBgn")
        .pipe(drop_not_expressed_genes)
        .reset_index("FBgn")
        .melt(id_vars="FBgn", var_name="cell_id", value_name="UMI")
        .assign(gene_on=lambda df: df.UMI > 0)
        .drop("UMI", axis=1)
    )


def drop_not_expressed_genes(df: pd.DataFrame) -> pd.DataFrame:
    mask = df.sum(axis=1) > 0
    return df[mask]


def genes_expressed_in_at_least_two_germline(
    prop_cells_on_by_cluster: pd.DataFrame
) -> set:
    genes_expressed_in_10pct_cells_by_cluster = (
        prop_cells_on_by_cluster > 0.1
    ).to_frame()
    num_germline_clusters_with_expression = (
        genes_expressed_in_10pct_cells_by_cluster.query("cluster in @GERMLINE")
        .groupby("FBgn")
        .sum()
        .squeeze()
    )

    mask = num_germline_clusters_with_expression >= 2
    return set(num_germline_clusters_with_expression[mask].index.to_list())


def genes_expressed_in_at_least_one_soma(prop_cells_on_by_cluster: pd.DataFrame) -> set:
    genes_expressed_in_10pct_cells_by_cluster = (
        prop_cells_on_by_cluster > 0.1
    ).to_frame()
    num_soma_clusters_with_expression = (
        genes_expressed_in_10pct_cells_by_cluster.query("cluster in @SOMA")
        .groupby("FBgn")
        .sum()
        .squeeze()
    )

    mask = num_soma_clusters_with_expression >= 1
    return set(num_soma_clusters_with_expression[mask].index.to_list())


if __name__ == "__main__":
    main()
