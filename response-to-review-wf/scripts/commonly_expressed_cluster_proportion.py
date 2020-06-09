"""Summary of the commonly expressed genes (>= 33% cells) """

import joblib
import pandas as pd


def _debug():
    from larval_gonad.mock import MockSnake

    snakemake = MockSnake(
        input=dict(
            annotation="../../references/gene_annotation_dmel_r6-26.feather",
            commonly_expressed="../../output/cellselection-wf/commonly_expressed_genes.pkl",
            clusters="../../output//seurat3-cluster-wf/combined_n3_clusters.feather",
            raw="../../output/cellselection-wf/raw.feather",
            tau="../../output/expression-atlas-wf/dmel_tau.feather",
            tsps="../../output/expression-atlas-wf/dmel_tsps.feather",
        )
    )


def main():
    expression = commonly_expressed_by_cluster()
    prop_cells_from_cluster = cluster_proportion(expression)
    fbgn_to_tau_tsps = load_tau_tsps()
    prop_cells_from_cluster.join(fbgn_to_tau_tsps).to_csv(snakemake.output[0], sep="\t")


def commonly_expressed_by_cluster() -> pd.DataFrame:
    annotation = (
        pd.read_feather(snakemake.input.annotation)
        .set_index("FBgn")
        .gene_symbol.to_frame()
        .reset_index()
    )

    common_fbgns = joblib.load(snakemake.input.commonly_expressed)
    clusters = pd.read_feather(snakemake.input.clusters)[["cell_id", "cluster"]]

    raw = (
        pd.read_feather(snakemake.input.raw)
        .query("FBgn in @common_fbgns")
        .melt(id_vars="FBgn", var_name="cell_id", value_name="UMI")
        .merge(clusters)
        .merge(annotation)
    )

    return raw


def cluster_proportion(expression: pd.DataFrame) -> pd.DataFrame:
    num_cells_with_expression_by_cluster = (
        expression.query("UMI > 0")
        .groupby(["FBgn", "gene_symbol"])
        .cluster.value_counts()
        .unstack()
        .loc[:, expression.cluster.cat.categories]
        .assign(total=lambda df: df.sum(axis=1))
        .set_index("total", append=True)
    )

    return num_cells_with_expression_by_cluster.div(
        num_cells_with_expression_by_cluster.index.get_level_values("total"), axis=0
    )


def load_tau_tsps() -> pd.DataFrame:
    tissue_specificity_scores = pd.concat(
        [
            pd.read_feather(snakemake.input.tau).set_index("FBgn").male_tau,
            pd.read_feather(snakemake.input.tsps).set_index("FBgn").male_tsps,
        ],
        axis=1,
    )

    return tissue_specificity_scores


if __name__ == "__main__":
    main()
