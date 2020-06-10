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
    prop_cells_on = proportion_cells_with_expression(expression)
    fbgn2symbol = read_annotation()
    fbgn_to_tau_tsps = load_tau_tsps()

    prop_cells_on.join(fbgn2symbol).join(fbgn_to_tau_tsps).set_index(
        "gene_symbol", append=True
    ).to_csv(snakemake.output[0], sep="\t")


def read_annotation() -> pd.Series:
    return pd.read_feather(snakemake.input.annotation).set_index("FBgn").gene_symbol


def commonly_expressed_by_cluster() -> pd.DataFrame:
    common_fbgns = joblib.load(snakemake.input.commonly_expressed)
    clusters = pd.read_feather(snakemake.input.clusters)[["cell_id", "cluster"]]

    raw = (
        pd.read_feather(snakemake.input.raw)
        .query("FBgn in @common_fbgns")
        .melt(id_vars="FBgn", var_name="cell_id", value_name="UMI")
        .merge(clusters)
    )

    return raw


def proportion_cells_with_expression(expression: pd.DataFrame) -> pd.DataFrame:
    # Total Cell Counts for colnames
    total_cells = (
        expression[["cell_id", "cluster"]]
        .drop_duplicates()
        .groupby("cluster")
        .size()
        .rename("num_cells")
        .map(lambda cell: f"(n={cell:,})")
        .reset_index()
    )
    column_names = pd.MultiIndex.from_frame(total_cells)

    # Calculate Prop Cells On (0 < UMI)
    prop_on = (
        expression.assign(gene_on=lambda df: df.UMI > 0)
        .groupby(["FBgn", "cluster"])
        .gene_on.mean()
        .rename("prop_on")
        .unstack()
    )
    prop_on.columns = column_names

    return prop_on


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
