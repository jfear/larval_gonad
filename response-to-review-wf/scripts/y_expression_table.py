import pandas as pd


def _debug():
    from larval_gonad.mock import MockSnake

    snakemake = MockSnake(
        input=dict(
            annotation="../../references/gene_annotation_dmel_r6-26.feather",
            expressed="../../output/cellselection-wf/expressed_genes.pkl",
            widely_expressed="../../output/cellselection-wf/commonly_expressed_genes.pkl",
            raw="../../output/cellselection-wf/raw.feather",
            clusters="../../output/seurat3-cluster-wf/combined_n3_clusters.feather",
        )
    )


def main():
    y_genes = get_y_genes()
    expression = load_counts(y_genes.index.to_list())
    expression.index = pd.MultiIndex.from_frame(
        y_genes.reindex(expression.index).reset_index()
    )
    expression.to_csv(snakemake.output[0], sep="\t")


def get_y_genes():
    return (
        pd.read_feather(snakemake.input.annotation)
        .query("FB_chrom == 'Y'")
        .set_index("FBgn")
        .gene_symbol
    )


def load_counts(y_fbgns):
    clusters = pd.read_feather(snakemake.input.clusters)
    return (
        pd.read_feather(snakemake.input.raw)
        .query("FBgn in @y_fbgns")
        .melt(id_vars="FBgn", var_name="cell_id", value_name="UMI")
        .merge(clusters)
        .groupby(["FBgn", "cluster"])
        .UMI.sum()
        .unstack()
        .pipe(drop_not_expressed)
    )


def drop_not_expressed(df: pd.DataFrame) -> pd.DataFrame:
    mask = df.sum(axis=1) > 0
    return df[mask]


if __name__ == "__main__":
    main()
