import pandas as pd

from larval_gonad.io import melt_cluster_matrix


def main():
    cell_annot = pd.read_feather(snakemake.input.clusters).set_index("cell_id")

    # remove genes that are in fewer than 3 cells
    raw = (
        pd.read_feather(snakemake.input.raw)
        .set_index("FBgn")
        .rename_axis(columns="cell_id")
        .pipe(lambda x: x[(x > 0).sum(axis=1) >= 3])
    )

    # Aggregate by cluster and rep
    agg = raw.T.join(cell_annot, how="right").groupby("cluster").sum().T

    # Save
    melt_cluster_matrix(agg, name="UMI").to_feather(snakemake.output[0])


if __name__ == "__main__":
    main()
