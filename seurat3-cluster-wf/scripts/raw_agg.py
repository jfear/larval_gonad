import pandas as pd

from larval_gonad.io import melt_cluster_rep_matrix

RAW = snakemake.input["raw"]
METADATA = snakemake.input["metadata"]
CLUSTER_ANNOT = snakemake.params["cluster_annot"]
CLUSTER_ORDER = snakemake.params["cluster_order"]

OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# RAW = "output/science_submission/raw.feather"
# METADATA = "output/seurat3-cluster-wf/combined_n3_metadata.feather"
# import yaml
# config = yaml.safe_load(open("config/common.yaml"))
# CLUSTER_ANNOT = config["cluster_annot"]
# CLUSTER_ORDER = config["cluster_order"]


def main():
    cell_annot = (
        pd.read_feather(METADATA, columns=["cell_id", "cluster"])
        .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
        .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
        .assign(rep=lambda df: df.cell_id.str.extract("(rep\d)", expand=False))
        .assign(
            rep=lambda df: pd.Categorical(df.rep, categories=sorted(df.rep.unique()))
        )
        .set_index("cell_id")
    )

    # remove genes that are in fewer than 3 cells
    raw = (
        pd.read_feather(RAW)
        .set_index("FBgn")
        .rename_axis(columns="cell_id")
        .pipe(lambda x: x[(x > 0).sum(axis=1) >= 3])
    )

    # Aggregate by cluster and rep
    agg = raw.T.join(cell_annot, how="right").groupby(["cluster", "rep"]).sum().T

    # Save
    melt_cluster_rep_matrix(agg, name="UMI").to_feather(OUTPUT_FILE)


if __name__ == "__main__":
    main()
