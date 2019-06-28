import matplotlib

matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

FILE_NAME = snakemake.input[0]

CLUSTER_ANNOT = snakemake.params.cluster_annot
CLUSTER_ORDER = snakemake.params.cluster_order

ONAME = snakemake.output[0]

# Debug Settings
# FILE_NAME = 'output/seurat3-cluster-wf/combined_n3_metadata.feather'
# import yaml
# config = yaml.safe_load(open('config/common.yaml'))
# CLUSTER_ANNOT = config['cluster_annot']
# CLUSTER_ORDER = config['cluster_order']


def main():
    df = (
        pd.read_feather(FILE_NAME, columns=["cell_id", "cluster"])
        .assign(cluster=lambda df: df.cluster.cat.rename_categories(CLUSTER_ANNOT))
        .assign(cluster=lambda df: df.cluster.cat.reorder_categories(CLUSTER_ORDER))
        .assign(rep=lambda df: df.cell_id.str.extract("(rep\d)", expand=False))
        .groupby(["cluster", "rep"])
        .size()
        .rename("num_cells")
        .reset_index()
    )

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, ax = plt.subplots()

    sns.barplot("cluster", "num_cells", hue="rep", data=df, ax=ax)
    ax.set(xlabel="", ylabel="Number of Cells")
    ax.set_axisbelow(True)
    ax.grid(axis="y")

    fig.savefig(ONAME, bbox_inches="tight")


if __name__ == "__main__":
    main()
