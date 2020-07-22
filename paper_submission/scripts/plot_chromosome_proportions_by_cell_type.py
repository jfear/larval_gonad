"""Plot Chromosome Expression By Cell Type.

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import larval_gonad.plotting  # pylint: disable=unused-import

plt.style.use("minimal")

COLORS = {
    k: v
    for k, v in zip(snakemake.params.cluster_order, snakemake.params.cluster_colors)
}


def main():
    df = pd.concat(
        [
            load_data(snakemake.input.all_genes, "All Expressed"),
            load_data(snakemake.input.tau_genes, "Tau"),
            load_data(snakemake.input.tsps_genes, "TSPS"),
            load_data(snakemake.input.widely_expressed_genes, "Widely Expressed"),
        ],
        sort=False,
        ignore_index=True,
    ).merge(pd.read_feather(snakemake.input.clusters), on="cell_id")

    g = sns.FacetGrid(
        df,
        row="gene_set",
        row_order=["All Expressed", "Tau", "TSPS", "Widely Expressed"],
        col="cluster",
        col_order=snakemake.params.cluster_order,
        hue="cluster",
        sharex=True,
        sharey="row",
    )
    g.map(
        boxplot,
        "chrom",
        "scaled_prop_reads",
        showfliers=False,
        notch=True,
        order=snakemake.params.chrom_order,
    )
    g.set_titles("{row_name}:{col_name}")
    g.set_xlabels("")
    g.set_ylabels("")
    tweak(g)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    g.savefig(snakemake.output[0])


def load_data(file_name: str, name: str):
    return pd.read_feather(file_name).assign(gene_set=name)


def boxplot(x, y, **kwargs):
    _ = kwargs.pop("color")
    label = kwargs.pop("label")
    sns.boxplot(x, y, color=COLORS[label], **kwargs)


def tweak(g: sns.FacetGrid):
    # Label Y-axis based on row
    for ax in g.axes[:, 0]:
        row, _ = ax.get_title().split(":")
        if row == "Tau":
            row = r"$\tau$"
        ax.set_ylabel(f"{row}\nScaled Proportion Reads")

    # Change titles on top row
    for ax in g.axes[0, :]:
        _, col = ax.get_title().split(":")
        title = snakemake.params.cluster_names[col]
        ax.set_title(title, fontdict=dict(fontsize=14, fontweight="bold"))

    # Remove titles everywhere else
    for ax in g.axes[1:, :].flat:
        ax.set_title("")

    return g


if __name__ == "__main__":
    main()
