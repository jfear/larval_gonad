"""Plot Chromosome Expression By Cell Type.

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import larval_gonad.plotting  # pylint: disable=unused-import

plt.style.use("minimal")

BOXPLOT_KWS = dict(showfliers=False, notch=True, order=snakemake.params.chrom_order)
FACET_KWS = dict(
    row="gene_set",
    col="cluster",
    hue="cluster",
    palette=snakemake.params.cluster_colors,
    row_order=["All Expressed", "Tau", "TSPS", "Widely Expressed"],
    col_order=snakemake.params.cluster_order,
    sharex=True,
    sharey="row",
    margin_titles=True,
    height=1.5
)


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

    g = sns.FacetGrid(df, **FACET_KWS)
    g.map(boxplot, "chrom", "scaled_prop_reads", **BOXPLOT_KWS)

    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    g.set_xlabels("")
    g.set_ylabels("Proportion Reads\n# Expressed Genes")
    g.despine(left=True)
    plt.subplots_adjust(hspace=0, wspace=0.05)

    g.savefig(snakemake.output[0])


def load_data(file_name: str, name: str):
    return pd.read_feather(file_name).assign(gene_set=name)


def boxplot(x, y, **kwargs):
    _ = kwargs.pop("label")
    sns.boxplot(x, y, **kwargs)


if __name__ == "__main__":
    main()
