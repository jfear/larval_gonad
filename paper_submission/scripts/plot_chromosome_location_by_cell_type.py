import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import larval_gonad.plotting  # pylint: disable=unused-import
from larval_gonad.normalization import zscore

plt.style.use("minimal")


SCATTER_KWS = dict(s=2, rasterized=True)
FACET_KWS = dict(
    row="cluster",
    col="chrom",
    hue="cluster",
    row_order=snakemake.params.cluster_order,
    col_order=snakemake.params.chrom_order,
    palette=snakemake.params.cluster_colors,
    sharey=True,
    sharex=True,
    height=0.8,
    aspect=4,
    margin_titles=True,
)


def main():
    gene_annot = load_gene_annot()
    df = load_tpm().pipe(zscore).pipe(melt).merge(gene_annot)

    g = sns.FacetGrid(df, **FACET_KWS)
    # Plot expressed genes along the chormosome (z-score)
    g.map(plot, "start", "z-score", **SCATTER_KWS)
    # Plot non-expressed genes along the chormosome (at -2)
    g.map(plot_missing_genes, "chrom", "FBgn", data=gene_annot, y_loc=-2, **SCATTER_KWS)

    g.set_titles(row_template="{row_name}", col_template="{col_name}")
    g.set_xlabels("")
    g.set_ylabels("Z-Score (TPM)")
    plt.subplots_adjust(hspace=0.05, wspace=0)

    g.savefig(snakemake.output[0])


def plot(x, y, **kwargs):
    _ = kwargs.pop("label")
    plt.scatter(x, y, **kwargs)
    ax = plt.gca()
    ax.axhline(0, color="k", ls="--")
    return ax


def plot_missing_genes(chroms, fbgns, **kwargs):
    _chrom = chroms.values[0]
    _fbgns = fbgns.values.tolist()

    data = kwargs.pop("data").copy()
    data["y"] = kwargs.pop("y_loc")
    df = data.query("FBgn not in @_fbgns and chrom == @_chrom")

    _ = kwargs.pop("label")
    _ = kwargs.pop("color")

    return plt.scatter(df.start, df.y, color="r", marker="^", **kwargs)


def load_tpm() -> pd.DataFrame:
    """Load TMP as a wide data set"""
    return pd.pivot(
        pd.read_feather(snakemake.input.tpm),
        index="FBgn",
        columns="cluster",
        values="TPM",
    )


def load_gene_annot() -> pd.DataFrame:
    return (
        pd.read_feather(
            snakemake.input.gene_annot, columns=["FBgn", "FB_chrom", "start"]
        )
        .rename(columns={"FB_chrom": "chrom"})
        .assign(
            chrom=lambda x: pd.Categorical(
                x.chrom, categories=snakemake.params.chrom_order, ordered=True
            )
        )
    )


def melt(df: pd.DataFrame) -> pd.DataFrame:
    return df.stack().rename("z-score").reset_index()


if __name__ == "__main__":
    main()
