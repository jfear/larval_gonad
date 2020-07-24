"""Plot heatmap of X, 4, and Y escapers."""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

import larval_gonad.plotting
from larval_gonad.plotting.biomarkers import _cleanup_xaxis

plt.style.use("minimal")

GERM_LINE = ["G", "EPS", "MPS", "LPS"]
SEX_CHROM = ["X", "Y", "4"]


def main():
    gene_annot = load_gene_annot()

    zscores = (
        pd.read_feather(snakemake.input.zscores)
        .merge(gene_annot)
        .query("chrom in @SEX_CHROM")
        .query("cluster in @GERM_LINE")
        .pipe(add_LPS_G_difference)
        .set_index(["FBgn", "gene_symbol", "chrom", "difference", "rep", "cluster"])
        .squeeze()
        .unstack()
        .unstack()
        .sort_index(level=["chrom", "difference"])
    )

    # Output table of escapers
    zscores.query("difference > 0").index.to_frame().reset_index(drop=True).to_csv(
        snakemake.output.escapers_table, sep="\t", index=False
    )

    # Output Summary Table
    pd.concat(
        [
            zscores.groupby("chrom").size().rename("num_genes"),
            zscores.groupby("chrom")
            .apply(lambda df: (df.index.get_level_values("difference") > 0).sum())
            .rename("num_escapers"),
        ],
        axis=1,
    ).assign(prop_escapers=lambda x: x.num_escapers / x.num_genes).to_csv(
        snakemake.output.escapers_summary, sep="\t"
    )

    # Plot all genes
    fig = plt.figure(figsize=plt.figaspect(2))
    plot_all_genes(zscores, fig)
    plt.savefig(snakemake.output.all_genes)
    plt.close()

    # Plot escapers
    _df = zscores.query("difference > 0").groupby(["chrom"]).tail(40)
    _df.index = _df.index.droplevel(["FBgn", "difference"])
    fig = plt.figure(figsize=(8, 8))
    plot_all_genes(_df, fig, wspace=1.4, yticklabels=True)
    plt.savefig(snakemake.output.escapers)
    plt.close()


def load_gene_annot():
    return (
        pd.read_feather(
            snakemake.input.gene_annot, columns=["FBgn", "gene_symbol", "FB_chrom"]
        )
        .rename(columns={"FB_chrom": "chrom"})
        .assign(
            chrom=lambda df: pd.Categorical(
                df.chrom, categories=snakemake.params.chrom_order, ordered=True
            )
        )
    )


def add_LPS_G_difference(df: pd.DataFrame) -> pd.DataFrame:
    differences = (
        pd.read_feather(snakemake.input.raw)
        .query("cluster in @GERM_LINE")
        .set_index(["FBgn", "cluster"])
        .squeeze()
        .unstack()
        .pipe(lambda x: x.LPS - x.G)
        .rename("difference")
        .reset_index()
    )
    return df.merge(differences)


def plot_all_genes(
    df: pd.DataFrame, fig: plt.Figure, wspace=0.1, **kwargs
) -> pd.DataFrame:
    heatmap_kws = dict(
        vmin=-3,
        vmax=3,
        cmap=snakemake.params.color,
        yticklabels=False,
        xticklabels=True,
        rasterized=True,
    )
    heatmap_kws.update(kwargs)

    # Build figure layout
    gs = GridSpec(2, 3, height_ratios=[1, 0.01], hspace=0.01, wspace=wspace)
    ax_X = fig.add_subplot(gs[0, 0])
    ax_Y = fig.add_subplot(gs[0, 1])
    ax_4 = fig.add_subplot(gs[0, 2])
    cax = fig.add_subplot(gs[1, :])

    # Add heatmaps
    sns.heatmap(
        df.query("chrom == 'X'"),
        ax=ax_X,
        cbar_ax=cax,
        cbar_kws=dict(orientation="horizontal"),
        **heatmap_kws
    )
    sns.heatmap(df.query("chrom == 'Y'"), ax=ax_Y, cbar=False, **heatmap_kws)
    sns.heatmap(df.query("chrom == '4'"), ax=ax_4, cbar=False, **heatmap_kws)

    # Tweak Axes
    _cleanup_xaxis(ax_X, ["G", "E1°", "M1°", "L1°"])
    _cleanup_xaxis(ax_Y, ["G", "E1°", "M1°", "L1°"])
    _cleanup_xaxis(ax_4, ["G", "E1°", "M1°", "L1°"])

    _cleanup_yaxis(ax_X)
    _cleanup_yaxis(ax_Y)
    _cleanup_yaxis(ax_4)

    ax_X.set_title("X")
    ax_Y.set_title("Y")
    ax_4.set_title("4")


def _cleanup_yaxis(ax: plt.Axes):
    ax.set(ylabel="")


if __name__ == "__main__":
    main()
