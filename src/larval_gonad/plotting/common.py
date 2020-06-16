import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_fbgn2chrom(gene_metadata: str) -> pd.Series:
    return (
        pd.read_feather(gene_metadata, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .squeeze()
        .rename("chrom")
    )


def get_fbgn2symbol(gene_metadata: str) -> pd.Series:
    return (
        pd.read_feather(gene_metadata, columns=["FBgn", "gene_symbol"]).set_index("FBgn").squeeze()
    )


def get_fbgn2length(gene_metadata: str) -> pd.Series:
    return pd.read_feather(gene_metadata, columns=["FBgn", "length"]).set_index("FBgn").squeeze()


def cluster_boxplot(y: str, data: pd.DataFrame, cluster_color: list, cluster_order: list, ax: plt.Axes, **kwargs):
    """Plot a boxplot with clusters on the x-axis.

    Parameters
    ----------
    y : str
        Value in data to use for y-axis.
    data : pd.DataFrame
        DataFrame containing at least two columns for x ("cluster") and y value.
    cluster_color : list
        List of colors to use
    cluster_order : list
        List ordering the x-axis
    ax : plt.Axes
    kwargs: dict
        Additional arguments passed to sns.boxplot.

    """
    args = dict(
        x="cluster",
        y=y,
        data=data,
        palette=cluster_color,
        order=cluster_order,
        notch=True,
        linewidth=.5,
        ax=ax
    )
    args.update(kwargs)
    sns.boxplot(**args)
    ax.set(xlabel="")
    return ax
