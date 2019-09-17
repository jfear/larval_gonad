import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_umap_panel(umap_data: str, cluster_data: str, colors: list, axes: list = None):
    """Plot panel of UMAP for each replicate.

    Example
    -------
    >>> umap_data = "output/seurat3-cluster-wf/combined_n3_umap.feather"
    >>> cluster_data = "output/seurat3-cluster-wf/combined_n3_clusters.feather"
    >>> from larval_gonad.config import read_config
    >>> color_config = read_config("config/colors.yaml")
    >>> colors = color_config["clusters"]
    """
    umap = pd.read_feather(umap_data).set_index("cell_id")
    clusters = pd.read_feather(cluster_data).set_index("cell_id")
    df = umap.join(clusters)

    if axes is None:
        _, axes = plt.subplots(1, 4, figsize=plt.figaspect(1 / 4), sharex=True, sharey=True)

    _make_scatter_panel(df, axes, colors)
    _add_labels(df, axes[-1])
    _cleanup_axes(axes)
    _add_axes_labels(plt.gcf())


def _make_scatter_panel(df, axes, colors):
    defaults = dict(
        x="UMAP_1",
        y="UMAP_2",
        hue="cluster",
        palette=colors,
        s=2,
        linewidth=0.02,
        edgecolor="k",
        rasterized=True,
        legend=False,
    )

    # Plot each replicate
    for (rep, subset), ax in zip(df.groupby("rep"), axes.flat[:-1]):
        sns.scatterplot(data=subset, ax=ax, **defaults)
        ax.set(title=rep)

    # Plot combined
    ax_combined = axes.flat[-1]
    sns.scatterplot(data=df, ax=ax_combined, **defaults)
    ax_combined.set(title="Combined")


def _add_labels(df, ax):
    """Add cluster labels to the combined panel"""
    for clus, row in df.groupby("cluster")[["UMAP_1", "UMAP_2"]].mean().iterrows():
        ax.text(
            row.UMAP_1,
            row.UMAP_2,
            clus,
            bbox=dict(facecolor=(1, 1, 1, 0.8), edgecolor="none", pad=0.2),
            ha="center",
            va="center",
            fontweight="bold",
        )


def _cleanup_axes(axes):
    for ax in axes:
        sns.despine(ax=ax)
        ax.set(xticks=[-10, 10], yticks=[-10, 10], xlabel="", ylabel="", aspect="equal")


def _add_axes_labels(fig):
    defaults = dict(transform=fig.transFigure, ha="center", va="center")
    fig.text(0.5, 0.05, "Profile Similarity (UMAP-1)", **defaults)
    fig.text(0.1, 0.5, "Profile Similarity (UMAP-2)", rotation=90, **defaults)
