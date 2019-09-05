"""Create UMAP panel of individual reps and combined"""
import os
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns


def main():
    # Get UMAP coords, cluster, and rep info
    umap = pd.read_feather(snakemake.input.umap).set_index("cell_id")
    clusters = pd.read_feather(snakemake.input.clusters).set_index("cell_id")
    df = umap.join(clusters)

    # Create figure panel with each rep separate and combined
    n_panels = len(df.rep.unique()) + 1
    fig, axes = plt.subplots(
        1,
        n_panels,
        figsize=plt.figaspect(1 / n_panels),
        sharex=True,
        sharey=True,
        gridspec_kw=dict(wspace=0.01),
    )
    plot_scatter(df, axes)
    plot_labels(df, axes.flat[-1])

    plt.savefig(snakemake.output[0])


def plot_scatter(df, axes):
    defaults = dict(
        x="UMAP_1",
        y="UMAP_2",
        hue="cluster",
        palette=snakemake.params.colors,
        s=3,
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


def plot_labels(df, ax_combined):
    """Add cluster labels to the combined panel"""
    for clus, row in df.groupby("cluster")[["UMAP_1", "UMAP_2"]].mean().iterrows():
        ax_combined.text(
            row.UMAP_1,
            row.UMAP_2,
            clus,
            bbox=dict(facecolor=(1, 1, 1, 0.8), edgecolor="none", pad=0.2),
            ha="center",
            va="center",
            fontweight="bold",
        )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
            ),
            params=dict(colors=sns.color_palette(n_colors=10)),
        )
        plt.style.use("../config/figure_styles.mplstyle")

    main()
