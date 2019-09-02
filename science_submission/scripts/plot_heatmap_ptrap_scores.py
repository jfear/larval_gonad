import matplotlib

matplotlib.use("Agg")

import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

PTRAP_FILE = snakemake.input["ptrap"]
GENE_ANNOTATION = snakemake.input["gene_annot"]
OUTPUT_FILE = snakemake.output[0]

# Debug Settings
# import os
# try:
#     os.chdir(os.path.join(os.getcwd(), 'science_submission/scripts'))
#     print(os.getcwd())
# except:
#     pass
# PTRAP_FILE = '../../output/science_submission/ptrap_scores.feather'
# GENE_ANNOTATION = '../../references/gene_annotation_dmel_r6-26.feather'


def main():
    fbgn2symbol = (
        pd.read_feather(GENE_ANNOTATION, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .squeeze()
    )

    df = pd.read_feather(PTRAP_FILE).groupby("FBgn").max().rename(fbgn2symbol)

    # Cluster
    link = linkage(df.values, "average")
    tree = dendrogram(link, no_plot=True)
    leaves = tree["leaves"]
    df_sorted = df.iloc[leaves, :]

    # plot
    plt.style.use("scripts/figure_styles.mplstyle")
    fig = plt.figure(figsize=(4, 8))
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(
        df_sorted,
        xticklabels=True,
        yticklabels=True,
        vmin=0,
        vmax=3,
        rasterized=True,
        cmap="inferno",
        ax=ax,
        cbar_ax=cax,
        cbar_kws=dict(
            label="Protein Expression (Arbitrary Score)",
            ticks=[0, 1, 2, 3],
            orientation="horizontal",
        ),
        linewidths=.5,
    )

    # Clean up X axis
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_tick_params(pad=0, length=2)

    # Clean up Y axis
    ax.set_ylabel("")
    ax.yaxis.set_tick_params(pad=0.1, length=2)
    plt.setp(
        ax.get_yticklabels(), fontstyle="italic", fontfamily="Helvetica", rotation=0, va="center"
    )

    # Clean up color bar
    cax.xaxis.set_tick_params(pad=0, length=2)

    fig.savefig(OUTPUT_FILE, bbox_inches="tight")


if __name__ == "__main__":
    main()
