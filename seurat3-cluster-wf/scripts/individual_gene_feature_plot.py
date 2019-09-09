import os
from pathlib import Path
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns


def main():
    zscores = pd.read_feather(snakemake.input.zscores).set_index("FBgn")
    zscore = zscores.loc[snakemake.wildcards.FBgn, :]
    umap = pd.read_feather(snakemake.input.umap).set_index("cell_id")

    df = umap.join(zscore.rename("zscore"))
    plot(df, snakemake.wildcards.symbol, snakemake.output[0])


def plot(umap, symbol, output_file):
    # set up cmap
    cmap = plt.get_cmap("viridis", 512)
    norm = matplotlib.colors.Normalize(-3, 3)

    # Plot
    ax = sns.scatterplot(
        x="UMAP_1",
        y="UMAP_2",
        data=umap.sort_values("zscore"),
        hue="zscore",
        hue_norm=norm,
        palette=cmap,
        s=8,
        linewidth=0,
        rasterized=True,
        legend=False,
    )
    sns.despine(ax=ax)
    ax.set_title(f"{symbol}", fontstyle="italic")
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.colorbar(sm)

    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="seurat3-cluster-wf",
            input=dict(
                zscores="../output/seurat3-cluster-wf/zscore_by_cell.feather",
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
            ),
            wildcards=dict(symbol="bol", FBgn="FBgn0011206")
        )

    main()
