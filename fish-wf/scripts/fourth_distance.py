import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection
from scipy.stats import wilcoxon

from larval_gonad import plotting

def main():
    plt.style.use(["1c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width)

    distance = (
        pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
        .query("Phase == 'middle/late'")
        .iloc[:, 6:8]
    )
    distance.columns = ["X", "4"]
    distance_melted = distance.melt(var_name="Chromosome", value_name="Distance (μm)")

    ax = sns.boxplot(
        "Chromosome",
        "Distance (μm)",
        data=distance_melted,
        palette=snakemake.params.colors,
        notch=True
    )

    # Clean up plot
    ax.set(xlabel="")
    sns.despine(ax=ax)

    # Test that not significant
    _, pval = wilcoxon(distance.X, distance["4"], alternative="two-sided", correction=True)
    if pval <= 0.05:
        # Extend axis and add NS.
        _max = distance.max().max() + 2
        ax.set_ylim(None, _max)

        ax.text(0.5, 0.99, f"p = {pval:.2E}", transform=ax.transAxes, va="top", ha="center")
        l = plt.Line2D([0.3, 0.7], [0.94, 0.94], transform=ax.transAxes, color="k", lw=0.8, ls="-")
        ax.add_line(l)

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../data/external/camila/fourth_distance_20200228.tsv",
            params=dict(colors=["red", "cyan"]),
        )

    main()
