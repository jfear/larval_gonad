import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection
from scipy.stats import mannwhitneyu


def main():
    volume = pd.read_csv(snakemake.input[0])
    ax = sns.boxplot(
        "chrom",
        "um3",
        data=volume.melt(var_name="chrom", value_name="um3"),
        palette=snakemake.params.colors,
    )

    # Clean up plot
    ax.set(ylabel=r"$\mu m^{3}$", xlabel="")
    sns.despine(ax=ax)

    # Test that not significant
    pval = np.round(mannwhitneyu(volume["X"], volume["2L"], alternative="less")[1], 3)
    if pval > 0.05:
        # Extend axis and add NS.
        _max = volume.max().max() + 2
        ax.set_ylim(None, _max)

        ax.text(0.5, 0.99, f"NS (p={pval})", transform=ax.transAxes, va="top", ha="center")
        l = plt.Line2D([0.3, 0.7], [0.94, 0.94], transform=ax.transAxes, color="k", lw=0.8, ls="-")
        ax.add_line(l)

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="fish-wf",
            input="../data/external/miriam/oligopaint_volumes.csv",
            params=dict(colors=["red", "grey"]),
        )
    plt.style.use("../config/figure_styles.mplstyle")

    main()
