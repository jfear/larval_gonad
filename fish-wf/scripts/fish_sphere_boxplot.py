import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt

from scipy.stats import ttest_rel

from larval_gonad import plotting

def main():
    plt.style.use(["1c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width)

    sphere = pd.read_csv(snakemake.input[0])
    ax = sns.boxplot(
        "chrom",
        "um3",
        data=sphere.melt(var_name="chrom", value_name="um3"),
        palette=snakemake.params.colors,
        notch=True
    )

    # Clean up plot
    ax.set(ylabel=r"$\Psi$", xlabel="")
    sns.despine(ax=ax)

    # Test that not significant
    pval = np.round(ttest_rel(sphere["X"], sphere["2L"])[1], 3)
    if pval <= 0.05:
        # Extend axis and add NS.
        _max = sphere.max().max() + 0.05
        ax.set_ylim(None, _max)

        ax.text(0.5, 0.99, f"p = {pval}", transform=ax.transAxes, va="top", ha="center")
        l = plt.Line2D([0.3, 0.7], [0.94, 0.94], transform=ax.transAxes, color="k", lw=0.8, ls="-")
        ax.add_line(l)

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="fish-wf",
            input="../data/external/miriam/oligopaint_sphere.csv",
            params=dict(colors=["red", "grey"]),
        )

    main()
