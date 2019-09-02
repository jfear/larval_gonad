"""Use tsps scores to define housekeeping genes.

1. Uses tsps scores to define housekeeping genes.
"""
import os
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from larval_gonad.io import pickle_dump


def main():
    tsps = pd.read_feather(snakemake.input[0]).set_index("YOgn").dropna(how="all")

    # Plot distribution
    ax = tsps.plot.kde()
    ax.axvline(snakemake.params[0], color="k", ls="--")
    ax.set(xlim=(-2, 2))
    plt.savefig(snakemake.output.svg)

    # Housekeeping genes
    flag_housekeeping = tsps.fillna(999) <= snakemake.params[0]

    # Save list of housekeeping genes
    pickle_dump(flag_housekeeping[flag_housekeeping.male].index.tolist(), snakemake.output.male)
    pickle_dump(flag_housekeeping[flag_housekeeping.female].index.tolist(), snakemake.output.female)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="expression-atlas-wf",
            input="../output/expression-atlas-wf/tsps/w1118.feather",
            params=1,
        )

    main()
