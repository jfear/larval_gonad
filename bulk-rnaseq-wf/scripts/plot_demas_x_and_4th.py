"""Demasculinization plot of whole body for all species."""
import os
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting import demasculinization, add_pvals


def main():
    db = shelve_load(snakemake.input[0])
    df = db["data"]
    male = db["male_qval"].fillna(1)  # If p-vals are NaN set to 1
    # Change female Y location to come from top.
    female = db["female_qval"].fillna(1)
    female.y = 1 - female.y

    fig, ax = plt.subplots(figsize=plt.figaspect(2))
    demasculinization(df, ax=ax, title=f"larval (GO)", color=snakemake.params.colors)
    add_pvals(male.x.values, male.y.values, male['q-value'].values, ax)
    add_pvals(female.x.values, female.y.values, female['q-value'].values, ax)
    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="bulk-rnaseq-wf",
            input="../output/bulk-rnaseq-wf/testis_bias_by_muller_x_and_4th.dat",
            params=dict(colors=["blue", "lightgray", "red"]),
        )
        plt.style.use("../config/figure_styles.mplstyle")

    main()
