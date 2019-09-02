"""Demasculinization plot of whole body for all species."""
import os
import pandas as pd

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.io import shelve_load
from larval_gonad.plotting import format_pval


def main():
    fig, axes = plt.subplots(1, len(snakemake.params.species), sharey=True, figsize=(8, 2))
    for species, ax in zip(snakemake.params.species, axes.flat):
        df, male, female = get_data(species)
        plot(df, species, ax)
        add_pval(male, ax)
        add_pval(female, ax)

    for ax in axes.flat[1:]:
        ax.yaxis.set_visible(False)

    plt.savefig(snakemake.output[0])


def get_data(species):
    """Loads shelve and get data ready for plotting"""
    db = shelve_load(snakemake.params.pattern.format(species=species))
    df = db["data"]
    male = db["male_qval"].fillna(1)  # If p-vals are NaN set to 1
    female = db["female_qval"].fillna(1)

    # Change female Y location to come from top.
    female.y = 1 - female.y

    return df, male, female


def plot(df, species, ax):
    """Make main plot"""
    df.plot.bar(
        stacked=True, ax=ax, color=snakemake.params.colors, width=0.9, edgecolor="k", linewidth=0.2
    )
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.set(title=species, xlabel="", ylim=(0, 1))
    ax.legend_ = None
    sns.despine(ax=ax, left=True)
    return ax


def add_pval(df, ax):
    """Add *** to plot"""
    for _, row in df.iterrows():
        format_pval(ax, row.x, row.y, row["q-value"], va="center", fontsize=10)
    return ax


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=[],
            params=dict(
                species=["w1118", "orgR", "dana", "dwil", "dyak"],
                pattern="../output/expression-atlas-wf/sex_bias_by_muller/{species}_WB.dat",
                colors=["blue", "lightgray", "red"],
            ),
        )

    main()
