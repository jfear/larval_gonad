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
    species = snakemake.wildcards.species
    tissue = snakemake.wildcards.tissue

    df, male, female = get_data(snakemake.input[0])

    fig, ax = plt.subplots(figsize=plt.figaspect(2))
    plot(df, ax, f"{species} ({tissue})", snakemake.params.colors)
    add_pval(male, ax)
    add_pval(female, ax)

    plt.savefig(snakemake.output[0])


def get_data(file_name: str):
    """Loads shelve and get data ready for plotting"""
    db = shelve_load(file_name)
    df = db["data"]
    male = db["male_qval"].fillna(1)  # If p-vals are NaN set to 1
    female = db["female_qval"].fillna(1)

    # Change female Y location to come from top.
    female.y = 1 - female.y

    return df, male, female


def plot(df, ax, title, colors):
    """Make main plot"""
    df.plot.bar(
        stacked=True, ax=ax, color=colors, width=0.9, edgecolor="k", linewidth=0.2
    )
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.set(title=title, xlabel="", ylim=(0, 1))
    ax.legend_ = None
    sns.despine(ax=ax, left=True, bottom=True)
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
            input='../output/expression-atlas-wf/sex_bias_by_muller_x_and_4th/dana_AC.dat',
            params=dict(
                colors=["blue", "lightgray", "red"],
            ),
            wildcards=dict(species="dana", tissue="AC")
        )
        plt.style.use("../config/figure_styles.mplstyle")

    main()
