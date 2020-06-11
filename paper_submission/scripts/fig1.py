""" Figure 1: Overall expression by chromosome arm

- x-axis: Chromosome Arm
- y-axis: number reads mapping per chrom / number of reads from each sample /
number expressed genes (> 5) / 1e3

- A: Adult bulk samples with testis and ovary data (prop_reads)
- B: L3 larval bulk samples with testis and ovary data (prop_reads)

"""
import os
from typing import List
from collections import ChainMap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad import plotting
from fig1_data import NAMES

PLOT_DEFAULTS = dict(x="chrom", order=["X", "2L", "2R", "3L", "3R", "4", "Y"], capsize=0.1)


def main():
    plt.style.use(["2c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width * 0.6)

    _, (ax1, ax2) = plt.subplots(1, 2, sharex=True, gridspec_kw=dict(hspace=0.05, wspace=0.05))
    plot_read_counts([ax1, ax2])

    for ax, label in zip([ax1, ax2], ["A", "B"]):
        ax.text(0.1, 0.99, label, fontweight="bold", transform=ax.transAxes)

    plt.savefig(snakemake.output[0])


def plot_read_counts(axes: List[plt.Axes]):
    kwargs = ChainMap({}, PLOT_DEFAULTS, dict(y=NAMES[0]))
    df = pd.read_feather(snakemake.input[0])
    plot_adult_bulk(axes[0], df, **kwargs)
    plot_L3_bulk(axes[1], df, **kwargs)

    for ax in axes:
        ax.set(ylim=(0, 110), yticks=[0, 50, 100], xlabel="")

    axes[1].set(ylabel="")


def plot_L3_bulk(ax: plt.Axes, df: pd.DataFrame, **kwargs):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'L3' and data_source == 'RNA-Seq'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **kwargs
    )
    ax.set_xlabel("")
    ax.set_title("L3 Bulk", y=0.9)
    ax.get_legend().remove()
    return ax


def plot_adult_bulk(ax: plt.Axes, df: pd.DataFrame, **kwargs):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'adult' and data_source == 'RNA-Seq'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **kwargs
    )
    ax.set_xlabel("")
    ax.set_title("Adult Bulk", y=0.9)
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        COLORS = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="paper_submission",
            input="../output/paper_submission/fig2_data_prop_reads.feather",
            params=dict(colors=COLORS),
        )

    main()
