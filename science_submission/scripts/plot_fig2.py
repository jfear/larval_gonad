""" Figure 2: Overall expression by chromosome arm

This is a set of 6 panels. Each panel has:

- x-axis: Chromosome Arm

The top row as:

- y-axis: number reads mapping per chrom / number of reads from each sample /
number expressed genes (> 5) / 1e3

The bottom row as:

- y-axis: percent gene expressed (i.e. > 5) per chrom / number expressed genes per arm * 100

Each panel plots a different set of data:

- A: Adult bulk samples with testis and ovary data (prop_reads)
- B: L3 larval bulk samples with testis and ovary data (prop_reads)
- C: L3 larval single cell testis data split by cell lineage (prop_reads)
- D: Adult bulk samples with testis and ovary data (pct_expressed)
- E: L3 larval bulk samples with testis and ovary data (pct_expressed)
- F: L3 larval single cell testis data split by cell lineage (pct_expressed)

"""
import os
from typing import List
from collections import ChainMap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad import plotting
from fig2_data import NAMES

PLOT_DEFAULTS = dict(x="chrom", order=["X", "2L", "2R", "3L", "3R", "4", "Y"], capsize=0.1)


def main():
    plt.style.use(["2c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width * 0.6)

    _, axes = plt.subplots(2, 3, sharex=True, gridspec_kw=dict(hspace=0.05, wspace=0.05))
    plot_read_counts(axes[0, :])
    plot_pct_expressed_genes(axes[1, :])

    for ax, label in zip(axes.flat, ["A", "B", "C", "D", "E", "F"]):
        ax.text(0.1, 0.99, label, fontweight="bold", transform=ax.transAxes)

    plt.savefig(snakemake.output[0])


def plot_read_counts(axes: List[plt.Axes]):
    kwargs = ChainMap({}, PLOT_DEFAULTS, dict(y=NAMES[0]))
    df = pd.read_feather(snakemake.input.prop_reads)
    plot_adult_bulk(axes[0], df, **kwargs)
    plot_L3_bulk(axes[1], df, **kwargs)
    plot_L3_scLineage(axes[2], df, **kwargs)

    for ax in axes:
        ax.set(ylim=(0, 0.15), yticks=[0, 0.1, 0.15], xlabel="")

    for ax in axes[1:]:
        ax.set(ylabel="")


def plot_pct_expressed_genes(axes: List[plt.Axes]):
    kwargs = ChainMap({}, PLOT_DEFAULTS, dict(y=NAMES[1]))
    df = pd.read_feather(snakemake.input.pct_expressed)
    plot_adult_bulk(axes[0], df, **kwargs)
    plot_L3_bulk(axes[1], df, **kwargs)
    plot_L3_scLineage(axes[2], df, **kwargs)

    for ax in axes:
        ax.set(ylim=(0, 100), yticks=[0, 50, 100], xlabel="Chromosome", title="")

    for ax in axes[1:]:
        ax.set(ylabel="")


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


def plot_L3_scLineage(ax: plt.Axes, df: pd.DataFrame, **kwargs):
    germline = snakemake.params.colors["germline"][0]
    soma = snakemake.params.colors["soma"][0]
    dat = df.query("stage == 'L3' and data_source == 'scRNA-Seq' and cell_type != 'None'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="cell_type",
        hue_order=["Germline", "Somatic"],
        palette=[germline, soma],
        **kwargs
    )
    ax.set(title="L3 scRNA-Seq Lineages", xlabel="")
    ax.set_xlabel("")
    ax.set_title("L3 scRNA-Seq Lineages", y=0.9)
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        COLORS = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                prop_reads="../output/science_submission/fig2_data_prop_reads.feather",
                pct_expressed="../output/science_submission/fig2_data_pct_expressed.feather",
            ),
            params=dict(colors=COLORS),
        )

    main()
