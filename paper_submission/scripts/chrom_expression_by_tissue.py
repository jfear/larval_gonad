""" Figure 1: Overall expression by chromosome arm

- x-axis: Chromosome Arm
- y-axis: Avg TPM Per Chromosome

"""
import os
from typing import List
from collections import ChainMap
from string import ascii_uppercase

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad import plotting

TISSUE_MAPPER = {
    "WB": "Adult Whole Body",
    "HD": "Adult Head",
    "TX": "Adult Thorax",
    "AC": "Adult Abdomen",
    "DG": "Adult Digestive System",
    "RE": "Adult Reproductive System",
    "GE": "Adult Genitalia",
    "GO": "Adult Gonad",
    "L3_GO": "Larval Gonad",
}
TISSUE_ORDER = [v for _, v in TISSUE_MAPPER.items()]

COLORS = [snakemake.params.colors["testis"][0], snakemake.params.colors["ovary"][0]]


def main():
    plt.style.use(["2c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width * 0.6)

    df = pd.read_feather(snakemake.input[0])
    df.tissue = df.tissue.map(TISSUE_MAPPER)

    g = sns.catplot(
        data=df,
        x="chrom",
        y="Avg TPM Per Chromosome",
        order=["X", "Y", "2L", "2R", "3L", "3R", "4"],
        hue="sex",
        hue_order=["m", "f"],
        palette=COLORS,
        col="tissue",
        col_order=TISSUE_ORDER,
        row="strain",
        row_order=["w1118", "OreR"],
        kind="bar",
        height=1.5,
        aspect=1.2,
        sharey=False,
        capsize=0.1,
        errwidth=1,
    )

    g.set_titles("{row_name}\n{col_name}")
    g.set_xlabels("Scaffold")

    labels = ascii_uppercase
    for ax, label in zip(g.axes.flat, labels):
        ax.text(0.01, 0.95, label, fontweight="bold", transform=ax.transAxes)

    plt.subplots_adjust(hspace=0.2, wspace=0.1)
    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    main()
