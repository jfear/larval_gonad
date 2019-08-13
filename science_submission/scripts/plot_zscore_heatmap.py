from itertools import chain

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

from larval_gonad.io import shelve_load


def main(snake):
    db = shelve_load(snake["input_file"], "data", "blocks", "annotation")
    yticks = db["data"].shape[0] < 100

    defaults = dict(
        xticklabels=True,
        yticklabels=False,
        rasterized=True,
        cmap=snake["cmap"],
        vmin=-3,
        vmax=3,
        cbar_kws=dict(label="Z-Score (TPM)", ticks=[-3, 0, 3], orientation="horizontal"),
    )

    if yticks:
        defaults.update(
            dict(yticklabels=True, linewidths=0.009, linecolor="gray", rasterized=False)
        )

    if db['data'].shape[0] < 20:
        figsize=(4, 4)
    elif db['data'].shape[0] < 50:
        figsize=(4, 6)
    else:
        figsize=(4, 10)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 1, height_ratios=[1, 0.01], hspace=0.01)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    sns.heatmap(db["data"], ax=ax, cbar_ax=cax, **defaults)
    ax.set(xlabel="", ylabel="")

    # Clean up X axis
    ax.xaxis.set_ticks_position("top")
    ax.set_xticklabels(
        list(chain.from_iterable([("", x, "") for x in db["data"].columns.levels[0]])),
        ha="center",
        va="bottom",
    )

    # Clean up Y axis
    if yticks:
        plt.setp(ax.get_yticklabels(), fontstyle="italic", fontsize=8)

    # Add lines separating cell types
    for i in range(1, 10):
        ax.axvline(i * 3, color="w")

    # Add lines separating lit genes
    if db["blocks"] is not None:
        for loc in db["blocks"][:-1]:
            ax.axhline(loc, color="w")

    fig.savefig(snake["output_file"])


if __name__ == "__main__":
    SNAKE = dict(
        input_file=snakemake.input[0], output_file=snakemake.output[0], cmap=snakemake.params[0]
    )

    plt.style.use("scripts/figure_styles.mplstyle")

    # Debug Settings
    # import os
    # try:
    #     os.chdir(os.path.join(os.getcwd(), "science_submission/scripts"))
    #     print(os.getcwd())
    # except:
    #     pass
    # from larval_gonad.config import read_config
    # snake = dict(
    #     input_file="../../output/science_submission/db/zscore_expressed_genes.dat",
    #     output_file="",
    #     cmap=read_config("../../config/colors.yaml")["heatmap"],
    # )
    # plt.style.use("figure_styles.mplstyle")

    main(SNAKE)
