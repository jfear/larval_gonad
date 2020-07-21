import os

import matplotlib
import matplotlib.pyplot as plt

import larval_gonad.plotting

plt.style.use("minimal")


def main():
    fig, ax = plt.subplots()

    for color, name in zip(
        snakemake.params.cluster_colors, snakemake.params.legend_names
    ):
        ax.scatter([], [], color=color, label=name)

    plt.legend(loc="center", ncol=1, frameon=False)
    ax.set_axis_off()

    fig.savefig(snakemake.output[0])


if __name__ == "__main__":
    main()
