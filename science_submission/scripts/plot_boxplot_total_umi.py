import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


FNAME = snakemake.input[0]

CLUSTER_COLORS = snakemake.params.cluster_colors

ONAME = snakemake.output[0]

# Debug Settings
# FNAME = 'output/science_submission/raw_by_cluster_rep.feather'
# import yaml
# config = yaml.safe_load(open('config/common.yaml'))
# CLUSTER_COLORS = yaml.full_load(open('config/colors.yaml'))['clusters']


def main():
    df = pd.read_feather(FNAME).groupby(["FBgn", "cluster"]).UMI.sum().to_frame().reset_index()

    plt.style.use("scripts/figure_styles.mplstyle")
    fig, ax = plt.subplots()

    sns.boxplot(
        "cluster",
        "UMI",
        data=df,
        palette=CLUSTER_COLORS,
        ax=ax,
        notch=True,
        linewidth=0.5,
        showfliers=False,
    )
    ax.set(xlabel="", ylabel="UMI Per Cell")

    ax.set_axisbelow(True)
    ax.grid(axis="y")

    fig.savefig(ONAME, bbox_inches="tight")


if __name__ == "__main__":
    main()
