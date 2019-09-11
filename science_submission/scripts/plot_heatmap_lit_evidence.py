import os
import matplotlib

matplotlib.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


def main():
    fbgn2symbol = (
        pd.read_feather(snakemake.input.gene_metadata, columns=["FBgn", "gene_symbol"])
        .set_index("FBgn")
        .to_dict()["gene_symbol"]
    )
    evidence = pd.read_csv(snakemake.input.lit_evidence, sep="\t", index_col=0).rename(fbgn2symbol)
    this_study = evidence.query('References == "This study"').index.tolist()

    colors = sns.color_palette(snakemake.params.cmap, n_colors=100)

    fig, ax = plt.subplots()
    sns.heatmap(
        evidence.iloc[:, :-1],
        cmap=["lightgray", colors[0], colors[-1]],
        yticklabels=True,
        xticklabels=True,
        cbar=False,
        square=True,
        ax=ax,
    )

    # Add legend
    low = mpatches.Patch(color=colors[0], label="Not Expressed")
    high = mpatches.Patch(color=colors[-1], label="Expressed")
    none = mpatches.Patch(color="lightgray", label="Not Tested")
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1], handles=[low, high, none])

    # Clean up X axis
    ax.set_xlabel("")
    ax.xaxis.set_ticks_position("top")
    plt.setp(ax.get_xticklabels(), rotation=90)

    # Clean up Y axis
    ax.set_ylabel("")
    labels = []
    for l in ax.get_yticklabels():
        l.set(fontsize=6, fontstyle="italic")
        if l.get_text() in this_study:
            l.set(fontsize=6, fontweight="bold")
        labels.append(l)
    ax.set_yticklabels(labels)

    fig.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        config = read_config("config/common.yaml")
        TAG = config["tag"]

        color_config = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_metadata=f"../references/gene_annotation_dmel_{TAG}.feather",
                lit_evidence="../data/external/miriam/lit_gene_dummy_vars.tsv",
            ),
            params=dict(cmap=color_config["heatmap"]),
        )

    plt.style.use("../config/figure_styles.mplstyle")
    plt.rcParams.update({
        "figure.figsize": (4, 6)
    })

    main()
