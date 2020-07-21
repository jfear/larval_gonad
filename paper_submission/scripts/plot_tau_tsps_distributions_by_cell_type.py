"""Violin plot showing tau and TSPS values by cell type."""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib

import larval_gonad.plotting  # pylint: disable=unused-import

plt.style.use("minimal")


def main():
    scores = load_scores()
    df = pd.concat(
        [
            create_gene_subset(
                scores, snakemake.input.expressed_fbgns, "All Expressed"
            ),
            create_gene_subset(
                scores, snakemake.input.widely_expressed_fbgns, "Widely Expressed"
            ),
        ],
        ignore_index=True,
    )

    g = sns.FacetGrid(
        row="metric",
        row_order=["tau", "tsps"],
        col="gene_set",
        col_order=["All Expressed", "Widely Expressed"],
        data=df,
        sharey="row",
        sharex=True,
        aspect=1.6,
    )
    g.map(
        sns.violinplot,
        "cluster",
        "score",
        order=snakemake.params.order,
        palette=snakemake.params.colors,
    )
    g.set_titles("{row_name}:{col_name}")
    tweak(g)

    g.savefig(snakemake.output[0])


def load_scores():
    return (
        pd.concat(
            [
                pd.read_feather(
                    snakemake.input.tau, columns=["FBgn", "male_tau"]
                ).set_index("FBgn"),
                pd.read_feather(
                    snakemake.input.tsps, columns=["FBgn", "male_tsps"]
                ).set_index("FBgn"),
            ],
            axis=1,
            sort=False,
        )
        .rename(columns={"male_tau": "tau", "male_tsps": "tsps"})
        .reset_index()
        .melt(id_vars="FBgn", var_name="metric", value_name="score")
        .merge(
            pd.read_feather(snakemake.input.clusters, columns=["FBgn", "cluster"]),
            on="FBgn",
            how="inner",
        )
    )


def create_gene_subset(df: pd.DataFrame, gene_set: str, name: str):
    fbgns = joblib.load(gene_set)  # pylint: disable=unused-variable
    return df.query("FBgn in @fbgns").assign(gene_set=name)


def tweak(g: sns.FacetGrid):
    g.set_xlabels("")

    ax: plt.Axes
    for ax in g.axes.flat:
        row, col = ax.get_title().split(":")
        ax.set_title(col)
        ax.set_ylabel(row)

    # Remove y-label on inner plots
    for ax in g.axes[:, 1]:
        ax.set_ylabel("")

    # Remove title on inner plots
    for ax in g.axes[1, :]:
        ax.set_title("")

    # Clean up xticklabels
    for ax in g.axes[1, :]:
        new_labels = [
            snakemake.params.names[label.get_text()] for label in ax.get_xticklabels()
        ]
        ax.set_xticklabels(new_labels, rotation=45)

    return g


if __name__ == "__main__":
    main()
