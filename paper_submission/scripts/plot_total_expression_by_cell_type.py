"""Plot total expression (TPM) for different gene sets."""
import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import larval_gonad.plotting  # pylint: disable=unused-import


plt.style.use("minimal")


def main():
    tpm = load_gene_expression()

    df = pd.concat(
        [
            create_gene_subset(tpm, snakemake.input.expressed_fbgns, "All Expressed"),
            create_gene_subset(
                tpm, snakemake.input.widely_expressed_fbgns, "Widely Expressed"
            ),
            create_gene_subset(tpm, snakemake.input.tau_fbgns, "Tau"),
            create_gene_subset(tpm, snakemake.input.tsps_fbgns, "TSPS"),
        ],
        ignore_index=True,
    )

    g = sns.FacetGrid(
        df,
        col="gene_set",
        col_order=["All Expressed", "Tau", "TSPS", "Widely Expressed"],
    )
    g.map(
        sns.boxplot,
        "cluster",
        "logTPM",
        order=snakemake.params.order,
        palette=snakemake.params.colors,
        notch=True,
        linewidth=0.5,
        showfliers=False,
    )
    tweak_axes(g)

    g.savefig(snakemake.output[0])


def load_gene_expression() -> pd.DataFrame:
    """Returns log10(TPM) with cluster and chromosome annotations."""
    gene_annot = pd.read_feather(
        snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"]
    ).rename(columns={"FB_chrom": "chrom"})

    return (
        pd.read_feather(snakemake.input.tpm)
        .merge(gene_annot, on="FBgn")
        .assign(logTPM=lambda df: np.log10(df.TPM + 0.01))
    )


def create_gene_subset(df: pd.DataFrame, gene_set: str, name: str):
    fbgns = joblib.load(gene_set)
    return df.query("FBgn in @fbgns").assign(gene_set=name)


def tweak_axes(g: sns.FacetGrid):
    g.set_titles("{col_name}")
    g.set_xlabels("")
    g.set_ylabels(r"$Log_{10}(TPM)$")

    ax: plt.Axes
    for ax in g.axes.flat:
        new_labels = [
            snakemake.params.names[label.get_text()] for label in ax.get_xticklabels()
        ]
        ax.set_xticklabels(new_labels, rotation=45)

    for ax in g.axes.flat[1:]:
        sns.despine(ax=ax, left=True)


if __name__ == "__main__":
    main()
