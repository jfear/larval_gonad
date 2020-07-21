"""Expression panel of the entire literature gene table."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad.normalization import tpm
import larval_gonad.plotting  # pylint: disable=unused-import

plt.style.use("minimal")

FBGNS = [
    "FBgn0039044",  # p53
    "FBgn0243486",  # rdo
    "FBgn0026573",  # ADD1
    "FBgn0083963",  # Nlg3
    "FBgn0011206",  # bol
    "FBgn0264953",  # Piezo
]

SYMBOLS = ["p53", "rdo", "ADD1", "Nlg3", "bol", "Piezo"]


def main():
    df = (
        pd.read_feather(snakemake.input.raw)
        .pipe(calculate_tpm)
        .query("FBgn in @FBGNS")
        .pipe(add_gene_symbol)
        .pipe(add_cluster_annotation)
        .assign(scaled=lambda data: scale_fbgn_by_max_tpm(data))
    )

    grid_defaults = dict(
        data=df, col="gene_symbol", col_wrap=2, col_order=SYMBOLS, sharey=True
    )

    bar_defaults = dict(
        palette=snakemake.params.colors,
        order=snakemake.params.cluster_order,
        capsize=0.1,
        errwidth=1,
    )

    # Plot Scaled TPMs
    g = sns.FacetGrid(**grid_defaults)
    g.map(sns.barplot, "cluster", "TPM", **bar_defaults)
    tweak_axes(g)
    g.savefig(snakemake.output.tpm)

    # Plot Scaled TPMs
    g = sns.FacetGrid(**grid_defaults)
    g.map(sns.barplot, "cluster", "scaled", **bar_defaults)
    tweak_axes(g)
    g.savefig(snakemake.output.scaled)


def calculate_tpm(df: pd.DataFrame) -> pd.DataFrame:
    gene_lengths = (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "length"])
        .set_index("FBgn")
        .squeeze()
        .reindex(df.FBgn)
    )
    return tpm(df.set_index("FBgn"), gene_lengths).reset_index()


def add_gene_symbol(df: pd.DataFrame) -> pd.DataFrame:
    fbgn2symbol = pd.read_feather(
        snakemake.input.gene_annot, columns=["FBgn", "gene_symbol"]
    )
    return df.merge(fbgn2symbol)


def add_cluster_annotation(df: pd.DataFrame) -> pd.DataFrame:
    cell_id2cluster = pd.read_feather(
        snakemake.input.clusters, columns=["cell_id", "cluster"]
    )
    return df.melt(
        id_vars=["FBgn", "gene_symbol"], var_name="cell_id", value_name="TPM"
    ).merge(cell_id2cluster)


def scale_fbgn_by_max_tpm(df: pd.DataFrame) -> np.ndarray:
    fbgn_max_tpm = df.groupby(["FBgn", "cluster"]).TPM.mean().groupby("FBgn").max()
    return df.set_index(["FBgn", "cell_id"]).TPM.div(fbgn_max_tpm, level="FBgn").values


def tweak_axes(g: sns.FacetGrid):
    g.set_titles("{col_name}", fontstyle="italic")
    g.set_xlabels("")
    for ax in g.axes.ravel():
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    return g


if __name__ == "__main__":
    main()
