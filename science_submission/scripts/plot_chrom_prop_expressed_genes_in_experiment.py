""" Percent of genes expressed per chromosome arm.

The FlyBase annotations have a little over 100 genes in the 4th and Y. I
think this number is too large. Here I count only genes that were expressed
in any of these experiments.

This is a set of 4 panels. Each panel has:

- x-axis: Chromosome Arm
- y-axis: percent gene expressed (i.e. > 5) per chrom / number expressed genes per arm * 100

Each panel plots a different set of data:

- A: L3 larval bulk samples with testis and ovary data
- B: Adult bulk samples with testis and ovary data
- C: L3 larval single cell testis data aggregated together
- D: L3 larval single cell testis data split by cell lineage

"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad import plotting

PLOT_DEFAULTS = dict(
    x="chrom",
    y="% Genes Expressed\n(All Genes)",
    order=["X", "2L", "2R", "3L", "3R", "4", "Y"],
    capsize=0.1,
)


def main():
    df = munge()

    plt.style.use(["2c", "science_base"])
    width = plt.rcParams["figure.figsize"][0]
    plt.rcParams["figure.figsize"] = (width, width * 0.8)

    fig, axes = plt.subplots(
        2, 2, sharex=True, sharey=True, gridspec_kw=dict(hspace=0.05, wspace=0.05)
    )
    plot_L3_bulk(axes[0, 0], df)
    plot_adult_bulk(axes[0, 1], df)
    plot_L3_scAgg(axes[1, 0], df)
    plot_L3_scLineage(axes[1, 1], df)
    axes[0, 0].set(ylim=(0, 100), yticks=[0, 50, 100])

    for ax in axes[:, 1]:
        ax.set(ylabel="")

    for ax in axes[1, :]:
        ax.set(xlabel="Chromosome")

    plt.savefig(snakemake.output[0])


def plot_L3_bulk(ax: plt.Axes, df: pd.DataFrame):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'L3' and data_source == 'RNA-Seq'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **PLOT_DEFAULTS
    )
    ax.set_xlabel("")
    ax.set_title("L3 Bulk", y=0.9)
    ax.get_legend().remove()
    return ax


def plot_adult_bulk(ax: plt.Axes, df: pd.DataFrame):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'adult' and data_source == 'RNA-Seq'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **PLOT_DEFAULTS
    )
    ax.set_xlabel("")
    ax.set_title("Adult Bulk", y=0.9)
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


def plot_L3_scAgg(ax: plt.Axes, df: pd.DataFrame):
    testis = snakemake.params.colors["testis"][0]
    dat = df.query("stage == 'L3' and data_source == 'scRNA-Seq' and cell_type == 'None'")
    sns.barplot(data=dat, ax=ax, palette=[testis], **PLOT_DEFAULTS)
    ax.set_xlabel("")
    ax.set_title("L3 scRNA-Seq", y=0.9)
    return ax


def plot_L3_scLineage(ax: plt.Axes, df: pd.DataFrame):
    germline = snakemake.params.colors["germline"][0]
    cyst = snakemake.params.colors["cyst"][0]
    soma = snakemake.params.colors["soma"][0]
    dat = df.query("stage == 'L3' and data_source == 'scRNA-Seq' and cell_type != 'None'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="cell_type",
        hue_order=["Germline", "Cyst Lineage", "Other Somatic"],
        palette=[germline, cyst, soma],
        **PLOT_DEFAULTS
    )
    ax.set(title="L3 scRNA-Seq Lineages", xlabel="")
    ax.set_xlabel("")
    ax.set_title("L3 scRNA-Seq Lineages", y=0.9)
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


def munge():
    counts = (
        pd.concat(
            [
                pd.read_feather(snakemake.input.adult_bulk),
                pd.read_feather(snakemake.input.L3_bulk),
                pd.read_feather(snakemake.input.L3_sc_agg),
                pd.read_feather(snakemake.input.L3_sc_by_class),
            ],
            sort=False,
        )
        .fillna("None")
        .assign(flag_on=lambda x: x.Count > 5)
    )

    expressed_FBgns = counts[counts.flag_on].FBgn.unique()

    fbgn2chrom = get_fbgn2chrom(expressed_FBgns).value_counts().rename("num_genes")
    fbgn2chrom.index.name = "chrom"

    pct_expressed = (
        counts.groupby(["chrom", "stage", "data_source", "cell_type", "rep", "tissue"])
        .flag_on.sum()
        .div(fbgn2chrom, axis=0)
        .mul(100)
        .rename(PLOT_DEFAULTS["y"])
    )

    return pct_expressed.to_frame().reset_index()


def get_fbgn2chrom(expressed):
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
        .reindex(expressed)
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug
        from larval_gonad.config import read_config

        COLORS = read_config("config/colors.yaml")

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                adult_bulk="../output/expression-atlas-wf/w1118_gene_counts.feather",
                L3_bulk="../output/bulk2-rnaseq-wf/testis_ovary_counts.feather",
                L3_sc_agg="../output/seurat3-cluster-wf/aggegated_gene_counts.feather",
                L3_sc_by_class="../output/seurat3-cluster-wf/aggegated_gene_counts_by_germ_soma.feather",
            ),
            params=dict(colors=COLORS),
        )

    main()
