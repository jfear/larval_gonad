import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from larval_gonad import plotting

plt.style.use(["2c", "science_base"])
figsize = (plt.rcParams["figure.figsize"][0], plt.rcParams["figure.figsize"][0])

plot_defaults = dict(
    x="chrom",
    y="Proportion Reads\n(Scaled by Number of Expressed Genes)",
    order=["X", "2L", "2R", "3L", "3R", "4", "Y"],
    capsize=0.1,
)


def main():
    df = munge()

    fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=True, sharey=True)
    plot_L3_bulk(axes[0, 0], df)
    plot_adult_bulk(axes[0, 1], df)
    plot_L3_scAgg(axes[1, 0], df)
    plot_L3_scLineage(axes[1, 1], df)
    axes[0, 0].set(ylim=(0, 0.15), yticks=[0, 0.06, 0.12])

    for ax in axes[:, 1]:
        ax.set(ylabel="")

    for ax in axes[1, :]:
        ax.set(xlabel="Chromosome")

    plt.savefig(snakemake.output[0])


def plot_L3_bulk(ax, df):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'L3' and data_source == 'RNA-Seq'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **plot_defaults
    )
    ax.set(title="L3 Bulk", xlabel="")
    ax.get_legend().remove()
    return ax


def plot_adult_bulk(ax, df):
    testis = snakemake.params.colors["testis"][0]
    ovary = snakemake.params.colors["ovary"][0]
    dat = df.query("stage == 'adult' and data_source == 'RNA-Seq'")
    sns.barplot(
        data=dat,
        ax=ax,
        hue="tissue",
        hue_order=["testis", "ovary"],
        palette=[testis, ovary],
        **plot_defaults
    )
    ax.set(title="Adult Bulk", xlabel="")
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


def plot_L3_scAgg(ax, df):
    testis = snakemake.params.colors["testis"][0]
    dat = df.query("stage == 'L3' and data_source == 'scRNA-Seq' and cell_type == 'None'")
    sns.barplot(data=dat, ax=ax, palette=[testis], **plot_defaults)
    ax.set(title="L3 scRNA-Seq", xlabel="")
    return ax


def plot_L3_scLineage(ax, df):
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
        **plot_defaults
    )
    ax.set(title="L3 scRNA-Seq Lineages", xlabel="")
    # ax.get_legend().remove()
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


def munge():
    counts = pd.concat(
        [
            pd.read_feather(snakemake.input.adult_bulk),
            pd.read_feather(snakemake.input.L3_bulk),
            pd.read_feather(snakemake.input.L3_sc_agg),
            pd.read_feather(snakemake.input.L3_sc_by_class),
        ],
        sort=False,
    ).fillna("None")

    expressed = counts.loc[counts.Count > 0, "FBgn"].unique().tolist()

    fbgn2chrom = get_fbgn2chrom(expressed).value_counts().rename("num_genes")
    fbgn2chrom.index.name = "chrom"

    total_counts = counts.groupby(
        ["stage", "data_source", "cell_type", "rep", "tissue"]
    ).Count.sum()

    prop_counts = (
        counts.groupby(["chrom", "stage", "data_source", "cell_type", "rep", "tissue"])
        .Count.sum()
        .div(total_counts)
        .div(fbgn2chrom / 1e3, axis=0)
        .rename(plot_defaults["y"])
    )

    return prop_counts.to_frame().reset_index()


def get_fbgn2chrom(fbgns):
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
        .reindex(fbgns)
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
