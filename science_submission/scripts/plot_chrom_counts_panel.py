import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plot_defaults = dict(x="chrom", y="Proportion Reads", order=["X", "2L", "2R", "3L", "3R", "4", "Y"])


def main():
    df = munge()

    fig, axes = plt.subplots(3, 2, figsize=plt.figaspect(1.5), sharex=True, sharey=True)
    plot_L3_bulk(axes[0, 0], df)
    plot_adult_bulk(axes[0, 1], df)
    plot_L3_scAgg(axes[1, 0], df)
    # plot_adult_scAgg(axes[1, 1], df)
    plot_L3_scLineage(axes[2, 0], df)
    # plot_adult_scLineage(axes[2, 1], df)


def plot_L3_bulk(ax, df):
    dat = df.query("stage == 'L3' and data_source == 'RNA-Seq'")
    sns.barplot(data=dat, ax=ax, hue="tissue", palette=["C3", "C0"], **plot_defaults)
    ax.set(title="L3 Bulk", xlabel="")
    ax.get_legend().remove()
    return ax


def plot_adult_bulk(ax, df):
    dat = df.query("stage == 'adult' and data_source == 'RNA-Seq'")
    sns.barplot(data=dat, ax=ax, hue="tissue", palette=["C3", "C0"], **plot_defaults)
    ax.set(title="Adult Bulk", xlabel="")
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax


def plot_L3_scAgg(ax, df):
    dat = df.query("stage == 'L3' and data_source == 'scRNA-Seq' and cell_type == 'None'")
    sns.barplot(data=dat, ax=ax, palette=["C0"], **plot_defaults)
    ax.set(title="L3 scRNA-Seq", xlabel="")
    return ax


def plot_adult_scAgg(ax, df):
    dat = df.query("stage == 'adult' and data_source == 'scRNA-Seq' and cell_type == 'None'")
    sns.barplot(data=dat, ax=ax, palette=["C0"], **plot_defaults)
    ax.set(title="Adult scRNA-Seq", xlabel="")
    return ax


def plot_L3_scLineage(ax, df):
    dat = df.query("stage == 'L3' and data_source == 'scRNA-Seq' and cell_type != 'None'")
    sns.barplot(data=dat, ax=ax, hue="cell_type", hue_order=["Germline", "Cyst Lineage", "Other Somatic"], **plot_defaults)
    ax.set(title="L3 scRNA-Seq Lineages", xlabel="")
    ax.get_legend().remove()
    return ax


def plot_adult_scLineage(ax, df):
    dat = df.query("stage == 'adult' and data_source == 'scRNA-Seq' and cell_type != 'None'")
    sns.barplot(data=dat, ax=ax, hue="cell_type", hue_order=["Germline", "Cyst Lineage", "Other Somatic"], **plot_defaults)
    ax.set(title="Adult scRNA-Seq Lineages", xlabel="")
    ax.legend(loc="upper left", bbox_to_anchor=[1, 1])
    return ax



def munge():
    fbgn2chrom = get_fbgn2chrom().value_counts()

    counts = pd.concat(
        [
            pd.read_feather(snakemake.input.adult_bulk),
            pd.read_feather(snakemake.input.L3_bulk),
            pd.read_feather(snakemake.input.L3_sc_agg),
            pd.read_feather(snakemake.input.L3_sc_by_class),
            # pd.read_feather(snakemake.input.adult_sc_agg),
            # pd.read_feather(snakemake.input.adult_sc_by_class),
        ],
        sort=False,
    ).fillna("None")

    total_counts = counts.groupby(
        ["stage", "data_source", "cell_type", "rep", "tissue"]
    ).Count.sum()

    prop_counts = (
        counts.groupby(["chrom", "stage", "data_source", "cell_type", "rep", "tissue"])
        .Count.sum()
        .div(total_counts)
        .rename("Proportion Reads")
    )
    return prop_counts.to_frame().reset_index()


def get_fbgn2chrom():
    return (
        pd.read_feather(snakemake.input.gene_annot, columns=["FBgn", "FB_chrom"])
        .set_index("FBgn")
        .query('FB_chrom in ["X", "2L", "2R", "3L", "3R", "4", "Y"]')
        .squeeze()
        .rename("chrom")
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
            input=dict(
                gene_annot="../references/gene_annotation_dmel_r6-26.feather",
                adult_bulk="../output/expression-atlas-wf/w1118_gene_counts.feather",
                L3_bulk="../output/bulk2-rnaseq-wf/testis_ovary_counts.feather",
                L3_sc_agg="../output/seurat3-cluster-wf/aggegated_gene_counts.feather",
                L3_sc_by_class="../output/seurat3-cluster-wf/aggegated_gene_counts_by_germ_soma.feather",
                adult_sc_agg="",
                adult_sc_by_class="",
            ),
        )

    main()
