import os

from more_itertools import flatten
import pandas as pd

from larval_gonad.constants import L3_SC


def main():
    df = (
        pd.concat(
            [
                read_cell_calls(),
                read_cluster(),
                read_umap(),
                read_ratios_all_genes(),
                read_ratios_common_genes(),
                read_gsea(),
            ],
            axis=1,
            sort=False,
        )
        .join(read_metadata(), how="inner")
        .assign(
            is_cell_used_in_study=lambda x: (x.is_cell ^ x.scrublet_is_multi)
            & (x.nFeature <= 5000)
        )
    )  # type: pd.DataFrame

    df.rename_axis("cell_id", inplace=True)
    df.reindex(
        columns=[
            "sample_id",
            "rep",
            "nUMI",
            "nFeature",
            "cellranger_is_cell",
            "droputils_is_cell",
            "is_cell",
            "scrublet_is_multi",
            "is_cell_used_in_study",
            "cluster",
            "UMAP_1",
            "UMAP_2",
            "male_bias_enrichment_score",
            "all_genes_x_to_a_ratio",
            "all_genes_fourth_to_a_ratio",
            "all_genes_y_to_a_ratio",
            "common_genes_x_to_a_ratio",
            "common_genes_fourth_to_a_ratio",
            "common_genes_y_to_a_ratio",
        ]
    ).to_csv(snakemake.output[0], sep="\t")


def read_metadata():
    return (
        pd.read_feather(snakemake.input.metadata)
        .assign(rep=lambda x: x.cell_id.str.extract(r"(rep\d)_.*", expand=False))
        .assign(sample_id=lambda x: x.rep.map(L3_SC))
        .loc[:, ["cell_id", "sample_id", "rep", "nUMI", "nFeature"]]
        .set_index("cell_id")
    )


def read_cluster():
    return pd.read_feather(snakemake.input.clusters).set_index("cell_id").cluster


def read_cell_calls():
    cells = pd.concat(list(map(read_cell_call, snakemake.input.cell_calls)))
    doublets = list(flatten(map(read_scrublet, snakemake.input.scrublets)))
    cells.loc[cells.index.isin(doublets), "scrublet_is_multi"] = True
    return cells


def read_cell_call(fname):
    return (
        pd.read_feather(fname)
        .assign(cellranger_is_cell=lambda x: x["cellranger3-wf"].astype(bool))
        .assign(droputils_is_cell=lambda x: x["droputils"].astype(bool))
        .assign(scrublet_is_multi=False)
        .loc[
            :,
            [
                "cell_id",
                "cellranger_is_cell",
                "droputils_is_cell",
                "is_cell",
                "scrublet_is_multi",
            ],
        ]
        .set_index("cell_id")
    )


def read_scrublet(fname):
    with open(fname) as fh:
        return fh.read().strip().split("\n")


def read_umap():
    return pd.read_feather(snakemake.input.umap).set_index("cell_id")


def read_ratios_all_genes():
    df = (
        pd.read_feather(snakemake.input.autosome_ratios_all_genes)
        .set_index("cell_id")
        .drop(["rep", "cluster"], axis=1)
    )
    cols = [f"all_genes_{x}" for x in df.columns]
    df.columns = cols
    return df


def read_ratios_common_genes():
    df = (
        pd.read_feather(snakemake.input.autosome_ratios_commonly_expressed_genes)
        .set_index("cell_id")
        .drop(["rep", "cluster"], axis=1)
    )
    cols = [f"common_genes_{x}" for x in df.columns]
    df.columns = cols
    return df


def read_gsea():
    return (
        pd.read_feather(snakemake.input.gsea)
        .rename(columns={"index": "cell_id", "male": "male_bias_enrichment_score"})
        .set_index("cell_id")
        .drop("cluster", axis=1)
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="paper_submission",
            input=dict(
                metadata="../output/cellselection-wf/cell_metadata.feather",
                clusters="../output/seurat3-cluster-wf/combined_n3_clusters.feather",
                cell_calls=[
                    "../output/cellselection-wf/testis1_combined_cell_calls.feather",
                    "../output/cellselection-wf/testis2_combined_cell_calls.feather",
                    "../output/cellselection-wf/testis3_combined_cell_calls.feather",
                ],
                scrublets=[
                    "../output/cellselection-wf/testis1_scrublet_dublets.txt",
                    "../output/cellselection-wf/testis2_scrublet_dublets.txt",
                    "../output/cellselection-wf/testis3_scrublet_dublets.txt",
                ],
                umap="../output/seurat3-cluster-wf/combined_n3_umap.feather",
                autosome_ratios_all_genes="../output/x-to-a-wf/autosome_ratios_expressed_by_cell.feather",
                autosome_ratios_commonly_expressed_genes="../output/x-to-a-wf/autosome_ratios_commonly_expressed_by_cell.feather",
            ),
        )

    main()
