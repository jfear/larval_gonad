import os

from more_itertools import flatten
import pandas as pd


def main():
    df = pd.concat([read_cell_calls(), read_cluster(), read_umap()], axis=1, sort=False).join(
        read_metadata(), how="inner"
    ).assign(flag_cell_used_in_study=lambda x: (x.is_cell ^ x.scrublet_is_multi) & (x.nFeature <= 5000))  # type: pd.DataFrame

    df.rename_axis("cell_id", inplace=True)
    df.reindex(
        columns=[
            "rep",
            "nUMI",
            "nFeature",
            "cellranger_is_cell",
            "droputils_is_cell",
            "is_cell",
            "scrublet_is_multi",
            "flag_cell_used_in_study",
            "cluster",
            "UMAP_1",
            "UMAP_2",
        ]
    ).to_csv(snakemake.output[0], sep="\t")


def read_metadata():
    return (
        pd.read_feather(snakemake.input.metadata)
        .assign(rep=lambda x: x.cell_id.str.extract(r"(rep\d)_.*", expand=False))
        .loc[:, ["cell_id", "rep", "nUMI", "nFeature"]]
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
            ["cell_id", "cellranger_is_cell", "droputils_is_cell", "is_cell", "scrublet_is_multi"],
        ]
        .set_index("cell_id")
    )


def read_scrublet(fname):
    with open(fname) as fh:
        return fh.read().strip().split("\n")


def read_umap():
    return pd.read_feather(snakemake.input.umap).set_index("cell_id")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            workdir="science_submission",
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
            ),
        )

    main()
