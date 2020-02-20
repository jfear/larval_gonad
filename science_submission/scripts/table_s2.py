import os
import pandas as pd


def main():
    with pd.ExcelWriter(snakemake.output[0]) as workbook:
        add_cell_level(workbook)
        add_gene_level(workbook)
        add_cluster_level(workbook)
        add_biomarkers(workbook)


def add_cell_level(workbook: pd.ExcelWriter):
    df = pd.read_csv(snakemake.input.cell, sep="\t", index_col=0)
    df.to_excel(workbook, sheet_name="Cell Level Data")


def add_gene_level(workbook):
    df = pd.read_csv(snakemake.input.gene, sep="\t", index_col=[0, 1, 2])
    df.to_excel(workbook, sheet_name="Gene Level Data")


def add_cluster_level(workbook):
    df = pd.read_csv(snakemake.input.cluster, sep="\t", index_col=0)
    df.to_excel(workbook, sheet_name="Cluster Level Data")


def add_biomarkers(workbook):
    df = pd.read_csv(snakemake.input.biomarkers, sep="\t", index_col=[0, 1, 2])
    df.to_excel(workbook, sheet_name="One vs Rest (Biomarkers)")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                cell="../../output/science_submission/cell_metadata.tsv",
                gene="../../output/science_submission/gene_metadata.tsv",
                cluster="../../output/science_submission/cluster_metadata.tsv",
                biomarkers="../../output/science_submission/one_vs_rest.tsv",
            ),
            output="../../output/science_submission/TableS2.xlsx",
        )

    main()
