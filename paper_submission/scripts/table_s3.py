import os
import pandas as pd

from larval_gonad.constants import L3_BULK, ADULT_BULK


def main():
    with pd.ExcelWriter(snakemake.output[0]) as workbook:
        add_counts(workbook)
        add_corr1(workbook)
        add_corr2(workbook)


def add_counts(workbook: pd.ExcelWriter):
    counts = pd.concat([load_larval_bulk(), load_adult_bulk()], axis=1, join="inner", sort=False)
    counts.to_excel(workbook, sheet_name="Gene Level")


def load_larval_bulk():
    return (
        pd.read_csv(snakemake.input.larval, sep="\t")
        .pipe(lambda x: x[~x.FBgn.str.contains("ERCC")])
        .rename(columns=L3_BULK)
        .sort_index(axis=1)
        .set_index(["FBgn", "gene_symbol", "Chromosome"])
    )


def load_adult_bulk():
    return (
        pd.read_csv(snakemake.input.adult, sep="\t")
        .rename(columns=ADULT_BULK)
        .sort_index(axis=1)
        .set_index(["FBgn", "gene_symbol", "Chromosome"])
    )


def add_corr1(workbook: pd.ExcelWriter):
    corr1 = pd.read_csv(snakemake.input.corr_reps, sep="\t", index_col=0)
    corr1.to_excel(workbook, sheet_name="Spearman Among Samples")


def add_corr2(workbook: pd.ExcelWriter):
    corr2 = pd.read_csv(snakemake.input.corr_clusters, sep="\t", index_col=0)
    corr2.to_excel(workbook, sheet_name="Spearman With Clusters")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from larval_gonad.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                larval="../../output/paper_submission/larval_bulk_RNA-Seq_gene_level_counts.tsv.gz",
                adult="../../output/paper_submission/adult_bulk_RNA-Seq_gene_level_counts.tsv.gz",
                corr_reps="../../output/paper_submission/spearman_corr_bulk_vs_sum_single_cell.tsv",
                corr_clusters="../../output/paper_submission/spearman_corr_bulk_vs_single_cell_clusters.tsv",
            ),
            output="../../output/paper_submission/TableS3.xlsx",
        )

    main()
